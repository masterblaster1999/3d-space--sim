#include "stellar/sim/MissionLogic.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Cargo.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/Reputation.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace stellar::sim {

struct CargoMissionSpec {
  econ::CommodityId commodity{econ::CommodityId::Food};
  int units{0};
  econ::MarketQuote originQ{};
  econ::MarketQuote destQ{};
};

// Pick a commodity + unit count for a cargo-style mission.
//
// Goals:
//  - Prefer goods that are in higher demand at the destination (dest mid > origin mid).
//  - Ensure the origin station has enough inventory for the chosen unit count.
//  - Avoid generating missions that would be contradictory to contraband rules.
static bool pickCargoMission(core::SplitMix64& rng,
                             core::u64 universeSeed,
                             const econ::StationEconomyState& originState,
                             const econ::StationEconomyModel& originModel,
                             const econ::StationEconomyState& destState,
                             const econ::StationEconomyModel& destModel,
                             core::u32 destFactionId,
                             bool requireLegalAtDest,
                             bool cargoProvided,
                             double maxCargoKg,
                             int minUnits,
                             int maxUnitsHard,
                             CargoMissionSpec& out) {
  struct Cand {
    econ::CommodityId id{econ::CommodityId::Food};
    double score{0.0};
    econ::MarketQuote o{};
    econ::MarketQuote d{};
    int maxUnits{0};
  };

  minUnits = std::max(1, minUnits);
  maxUnitsHard = std::max(minUnits, maxUnitsHard);
  maxCargoKg = std::max(0.2, maxCargoKg);

  // Provided cargo should be limited to avoid draining the station economy too aggressively.
  const double invFrac = cargoProvided ? 0.28 : 0.80;

  std::vector<Cand> cands;
  cands.reserve(econ::kCommodityCount);

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const auto cid = static_cast<econ::CommodityId>(i);

    if (requireLegalAtDest && isIllegalCommodity(universeSeed, destFactionId, cid)) continue;

    const auto oq = econ::quote(originState, originModel, cid, 0.10);
    const auto dq = econ::quote(destState, destModel, cid, 0.10);

    if (oq.inventory < 1.0) continue;

    const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
    const int maxByKg = (int)std::floor(maxCargoKg / massKg);
    const int maxByInv = (int)std::floor(std::max(0.0, oq.inventory) * invFrac);
    const int maxUnits = std::min(maxUnitsHard, std::min(maxByKg, maxByInv));

    if (maxUnits < minUnits) continue;

    // Score cargo by demand: destination mid-price relative to origin mid-price.
    // Also bias slightly towards higher-value goods so boards aren't always dominated by
    // ultra-bulky low-value commodities.
    const double ratio = dq.mid / std::max(1e-6, oq.mid);
    const double valueBias = 0.65 + 0.35 * std::min(1.0, dq.mid / 120.0);
    const double score = ratio * valueBias;

    cands.push_back(Cand{cid, score, oq, dq, maxUnits});
  }

  if (cands.empty()) return false;

  std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
    if (a.score != b.score) return a.score > b.score;
    return (int)a.id < (int)b.id;
  });

  const int top = std::min<int>(3, (int)cands.size());
  const int pick = (top <= 1) ? 0 : rng.range(0, top - 1);
  const Cand& c = cands[(std::size_t)pick];

  const int units = rng.range(minUnits, c.maxUnits);
  out.commodity = c.id;
  out.units = units;
  out.originQ = c.o;
  out.destQ = c.d;
  return true;
}

static int activePassengerCount(const SaveGame& s) {
  int n = 0;
  for (const auto& m : s.missions) {
    if (m.completed || m.failed) continue;
    if (m.type != MissionType::Passenger) continue;
    n += std::max(0, (int)std::llround(m.units));
  }
  return n;
}

static core::u64 missionBoardSeed(StationId stationId, int dayStamp, core::u32 factionId) {
  // Similar to the in-game seeding approach: mix stable identifiers.
  // Uses large odd constants to improve bit diffusion.
  return (core::u64)stationId * 1469598103934665603ull ^ (core::u64)dayStamp * 1099511628211ull ^ (core::u64)factionId;
}

void refreshMissionOffers(Universe& universe,
                          const StarSystem& currentSystem,
                          const Station& dockedStation,
                          double timeDays,
                          double playerRep,
                          SaveGame& ioSave,
                          const MissionBoardParams& params) {
  const int dayStamp = (int)std::floor(timeDays);
  if (ioSave.missionOffersStationId == dockedStation.id && ioSave.missionOffersDayStamp == dayStamp) return;

  ioSave.missionOffersStationId = dockedStation.id;
  ioSave.missionOffersDayStamp = dayStamp;
  ioSave.missionOffers.clear();

  core::SplitMix64 mrng(missionBoardSeed(dockedStation.id, dayStamp, dockedStation.factionId));

  auto candidates = universe.queryNearby(currentSystem.stub.posLy, params.searchRadiusLy, params.maxCandidates);
  std::vector<SystemStub> dests;
  dests.reserve(candidates.size());
  for (const auto& s : candidates) {
    if (s.id == currentSystem.stub.id) continue;
    if (s.stationCount <= 0) continue;
    dests.push_back(s);
  }

  const int baseCount = params.baseOfferCount
                      + (playerRep >= params.repThreshold1 ? 1 : 0)
                      + (playerRep >= params.repThreshold2 ? 1 : 0);
  const int offerCount = std::clamp(baseCount, params.minOfferCount, params.maxOfferCount);

  const double repScale = 1.0 + std::clamp(playerRep, 0.0, 100.0) / 100.0 * params.repRewardBonus;

  // Used for pricing delivery missions that require purchasing cargo.
  const double feeEff = applyReputationToFeeRate(dockedStation.feeRate, playerRep);

  // Cache origin economy snapshot for this refresh.
  auto& originEcon = universe.stationEconomy(dockedStation, timeDays);

  // Ship-fit: tailor offer payload sizes to the player's current ship capacity.
  // For cargo missions that grant cargo immediately on acceptance, we use the *free*
  // cargo space so the board avoids generating offers the player can't accept right now.
  const double cargoFreeKg = std::max(0.0, ioSave.cargoCapacityKg - econ::cargoMassKg(ioSave.cargo));
  const double cargoCapKg = std::max(0.0, ioSave.cargoCapacityKg);
  const int freeSeats = std::max(0, ioSave.passengerSeats - activePassengerCount(ioSave));

    // Precompute cumulative weight thresholds once (the RNG is deterministic per board seed).
  const double tCourier = params.wCourier;
  const double tDelivery = tCourier + params.wDelivery;
  const double tMulti = tDelivery + params.wMultiDelivery;
  const double tSalvage = tMulti + params.wSalvage;
  const double tPassenger = tSalvage + params.wPassenger;
  const double tSmuggle = tPassenger + params.wSmuggle;
  const double tBountyScan = tSmuggle + params.wBountyScan;

  for (int i = 0; i < offerCount; ++i) {
    const double r = mrng.nextUnit();

    MissionType type = MissionType::BountyKill;
    if (r < tCourier) type = MissionType::Courier;
    else if (r < tDelivery) type = MissionType::Delivery;
    else if (r < tMulti) type = MissionType::MultiDelivery;
    else if (r < tSalvage) type = MissionType::Salvage;
    else if (r < tPassenger) type = MissionType::Passenger;
    else if (r < tSmuggle) type = MissionType::Smuggle;
    else if (r < tBountyScan) type = MissionType::BountyScan;
    else type = MissionType::BountyKill;

    // Multi-hop deliveries need at least two distinct candidate systems.
    if (type == MissionType::MultiDelivery && dests.size() < 2) {
      type = MissionType::Delivery;
    }

    // If we don't have any destination candidates, we can still offer local Salvage jobs.
    if (type != MissionType::Salvage && dests.empty()) {
      type = MissionType::Salvage;
    }

    Mission m{};
    m.id = 0; // offers are not persisted as "accepted" missions
    m.factionId = dockedStation.factionId;
    m.fromSystem = currentSystem.stub.id;
    m.fromStation = dockedStation.id;

    // Default to a local job that returns to the issuing station.
    m.toSystem = currentSystem.stub.id;
    m.toStation = dockedStation.id;

    // Default deadline (overridden per type)
    m.deadlineDay = timeDays + 1.0;

    double distLy = 0.0;
    const Station* destSt = nullptr;
    sim::SystemStub destStub{};

    if (type != MissionType::Salvage) {
      // Pick a random destination system/station from nearby candidates.
      const int pick = (int)(mrng.nextU32() % (core::u32)dests.size());
      destStub = dests[(std::size_t)pick];
      const auto& destSys = universe.getSystem(destStub.id, &destStub);
      if (destSys.stations.empty()) continue;
      destSt = &destSys.stations[(std::size_t)(mrng.nextU32() % (core::u32)destSys.stations.size())];

      m.toSystem = destStub.id;
      m.toStation = destSt->id;

      distLy = (destStub.posLy - currentSystem.stub.posLy).length();
      m.deadlineDay = timeDays + 1.0 + distLy / 20.0;
    } else {
      // Salvage is a local in-system job: make it snappy.
      m.deadlineDay = timeDays + 0.35 + mrng.range(0.05, 0.20);
    }

    if (type == MissionType::Courier) {
      m.type = MissionType::Courier;
      m.reward = 350.0 + distLy * 110.0;
    } else if (type == MissionType::Delivery) {
      m.type = MissionType::Delivery;
      m.cargoProvided = (mrng.nextUnit() < 0.25);

      // Ship-fit: cargo missions grant cargo immediately on acceptance, so only generate
      // jobs that fit into the player's current free cargo space.
      const double freeKg = cargoFreeKg;
      if (freeKg + 1e-9 < 0.2) {
        // Hold is effectively full (or too small to carry any commodity); keep the board usable.
        m.type = MissionType::Courier;
        m.reward = 330.0 + distLy * 115.0;
      } else {
        const double maxCargoKg = std::min(260.0, freeKg * 0.95);
        const int minUnits = std::clamp((int)std::floor(freeKg / 30.0), 1, 12);
        const int maxUnitsHard = std::clamp((int)std::floor(freeKg / 1.2), 6, 180);

        // Economy-aware cargo selection. Ensure we don't ask the player to deliver contraband
        // as a "legal" delivery job.
        CargoMissionSpec spec{};
        bool ok = false;
        if (destSt) {
          auto& destEcon = universe.stationEconomy(*destSt, timeDays);
          ok = pickCargoMission(mrng,
                                universe.seed(),
                                originEcon,
                                dockedStation.economyModel,
                                destEcon,
                                destSt->economyModel,
                                destSt->factionId,
                                /*requireLegalAtDest=*/true,
                                /*cargoProvided=*/m.cargoProvided,
                                /*maxCargoKg=*/maxCargoKg,
                                /*minUnits=*/minUnits,
                                /*maxUnitsHard=*/maxUnitsHard,
                                spec);
        }

        if (!ok) {
          // Fallback: if the economy can't support a cargo job right now, degrade to courier.
          m.type = MissionType::Courier;
          m.reward = 330.0 + distLy * 115.0;
        } else {
          m.commodity = spec.commodity;
          m.units = (double)spec.units;

          const double base = 250.0 + distLy * 58.0 + mrng.range(0.0, 140.0);

          if (m.cargoProvided) {
            // Provided cargo: payout is more "service" oriented; scale with destination price.
            m.reward = base + spec.destQ.bid * (double)spec.units * 0.55;
          } else {
            // Purchased cargo: ensure missions are meaningfully profitable even for high-value goods.
            const double cost = spec.originQ.ask * (double)spec.units * (1.0 + feeEff);
            const double margin = 0.18 + mrng.range(0.0, 0.08);
            m.reward = base + cost * (1.0 + margin);
          }
        }
      }
    } else if (type == MissionType::MultiDelivery) {
      m.type = MissionType::MultiDelivery;
      m.cargoProvided = (mrng.nextUnit() < 0.35);

      // Multi-hop deliveries grant cargo immediately on acceptance; size them to the player's
      // current free cargo capacity to avoid impossible offers.
      const double freeKg = cargoFreeKg;
      if (freeKg + 1e-9 < 0.2) {
        m.type = MissionType::Courier;
        m.reward = 360.0 + distLy * 120.0;
      } else {
        const double maxCargoKg = std::min(240.0, freeKg * 0.90);
        const int minUnits = std::clamp((int)std::floor(freeKg / 35.0), 1, 10);
        const int maxUnitsHard = std::clamp((int)std::floor(freeKg / 1.4), 6, 160);

        // Pick a cargo load that makes sense for the final destination's market.
        CargoMissionSpec spec{};
        bool ok = false;
        if (destSt) {
          auto& destEcon = universe.stationEconomy(*destSt, timeDays);
          ok = pickCargoMission(mrng,
                                universe.seed(),
                                originEcon,
                                dockedStation.economyModel,
                                destEcon,
                                destSt->economyModel,
                                destSt->factionId,
                                /*requireLegalAtDest=*/true,
                                /*cargoProvided=*/m.cargoProvided,
                                /*maxCargoKg=*/maxCargoKg,
                                /*minUnits=*/minUnits,
                                /*maxUnitsHard=*/maxUnitsHard,
                                spec);
        }

        if (!ok) {
          // Fallback: keep the board populated.
          m.type = MissionType::Courier;
          m.reward = 360.0 + distLy * 120.0;
        } else {
          m.commodity = spec.commodity;
          m.units = (double)spec.units;

          // Pick a via system/station. Prefer routes that aren't extreme detours.
          SystemStub viaStub{};
          const Station* viaSt = nullptr;
          double routeLy = 0.0;

          const double directLy = std::max(1e-6, distLy);
          for (int tries = 0; tries < 18; ++tries) {
            const auto candStub = dests[(std::size_t)(mrng.nextU32() % (core::u32)dests.size())];
            if (candStub.id == destStub.id) continue;

            const double d1 = (candStub.posLy - currentSystem.stub.posLy).length();
            const double d2 = (destStub.posLy - candStub.posLy).length();
            const double route = d1 + d2;
            const double detour = route / directLy;
            if (detour > 2.25) continue;

            const auto& candSys = universe.getSystem(candStub.id, &candStub);
            if (candSys.stations.empty()) continue;
            const auto& candSt = candSys.stations[(std::size_t)(mrng.nextU32() % (core::u32)candSys.stations.size())];

            if (!viaSt || route < routeLy) {
              viaStub = candStub;
              viaSt = &candSt;
              routeLy = route;
            }
          }

          if (!viaSt) {
            // Couldn't find a reasonable intermediate stop; downgrade to a normal delivery.
            m.type = MissionType::Delivery;

            const double base = 270.0 + distLy * 60.0 + mrng.range(0.0, 160.0);
            if (m.cargoProvided) {
              m.reward = base + spec.destQ.bid * (double)spec.units * 0.55;
            } else {
              const double cost = spec.originQ.ask * (double)spec.units * (1.0 + feeEff);
              const double margin = 0.18 + mrng.range(0.0, 0.08);
              m.reward = base + cost * (1.0 + margin);
            }

            // Keep the direct-leg deadline picked above.
          } else {
            m.viaSystem = viaStub.id;
            m.viaStation = viaSt->id;

            // Reward/deadline scale with the actual route distance (origin->via + via->dest),
            // plus a small "extra stop" premium.
            const double base = 380.0 + routeLy * 70.0 + mrng.range(0.0, 220.0);

            if (m.cargoProvided) {
              m.reward = base + spec.destQ.bid * (double)spec.units * 0.65;
            } else {
              const double cost = spec.originQ.ask * (double)spec.units * (1.0 + feeEff);
              const double margin = 0.22 + mrng.range(0.0, 0.10);
              m.reward = base + cost * (1.0 + margin);
            }

            m.reward += 140.0 + mrng.range(0.0, 70.0);
            m.deadlineDay = timeDays + 2.05 + routeLy / 18.0 + mrng.range(0.0, 0.25);
          }
        }
      }
    } else if (type == MissionType::Salvage) {
      m.type = MissionType::Salvage;

      // Salvage jobs are local: recover loose goods from a mission derelict signal and return here.
      static const std::array<econ::CommodityId, 4> kSalvage = {
        econ::CommodityId::Machinery,
        econ::CommodityId::Electronics,
        econ::CommodityId::Metals,
        econ::CommodityId::Luxury,
      };

      // Ensure the requested salvage can fit in the player's cargo hold (at completion time).
      const double capKg = cargoCapKg;
      bool ok = false;

      for (int tries = 0; tries < 8; ++tries) {
        m.commodity = kSalvage[(std::size_t)(mrng.nextU32() % (core::u32)kSalvage.size())];
        const double massKg = std::max(1e-6, econ::commodityDef(m.commodity).massKg);
        const int maxByKg = (int)std::floor(capKg / massKg + 1e-9);
        if (maxByKg >= 1) {
          const int maxUnits = std::clamp(maxByKg, 1, 18);
          const int minUnits = std::min(3, maxUnits);
          m.units = (double)mrng.range(minUnits, maxUnits);
          ok = true;
          break;
        }
      }

      if (!ok) {
        // If the ship can't carry any salvage at all, keep the board usable.
        m.type = MissionType::Courier;
        m.reward = 360.0;
      } else {
        const double base = econ::commodityDef(m.commodity).basePrice;
        m.reward = 520.0 + (double)m.units * base * 1.05 + mrng.range(0.0, 180.0);

        // Reuse targetNpcId as a stable "signal id" so the game can spawn/track the mission site.
        m.targetNpcId = (mrng.nextU64() | 0x8000000000000000ull);
        m.cargoProvided = false;
      }
    } else if (type == MissionType::Passenger) {
      // Ship-fit: only offer passenger party sizes that the player can accept right now.
      const int maxParty = std::min(6, freeSeats);
      if (maxParty <= 0) {
        m.type = MissionType::Courier;
        m.reward = 340.0 + distLy * 118.0;
      } else {
        m.type = MissionType::Passenger;
        // units is interpreted as "passenger count" for Passenger missions.
        m.units = (double)mrng.range(1, maxParty);
        // Reward scales with distance and party size.
        m.reward = 450.0 + distLy * 125.0 + (double)m.units * 85.0;
        // Slightly tighter deadline than courier to encourage routing decisions.
        m.deadlineDay = timeDays + 0.9 + distLy / 24.0;
        m.cargoProvided = false;
      }
    } else if (type == MissionType::Smuggle) {
      m.type = MissionType::Smuggle;

      // Smuggling jobs should only exist where there is actual contraband.
      const core::u32 destFaction = destSt ? destSt->factionId : 0;
      const core::u32 mask = illegalCommodityMask(universe.seed(), destFaction);

      std::vector<econ::CommodityId> illegal;
      illegal.reserve(econ::kCommodityCount);
      for (std::size_t c = 0; c < econ::kCommodityCount; ++c) {
        if ((mask & ((core::u32)1u << (core::u32)c)) != 0u) illegal.push_back((econ::CommodityId)c);
      }

      if (destFaction == 0 || illegal.empty()) {
        // No contraband rules here; downgrade to a normal courier job.
        m.type = MissionType::Courier;
        m.reward = 340.0 + distLy * 118.0;
      } else {
        // Ship-fit: smuggling grants cargo immediately on acceptance, so filter to commodities that
        // can actually fit in the player's *current* free hold.
        const double freeKg = cargoFreeKg;

        std::vector<econ::CommodityId> fitIllegal;
        fitIllegal.reserve(illegal.size());
        for (const auto cid : illegal) {
          const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
          if (massKg <= freeKg * 0.95 + 1e-9) fitIllegal.push_back(cid);
        }

        if (fitIllegal.empty()) {
          // Not enough free capacity to accept even 1 unit of any contraband.
          m.type = MissionType::Courier;
          m.reward = 340.0 + distLy * 118.0;
        } else {
          m.commodity = fitIllegal[(std::size_t)(mrng.nextU32() % (core::u32)fitIllegal.size())];

          const double massKg = std::max(1e-6, econ::commodityDef(m.commodity).massKg);
          const int maxByKg = (int)std::floor((freeKg * 0.95) / massKg + 1e-9);
          const int maxUnits = std::clamp(maxByKg, 1, 55);

          const int minUnits = std::min(8, maxUnits);
          const int units = mrng.range(minUnits, maxUnits);
          m.units = (double)units;

          // Smuggling cargo is provided by the contact (not drawn from open-market inventory).
          m.cargoProvided = true;

          // Reward scales with destination price (risk/interest) + distance.
          double destMid = econ::commodityDef(m.commodity).basePrice;
          if (destSt) {
            auto& destEcon = universe.stationEconomy(*destSt, timeDays);
            destMid = econ::quote(destEcon, destSt->economyModel, m.commodity, 0.10).mid;
          }

          m.reward = 720.0
                   + distLy * 170.0
                   + destMid * (double)units * 0.28
                   + mrng.range(0.0, 220.0);
          m.deadlineDay = timeDays + 1.15 + distLy / 20.0;
        }
      }
    } else if (type == MissionType::BountyScan) {
      m.type = MissionType::BountyScan;
      m.targetNpcId = std::max<core::u64>(1, mrng.nextU64());
      m.reward = 900.0 + distLy * 80.0;
      m.deadlineDay = timeDays + 1.5 + distLy / 22.0;
    } else {
      m.type = MissionType::BountyKill;
      m.targetNpcId = std::max<core::u64>(1, mrng.nextU64());
      m.reward = 1400.0 + distLy * 90.0;
      m.deadlineDay = timeDays + 1.8 + distLy / 20.0;
    }

    m.reward *= repScale;
    ioSave.missionOffers.push_back(std::move(m));
  }
}

bool acceptMissionOffer(Universe& universe,
                        const Station& dockedStation,
                        double timeDays,
                        double playerRep,
                        SaveGame& ioSave,
                        std::size_t offerIndex,
                        std::string* outError) {
  if (offerIndex >= ioSave.missionOffers.size()) {
    if (outError) *outError = "Offer index out of range.";
    return false;
  }

  Mission m = ioSave.missionOffers[offerIndex];

  const double feeEff = applyReputationToFeeRate(dockedStation.feeRate, playerRep);
  auto& stEcon = universe.stationEconomy(dockedStation, timeDays);

  // Delivery-like missions need cargo acquisition.
  const bool isDelivery = (m.type == MissionType::Delivery || m.type == MissionType::MultiDelivery || m.type == MissionType::Smuggle);
  if (isDelivery && m.units > 0.0) {
    const econ::CommodityId cid = m.commodity;
    const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
    const double needKg = massKg * (double)m.units;
    if (econ::cargoMassKg(ioSave.cargo) + needKg > ioSave.cargoCapacityKg + 1e-6) {
      if (outError) *outError = "Not enough cargo capacity for this delivery.";
      return false;
    }

    if (m.type == MissionType::Smuggle) {
      // Smuggling cargo is provided by the job contact (not pulled from the legal station market).
      ioSave.cargo[(std::size_t)cid] += (double)m.units;
    } else if (m.cargoProvided) {
      // Provided cargo should come from real station inventory.
      const double taken = econ::takeInventory(stEcon, dockedStation.economyModel, cid, (double)m.units);
      if (taken + 1e-6 < (double)m.units) {
        // Restore any partial transfer.
        econ::addInventory(stEcon, dockedStation.economyModel, cid, taken);
        if (outError) *outError = "Station can't provide that much cargo right now.";
        return false;
      }
      ioSave.cargo[(std::size_t)cid] += taken;
    } else {
      const auto q = econ::quote(stEcon, dockedStation.economyModel, cid, 0.10);
      const double totalCost = q.ask * (double)m.units * (1.0 + feeEff);
      if (q.inventory + 1e-6 < (double)m.units) {
        if (outError) *outError = "Station doesn't have enough inventory to stock this delivery.";
        return false;
      }
      if (ioSave.credits + 1e-6 < totalCost) {
        if (outError) *outError = "Not enough credits to buy the delivery cargo.";
        return false;
      }
      const auto tr = econ::buy(stEcon, dockedStation.economyModel, cid, (double)m.units, ioSave.credits, 0.10, feeEff);
      if (!tr.ok) {
        if (outError) *outError = tr.reason ? tr.reason : "Trade failed.";
        return false;
      }
      ioSave.cargo[(std::size_t)cid] += tr.unitsDelta;
    }
  }

  // Salvage missions don't grant cargo up-front, but we should still ensure the requested salvage will
  // fit in the player's hold (otherwise the job is impossible to complete).
  if (m.type == MissionType::Salvage && m.units > 0.0) {
    const econ::CommodityId cid = m.commodity;
    const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
    const double needKg = massKg * (double)m.units;
    if (needKg > ioSave.cargoCapacityKg + 1e-6) {
      if (outError) *outError = "This salvage won't fit in your cargo hold.";
      return false;
    }
  }

  // Passenger missions require cabin capacity.
  if (m.type == MissionType::Passenger) {
    const int party = std::max(0, (int)std::llround(m.units));
    const int used = activePassengerCount(ioSave);
    const int cap = std::max(0, ioSave.passengerSeats);
    if (party <= 0) {
      if (outError) *outError = "Passenger party size is invalid.";
      return false;
    }
    if (used + party > cap) {
      if (outError) *outError = "Not enough passenger seats for this job.";
      return false;
    }
  }

  m.id = ioSave.nextMissionId++;
  ioSave.missions.push_back(m);
  ioSave.missionOffers.erase(ioSave.missionOffers.begin() + (std::ptrdiff_t)offerIndex);
  return true;
}

MissionTickResult tickMissionDeadlines(SaveGame& ioSave, double timeDays, double repPenaltyOnFail) {
  MissionTickResult r{};
  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;
    if (m.deadlineDay > 0.0 && timeDays > m.deadlineDay) {
      m.failed = true;
      addReputation(ioSave, m.factionId, repPenaltyOnFail);
      ++r.failed;
    }
  }
  return r;
}

MissionDockResult tryCompleteMissionsAtDock(Universe& universe,
                                            const StarSystem& currentSystem,
                                            const Station& dockedStation,
                                            double timeDays,
                                            SaveGame& ioSave,
                                            double repRewardOnComplete) {
  MissionDockResult r{};
  const SystemId sysId = currentSystem.stub.id;
  const StationId here = dockedStation.id;

  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;

    const bool atFinal = (sysId == m.toSystem && here == m.toStation);
    const bool atVia = (m.viaSystem != 0 && sysId == m.viaSystem && here == m.viaStation);

    if (m.type == MissionType::Courier) {
      if (atFinal) {
        m.completed = true;
        ioSave.credits += m.reward;
        addReputation(ioSave, m.factionId, repRewardOnComplete);
        ++r.completed;
      }
    } else if (m.type == MissionType::Passenger) {
      if (atFinal) {
        m.completed = true;
        ioSave.credits += m.reward;
        addReputation(ioSave, m.factionId, repRewardOnComplete);
        ++r.completed;
      }
    } else if (m.type == MissionType::Salvage) {
      // Salvage jobs require that the player has actually visited the mission site (m.scanned),
      // then returned with the requested goods.
      if (atFinal && m.scanned) {
        const econ::CommodityId cid = m.commodity;
        const double have = ioSave.cargo[(std::size_t)cid];
        if (have + 1e-6 >= m.units) {
          ioSave.cargo[(std::size_t)cid] -= m.units;
          if (ioSave.cargo[(std::size_t)cid] < 1e-6) ioSave.cargo[(std::size_t)cid] = 0.0;

          // Feed recovered salvage into the local market (best-effort).
          auto& econHere = universe.stationEconomy(dockedStation, timeDays);
          econ::addInventory(econHere, dockedStation.economyModel, cid, m.units);

          m.completed = true;
          ioSave.credits += m.reward;
          addReputation(ioSave, m.factionId, repRewardOnComplete);
          ++r.completed;
        }
      }
    } else if (m.type == MissionType::Delivery || m.type == MissionType::MultiDelivery || m.type == MissionType::Smuggle) {
      if (m.type == MissionType::MultiDelivery && m.viaSystem != 0 && m.leg == 0 && atVia) {
        m.leg = 1;
        ++r.progressedMultiLeg;
      } else if (atFinal && (m.viaSystem == 0 || m.leg >= 1)) {
        const econ::CommodityId cid = m.commodity;
        const double have = ioSave.cargo[(std::size_t)cid];
        if (have + 1e-6 >= m.units) {
          ioSave.cargo[(std::size_t)cid] -= m.units;
          if (ioSave.cargo[(std::size_t)cid] < 1e-6) ioSave.cargo[(std::size_t)cid] = 0.0;

          if (m.type != MissionType::Smuggle) {
            // Feed delivered cargo back into the destination market (best-effort).
            auto& econHere = universe.stationEconomy(dockedStation, timeDays);
            econ::addInventory(econHere, dockedStation.economyModel, cid, m.units);
          }

          m.completed = true;
          ioSave.credits += m.reward;
          addReputation(ioSave, m.factionId, repRewardOnComplete);
          ++r.completed;
        }
      }
    }
  }

  return r;
}


MissionEventResult tryCompleteBountyScan(SaveGame& ioSave,
                                         SystemId currentSystemId,
                                         core::u64 targetNpcId,
                                         double repRewardOnComplete) {
  MissionEventResult r{};
  if (currentSystemId == 0 || targetNpcId == 0) return r;

  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;
    if (m.type != MissionType::BountyScan) continue;
    if (m.toSystem != currentSystemId) continue;
    if (m.targetNpcId != targetNpcId) continue;

    m.scanned = true;
    m.completed = true;
    ioSave.credits += m.reward;
    addReputation(ioSave, m.factionId, repRewardOnComplete);
    ++r.completed;
    r.rewardCr += m.reward;
  }

  return r;
}

MissionEventResult tryCompleteBountyKill(SaveGame& ioSave,
                                         SystemId currentSystemId,
                                         core::u64 targetNpcId,
                                         double repRewardOnComplete) {
  MissionEventResult r{};
  if (currentSystemId == 0 || targetNpcId == 0) return r;

  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;
    if (m.type != MissionType::BountyKill) continue;
    if (m.toSystem != currentSystemId) continue;
    if (m.targetNpcId != targetNpcId) continue;

    m.completed = true;
    ioSave.credits += m.reward;
    addReputation(ioSave, m.factionId, repRewardOnComplete);
    ++r.completed;
    r.rewardCr += m.reward;
  }

  return r;
}

} // namespace stellar::sim
