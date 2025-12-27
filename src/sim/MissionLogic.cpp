#include "stellar/sim/MissionLogic.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Cargo.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/Reputation.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace stellar::sim {

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

  for (int i = 0; i < offerCount && !dests.empty(); ++i) {
    const int pick = (int)(mrng.nextU32() % (core::u32)dests.size());
    const auto destStub = dests[(std::size_t)pick];
    const auto& destSys = universe.getSystem(destStub.id, &destStub);
    const auto& destSt = destSys.stations[(std::size_t)(mrng.nextU32() % (core::u32)destSys.stations.size())];

    const double distLy = (destStub.posLy - currentSystem.stub.posLy).length();

    Mission m{};
    m.id = 0; // offers are not persisted as "accepted" missions
    m.factionId = dockedStation.factionId;
    m.fromSystem = currentSystem.stub.id;
    m.fromStation = dockedStation.id;
    m.toSystem = destStub.id;
    m.toStation = destSt.id;
    m.deadlineDay = timeDays + 1.0 + distLy / 20.0;

    const double r = mrng.nextUnit();
    const double tCourier = params.wCourier;
    const double tDelivery = tCourier + params.wDelivery;
    const double tMulti = tDelivery + params.wMultiDelivery;
    const double tPassenger = tMulti + params.wPassenger;
    const double tSmuggle = tPassenger + params.wSmuggle;
    const double tBountyScan = tSmuggle + params.wBountyScan;

    if (r < tCourier) {
      m.type = MissionType::Courier;
      m.reward = 350.0 + distLy * 110.0;
    } else if (r < tDelivery) {
      m.type = MissionType::Delivery;
      m.commodity = (econ::CommodityId)(mrng.nextU32() % (core::u32)econ::kCommodityCount);
      m.units = 25 + (mrng.nextU32() % 120);
      m.reward = 250.0 + distLy * 120.0 + (double)m.units * 6.0;
      m.cargoProvided = (mrng.nextUnit() < 0.25);
    } else if (r < tMulti && dests.size() >= 2) {
      m.type = MissionType::MultiDelivery;
      m.commodity = (econ::CommodityId)(mrng.nextU32() % (core::u32)econ::kCommodityCount);
      m.units = 20 + (mrng.nextU32() % 100);
      m.reward = 400.0 + distLy * 150.0 + (double)m.units * 8.0;
      m.cargoProvided = (mrng.nextUnit() < 0.35);

      // Pick a via system/station (try to avoid picking the same destination).
      for (int tries = 0; tries < 8; ++tries) {
        const auto viaStub = dests[(std::size_t)(mrng.nextU32() % (core::u32)dests.size())];
        if (viaStub.id == destStub.id) continue;
        const auto& viaSys = universe.getSystem(viaStub.id, &viaStub);
        if (viaSys.stations.empty()) continue;
        const auto& viaSt = viaSys.stations[(std::size_t)(mrng.nextU32() % (core::u32)viaSys.stations.size())];
        m.viaSystem = viaStub.id;
        m.viaStation = viaSt.id;
        break;
      }
      // Slightly longer default deadline for multi-hop.
      m.deadlineDay = timeDays + 2.0 + distLy / 18.0;
    } else if (r < tPassenger) {
      m.type = MissionType::Passenger;
      // units is interpreted as "passenger count" for Passenger missions.
      m.units = 1 + (mrng.nextU32() % 8);
      // Reward scales with distance and party size.
      m.reward = 450.0 + distLy * 125.0 + (double)m.units * 85.0;
      // Slightly tighter deadline than courier to encourage routing decisions.
      m.deadlineDay = timeDays + 0.9 + distLy / 24.0;
      m.cargoProvided = false;
    } else if (r < tSmuggle) {
      m.type = MissionType::Smuggle;

      // Pick a commodity that is illegal in the destination's jurisdiction.
      // Smuggling jobs are treated as "cargo provided" by the contact (not drawn from open-market inventory).
      const core::u32 destFaction = destSt.factionId;
      const core::u32 mask = illegalCommodityMask(universe.seed(), destFaction);

      std::vector<econ::CommodityId> illegal;
      illegal.reserve(econ::kCommodityCount);
      for (std::size_t c = 0; c < econ::kCommodityCount; ++c) {
        if ((mask & ((core::u32)1u << (core::u32)c)) != 0u) illegal.push_back((econ::CommodityId)c);
      }

      if (!illegal.empty() && destFaction != 0) {
        m.commodity = illegal[(std::size_t)(mrng.nextU32() % (core::u32)illegal.size())];
      } else {
        // Fallback: if the destination has no contraband rules, pick any commodity.
        m.commodity = (econ::CommodityId)(mrng.nextU32() % (core::u32)econ::kCommodityCount);
      }

      // Contraband shipments are smaller but higher margin.
      m.units = 10 + (mrng.nextU32() % 45);
      m.reward = 650.0 + distLy * 155.0 + (double)m.units * 18.0;
      m.deadlineDay = timeDays + 1.2 + distLy / 19.0;
      m.cargoProvided = true;
    } else if (r < tBountyScan) {
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

} // namespace stellar::sim
