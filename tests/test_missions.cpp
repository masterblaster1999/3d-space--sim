#include "stellar/sim/MissionLogic.h"

#include "stellar/econ/Cargo.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Reputation.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <iostream>

int test_missions() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  Universe u(1337);

  // Find a nearby system with at least one station.
  const auto stubs = u.queryNearby({0, 0, 0}, 200.0, 64);
  if (stubs.empty()) {
    std::cerr << "[test_missions] no systems returned\n";
    return 1;
  }

  const SystemStub* chosenStub = nullptr;
  for (const auto& s : stubs) {
    if (s.stationCount > 0) { chosenStub = &s; break; }
  }
  if (!chosenStub) {
    std::cerr << "[test_missions] no systems with stations in query\n";
    return 1;
  }

  const auto& sys = u.getSystem(chosenStub->id, chosenStub);
  if (sys.stations.empty()) {
    std::cerr << "[test_missions] stub said stations but system had none\n";
    return 1;
  }

  const auto& station = sys.stations.front();

  SaveGame s{};
  s.seed = 1337;
  s.timeDays = 12.25;
  s.credits = 50000.0;
  s.cargoCapacityKg = 420.0;

  // Deterministic mission board generation.
  {
    SaveGame a = s;
    SaveGame b = s;

    refreshMissionOffers(u, sys, station, a.timeDays, /*rep*/0.0, a);
    refreshMissionOffers(u, sys, station, b.timeDays, /*rep*/0.0, b);

    if (a.missionOffers.empty()) {
      std::cerr << "[test_missions] expected some offers\n";
      ++fails;
    }

    if (a.missionOffers.size() != b.missionOffers.size()) {
      std::cerr << "[test_missions] deterministic offers size mismatch\n";
      ++fails;
    } else if (!a.missionOffers.empty()) {
      // Spot-check a few stable fields.
      const auto& ma = a.missionOffers.front();
      const auto& mb = b.missionOffers.front();
      if (ma.type != mb.type || ma.toSystem != mb.toSystem || ma.toStation != mb.toStation) {
        std::cerr << "[test_missions] deterministic offers content mismatch\n";
        ++fails;
      }
    }

    // Basic sanity checks on generated offers (legality + profitability for buy-to-deliver jobs).
    // These are intentionally lightweight: they catch regressions where mission rewards end up
    // below the expected acquisition cost for high-value cargo.
    auto findStationById = [](const StarSystem& sys, sim::StationId id) -> const sim::Station* {
      for (const auto& st : sys.stations) {
        if (st.id == id) return &st;
      }
      return nullptr;
    };

    auto& originEcon = u.stationEconomy(station, a.timeDays);
    const double feeEff = applyReputationToFeeRate(station.feeRate, /*rep*/0.0);

    double freeKg = a.cargoCapacityKg - econ::cargoMassKg(a.cargo);
    if (freeKg < 0.0) freeKg = 0.0;

    for (const auto& m : a.missionOffers) {
      // Ship-fit sanity: generated offers should be accept-able given the player's current ship capacity.
      if (m.type == MissionType::Passenger) {
        const int party = std::max(0, (int)std::llround(m.units));
        if (party > a.passengerSeats) {
          std::cerr << "[test_missions] passenger offer exceeds seat capacity\n";
          ++fails;
        }
      }

      if ((m.type == MissionType::Delivery || m.type == MissionType::MultiDelivery || m.type == MissionType::Smuggle) && m.units > 0.0) {
        const double massKg = econ::commodityDef(m.commodity).massKg;
        const double needKg = massKg * m.units;
        if (needKg > freeKg + 1e-6) {
          std::cerr << "[test_missions] cargo offer exceeds free cargo capacity\n";
          ++fails;
        }
      }

      if (m.type == MissionType::Salvage && m.units > 0.0) {
        const double massKg = econ::commodityDef(m.commodity).massKg;
        const double needKg = massKg * m.units;
        if (needKg > a.cargoCapacityKg + 1e-6) {
          std::cerr << "[test_missions] salvage offer exceeds cargo capacity\n";
          ++fails;
        }
      }

      if (m.type == MissionType::Smuggle) {
        const auto& ds = u.getSystem(m.toSystem);
        const auto* dst = findStationById(ds, m.toStation);
        if (dst && dst->factionId != 0) {
          if (!isIllegalCommodity(u.seed(), dst->factionId, m.commodity)) {
            std::cerr << "[test_missions] smuggle mission commodity not illegal at destination\n";
            ++fails;
          }
        }
      }

      if (m.type == MissionType::Delivery || m.type == MissionType::MultiDelivery) {
        const auto& ds = u.getSystem(m.toSystem);
        const auto* dst = findStationById(ds, m.toStation);
        if (dst && dst->factionId != 0) {
          if (isIllegalCommodity(u.seed(), dst->factionId, m.commodity)) {
            std::cerr << "[test_missions] delivery mission requests illegal cargo at destination\n";
            ++fails;
          }
        }

        if (!m.cargoProvided) {
          const auto q = econ::quote(originEcon, station.economyModel, m.commodity, 0.10);
          const double cost = q.ask * m.units * (1.0 + feeEff);
          if (m.reward + 1e-6 < cost * 1.05) {
            std::cerr << "[test_missions] delivery mission reward too low vs acquisition cost\n";
            ++fails;
          }
        }
      }
    }
  }

  // Ship-fit mission boards: passenger offers should not appear if the ship has zero seats.
  {
    SaveGame p = s;
    p.passengerSeats = 0;

    refreshMissionOffers(u, sys, station, p.timeDays, /*rep*/0.0, p);
    for (const auto& m : p.missionOffers) {
      if (m.type == MissionType::Passenger) {
        std::cerr << "[test_missions] expected no passenger offers when passengerSeats == 0\n";
        ++fails;
      }
    }
  }

  // Ship-fit mission boards: cargo-bearing offers should respect small cargo capacities.
  {
    SaveGame p = s;
    p.cargoCapacityKg = 12.0;
    // Keep cargo empty so free space == capacity.

    refreshMissionOffers(u, sys, station, p.timeDays, /*rep*/0.0, p);

    double freeKg = p.cargoCapacityKg - econ::cargoMassKg(p.cargo);
    if (freeKg < 0.0) freeKg = 0.0;

    for (const auto& m : p.missionOffers) {
      if ((m.type == MissionType::Delivery || m.type == MissionType::MultiDelivery || m.type == MissionType::Smuggle) && m.units > 0.0) {
        const double massKg = econ::commodityDef(m.commodity).massKg;
        const double needKg = massKg * m.units;
        if (needKg > freeKg + 1e-6) {
          std::cerr << "[test_missions] small-hold cargo offer exceeds free cargo capacity\n";
          ++fails;
        }
      }

      if (m.type == MissionType::Salvage && m.units > 0.0) {
        const double massKg = econ::commodityDef(m.commodity).massKg;
        const double needKg = massKg * m.units;
        if (needKg > p.cargoCapacityKg + 1e-6) {
          std::cerr << "[test_missions] small-hold salvage offer exceeds cargo capacity\n";
          ++fails;
        }
      }

      if (m.type == MissionType::Passenger) {
        const int party = std::max(0, (int)std::llround(m.units));
        if (party > p.passengerSeats) {
          std::cerr << "[test_missions] small-hold passenger offer exceeds seat capacity\n";
          ++fails;
        }
      }
    }
  }

  // Accept/complete a fabricated delivery mission to exercise inventory/credits paths.
  {
    SaveGame p = s;
    setReputation(p, station.factionId, 10.0);
    const double rep = getReputation(p, station.factionId);

    // Ensure station has enough inventory.
    auto& econState = u.stationEconomy(station, p.timeDays);
    econ::addInventory(econState, station.economyModel, econ::CommodityId::Food, 1000.0);

    Mission offer{};
    offer.type = MissionType::Delivery;
    offer.factionId = station.factionId;
    offer.fromSystem = sys.stub.id;
    offer.fromStation = station.id;
    offer.toSystem = sys.stub.id;
    offer.toStation = station.id;
    offer.commodity = econ::CommodityId::Food;
    offer.units = 10.0;
    offer.reward = 1234.0;
    offer.deadlineDay = p.timeDays + 10.0;
    offer.cargoProvided = true;
    p.missionOffers = {offer};
    p.missionOffersStationId = station.id;
    p.missionOffersDayStamp = (int)std::floor(p.timeDays);

    std::string err;
    if (!acceptMissionOffer(u, station, p.timeDays, rep, p, 0, &err)) {
      std::cerr << "[test_missions] acceptMissionOffer failed: " << err << "\n";
      ++fails;
    } else {
      if (p.missions.size() != 1 || !p.missionOffers.empty()) {
        std::cerr << "[test_missions] acceptance did not move offer -> missions\n";
        ++fails;
      }
      if (p.cargo[(std::size_t)econ::CommodityId::Food] + 1e-6 < 10.0) {
        std::cerr << "[test_missions] acceptance did not grant cargo\n";
        ++fails;
      }
    }

    // Dock completion should consume cargo and reward credits.
    const double beforeCredits = p.credits;
    const double beforeFood = p.cargo[(std::size_t)econ::CommodityId::Food];
    const auto res = tryCompleteMissionsAtDock(u, sys, station, p.timeDays, p);

    if (res.completed != 1) {
      std::cerr << "[test_missions] expected completion count 1\n";
      ++fails;
    }
    if (p.credits <= beforeCredits) {
      std::cerr << "[test_missions] expected credits to increase on completion\n";
      ++fails;
    }
    if (p.cargo[(std::size_t)econ::CommodityId::Food] >= beforeFood - 1e-6) {
      std::cerr << "[test_missions] expected cargo to be consumed on completion\n";
      ++fails;
    }
  }


  // Smuggling missions: cargo is granted without touching the station's legal inventory.
  // Completion consumes cargo but does not feed inventory back into the destination market.
  {
    SaveGame p = s;
    setReputation(p, station.factionId, 10.0);
    const double rep = getReputation(p, station.factionId);

    auto& econState = u.stationEconomy(station, p.timeDays);
    const double invBefore = econState.inventory[(int)econ::CommodityId::Luxury];

    Mission offer{};
    offer.type = MissionType::Smuggle;
    offer.factionId = station.factionId;
    offer.fromSystem = sys.stub.id;
    offer.fromStation = station.id;
    offer.toSystem = sys.stub.id;
    offer.toStation = station.id;
    offer.commodity = econ::CommodityId::Luxury;
    offer.units = 5.0;
    offer.reward = 4321.0;
    offer.deadlineDay = p.timeDays + 10.0;
    offer.cargoProvided = true;
    p.missionOffers = {offer};
    p.missionOffersStationId = station.id;
    p.missionOffersDayStamp = (int)std::floor(p.timeDays);

    std::string err;
    const double creditsBefore = p.credits;
    if (!acceptMissionOffer(u, station, p.timeDays, rep, p, 0, &err)) {
      std::cerr << "[test_missions] acceptMissionOffer (Smuggle) failed: " << err << "\n";
      ++fails;
    } else {
      if (p.cargo[(std::size_t)econ::CommodityId::Luxury] + 1e-6 < 5.0) {
        std::cerr << "[test_missions] smuggle acceptance did not grant cargo\n";
        ++fails;
      }
      if (std::abs(p.credits - creditsBefore) > 1e-6) {
        std::cerr << "[test_missions] smuggle acceptance unexpectedly changed credits\n";
        ++fails;
      }
      const double invAfter = econState.inventory[(int)econ::CommodityId::Luxury];
      if (std::abs(invAfter - invBefore) > 1e-6) {
        std::cerr << "[test_missions] smuggle acceptance unexpectedly changed legal inventory\n";
        ++fails;
      }
    }

    // Dock completion should consume cargo and reward credits, but not alter station inventory.
    const double invBeforeComplete = econState.inventory[(int)econ::CommodityId::Luxury];
    const double beforeCredits = p.credits;
    const auto res = tryCompleteMissionsAtDock(u, sys, station, p.timeDays, p);

    if (res.completed != 1) {
      std::cerr << "[test_missions] expected smuggle completion count 1\n";
      ++fails;
    }
    if (p.credits <= beforeCredits) {
      std::cerr << "[test_missions] expected credits to increase on smuggle completion\n";
      ++fails;
    }
    if (p.cargo[(std::size_t)econ::CommodityId::Luxury] > 1e-6) {
      std::cerr << "[test_missions] expected smuggle cargo to be fully consumed on completion\n";
      ++fails;
    }

    const double invAfterComplete = econState.inventory[(int)econ::CommodityId::Luxury];
    if (std::abs(invAfterComplete - invBeforeComplete) > 1e-6) {
      std::cerr << "[test_missions] smuggle completion unexpectedly changed legal inventory\n";
      ++fails;
    }
  }

  
  // Salvage missions: require visiting the mission site (m.scanned) before completing at the station.
  {
    SaveGame p = s;
    p.credits = 0.0;
    setReputation(p, station.factionId, 0.0);
    p.cargo[(std::size_t)econ::CommodityId::Electronics] = 5.0;

    Mission m{};
    m.id = 200;
    m.type = MissionType::Salvage;
    m.factionId = station.factionId;
    m.toSystem = sys.stub.id;
    m.toStation = station.id;
    m.commodity = econ::CommodityId::Electronics;
    m.units = 3.0;
    m.reward = 777.0;
    m.deadlineDay = p.timeDays + 1.0;
    m.scanned = false;
    p.missions = {m};

    const double creditsBefore = p.credits;
    const double cargoBefore = p.cargo[(std::size_t)econ::CommodityId::Electronics];

    const auto res0 = tryCompleteMissionsAtDock(u, sys, station, p.timeDays, p);
    if (res0.completed != 0 || p.missions.front().completed) {
      std::cerr << "[test_missions] expected salvage not to complete before site visit\n";
      ++fails;
    }

    // Mark the salvage site as visited and try again.
    p.missions.front().scanned = true;
    const auto res1 = tryCompleteMissionsAtDock(u, sys, station, p.timeDays, p);
    if (res1.completed != 1 || !p.missions.front().completed) {
      std::cerr << "[test_missions] expected salvage completion after site visit\n";
      ++fails;
    }
    if (p.credits <= creditsBefore) {
      std::cerr << "[test_missions] expected credits to increase on salvage completion\n";
      ++fails;
    }
    if (p.cargo[(std::size_t)econ::CommodityId::Electronics] >= cargoBefore - 1e-6) {
      std::cerr << "[test_missions] expected cargo to be consumed on salvage completion\n";
      ++fails;
    }
  }


// Passenger missions: seat capacity checks + completion path.
  {
    SaveGame p = s;
    p.passengerSeats = 2;
    setReputation(p, station.factionId, 5.0);
    const double rep = getReputation(p, station.factionId);

    Mission offer{};
    offer.type = MissionType::Passenger;
    offer.factionId = station.factionId;
    offer.fromSystem = sys.stub.id;
    offer.fromStation = station.id;
    offer.toSystem = sys.stub.id;
    offer.toStation = station.id;
    offer.units = 3.0; // party size
    offer.reward = 777.0;
    offer.deadlineDay = p.timeDays + 5.0;
    p.missionOffers = {offer};

    std::string err;
    if (acceptMissionOffer(u, station, p.timeDays, rep, p, 0, &err)) {
      std::cerr << "[test_missions] expected passenger accept to fail (not enough seats)\n";
      ++fails;
    }

    // Now accept a valid party size.
    offer.units = 2.0;
    p.missionOffers = {offer};
    if (!acceptMissionOffer(u, station, p.timeDays, rep, p, 0, &err)) {
      std::cerr << "[test_missions] passenger acceptMissionOffer failed: " << err << "\n";
      ++fails;
    } else {
      const double before = p.credits;
      const auto res2 = tryCompleteMissionsAtDock(u, sys, station, p.timeDays, p);
      if (res2.completed != 1) {
        std::cerr << "[test_missions] expected passenger completion count 1\n";
        ++fails;
      }
      if (p.credits <= before) {
        std::cerr << "[test_missions] expected credits to increase on passenger completion\n";
        ++fails;
      }
    }
  }



  // Deadline ticking: missions should fail when timeDays passes deadline (and apply rep penalty).
  {
    SaveGame p = s;
    setReputation(p, station.factionId, 0.0);

    Mission m{};
    m.id = 42;
    m.type = MissionType::Courier;
    m.factionId = station.factionId;
    m.toSystem = sys.stub.id;
    m.toStation = station.id;
    m.reward = 100.0;
    m.deadlineDay = p.timeDays - 0.10; // already overdue
    p.missions = {m};

    const double repBefore = getReputation(p, station.factionId);
    const auto res = tickMissionDeadlines(p, p.timeDays, -4.0);

    if (res.failed != 1) {
      std::cerr << "[test_missions] expected deadline failure count 1\n";
      ++fails;
    }
    if (!p.missions.front().failed) {
      std::cerr << "[test_missions] expected mission to be marked failed on deadline\n";
      ++fails;
    }
    const double repAfter = getReputation(p, station.factionId);
    if (repAfter > repBefore - 3.5) { // allow small float fuzz, but should be ~-4
      std::cerr << "[test_missions] expected rep to decrease on deadline\n";
      ++fails;
    }
  }

  // Multi-hop deliveries: progress at via station, then complete at final station.
  {
    const SystemStub* destStub = nullptr;
    for (const auto& ss : stubs) {
      if (ss.id != sys.stub.id && ss.stationCount > 0) { destStub = &ss; break; }
    }

    if (!destStub) {
      std::cerr << "[test_missions] no second system with stations to test multi-hop\n";
      ++fails;
    } else {
      const auto& destSys = u.getSystem(destStub->id, destStub);
      const auto& destStation = destSys.stations.front();

      SaveGame p = s;
      p.credits = 0.0;
      setReputation(p, station.factionId, 0.0);
      p.cargo[(std::size_t)econ::CommodityId::Food] = 10.0;

      Mission m{};
      m.id = 99;
      m.type = MissionType::MultiDelivery;
      m.factionId = station.factionId;
      m.fromSystem = sys.stub.id;
      m.fromStation = station.id;
      m.viaSystem = sys.stub.id;
      m.viaStation = station.id;
      m.toSystem = destSys.stub.id;
      m.toStation = destStation.id;
      m.commodity = econ::CommodityId::Food;
      m.units = 5.0;
      m.reward = 999.0;
      m.deadlineDay = p.timeDays + 5.0;
      m.cargoProvided = true;
      p.missions = {m};

      const auto viaRes = tryCompleteMissionsAtDock(u, sys, station, p.timeDays, p);
      if (viaRes.progressedMultiLeg != 1) {
        std::cerr << "[test_missions] expected progressedMultiLeg == 1\n";
        ++fails;
      }
      if (p.missions.front().leg != 1 || p.missions.front().completed) {
        std::cerr << "[test_missions] expected multi-hop to advance to leg 1 without completing\n";
        ++fails;
      }

      const double creditsBefore = p.credits;
      const double foodBefore = p.cargo[(std::size_t)econ::CommodityId::Food];

      const auto finRes = tryCompleteMissionsAtDock(u, destSys, destStation, p.timeDays, p);
      if (finRes.completed != 1 || !p.missions.front().completed) {
        std::cerr << "[test_missions] expected multi-hop completion at final station\n";
        ++fails;
      }
      if (p.credits <= creditsBefore) {
        std::cerr << "[test_missions] expected credits to increase on multi-hop completion\n";
        ++fails;
      }
      if (p.cargo[(std::size_t)econ::CommodityId::Food] >= foodBefore - 1e-6) {
        std::cerr << "[test_missions] expected cargo to be consumed on multi-hop completion\n";
        ++fails;
      }
    }
  }

  // Bounty missions: complete on scan or kill events (shared logic).
  {
    SaveGame p = s;
    p.credits = 0.0;
    setReputation(p, station.factionId, 0.0);

    Mission m{};
    m.id = 123;
    m.type = MissionType::BountyScan;
    m.factionId = station.factionId;
    m.toSystem = sys.stub.id;
    m.targetNpcId = 555;
    m.reward = 321.0;
    m.deadlineDay = p.timeDays + 2.0;
    p.missions = {m};

    const auto res = tryCompleteBountyScan(p, sys.stub.id, 555, +2.0);
    if (res.completed != 1) {
      std::cerr << "[test_missions] expected bounty scan completion count 1\n";
      ++fails;
    }
    if (!p.missions.front().completed || !p.missions.front().scanned) {
      std::cerr << "[test_missions] expected bounty scan to mark mission completed+scanned\n";
      ++fails;
    }
    if (p.credits < 321.0 - 1e-6) {
      std::cerr << "[test_missions] expected credits to increase on bounty scan completion\n";
      ++fails;
    }
    if (getReputation(p, station.factionId) < 1.5) {
      std::cerr << "[test_missions] expected rep to increase on bounty scan completion\n";
      ++fails;
    }
  }

  {
    SaveGame p = s;
    p.credits = 0.0;
    setReputation(p, station.factionId, 0.0);

    Mission m{};
    m.id = 124;
    m.type = MissionType::BountyKill;
    m.factionId = station.factionId;
    m.toSystem = sys.stub.id;
    m.targetNpcId = 777;
    m.reward = 500.0;
    m.deadlineDay = p.timeDays + 2.0;
    p.missions = {m};

    const auto res = tryCompleteBountyKill(p, sys.stub.id, 777, +2.0);
    if (res.completed != 1) {
      std::cerr << "[test_missions] expected bounty kill completion count 1\n";
      ++fails;
    }
    if (!p.missions.front().completed) {
      std::cerr << "[test_missions] expected bounty kill to mark mission completed\n";
      ++fails;
    }
    if (p.credits < 500.0 - 1e-6) {
      std::cerr << "[test_missions] expected credits to increase on bounty kill completion\n";
      ++fails;
    }
    if (getReputation(p, station.factionId) < 1.5) {
      std::cerr << "[test_missions] expected rep to increase on bounty kill completion\n";
      ++fails;
    }
  }


  if (fails == 0) std::cout << "[test_missions] pass\n";
  return fails;
}
