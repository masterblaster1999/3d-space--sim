#include "stellar/sim/MissionLogic.h"

#include "stellar/econ/Cargo.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/Reputation.h"
#include "stellar/sim/Universe.h"

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

  if (fails == 0) std::cout << "[test_missions] pass\n";
  return fails;
}
