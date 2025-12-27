#include "stellar/sim/SaveGame.h"

#include "stellar/econ/Commodity.h"

#include <cmath>
#include <filesystem>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_savegame() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  SaveGame s{};
  s.version = 10;
  s.seed = 123456789ull;
  s.timeDays = 42.5;
  s.currentSystem = 999001;
  s.dockedStation = 1002003;
  s.credits = 9876.5;
  s.shipPosKm = {1.0, 2.0, 3.0};
  s.shipVelKmS = {0.1, 0.0, -0.2};
  s.fuel = 12.3;
  s.fuelMax = 20.0;
  s.cargoCapacityKg = 240.0;
  s.passengerSeats = 6;
  s.smuggleHoldMk = 2;
  s.hull = 0.75;
  s.shield = 0.15;
  s.cargo[static_cast<std::size_t>(econ::CommodityId::Metals)] = 5.5;
  s.cargo[static_cast<std::size_t>(econ::CommodityId::Food)] = 12.0;
  s.scannedKeys = {111, 222, 333};

  // One accepted mission.
  {
    Mission m{};
    m.id = 1;
    m.type = MissionType::Delivery;
    m.factionId = 7;
    m.fromSystem = 10;
    m.fromStation = 11;
    m.toSystem = 20;
    m.toStation = 21;
    m.viaSystem = 15;
    m.viaStation = 16;
    m.leg = 2;
    m.commodity = econ::CommodityId::Food;
    m.units = 8.0;
    m.reward = 1234.0;
    m.deadlineDay = 55.0;
    m.cargoProvided = true;
    m.scanned = true;
    s.missions.push_back(m);
  }

  // Mission offers (board persistence).
  s.missionOffersStationId = 11;
  s.missionOffersDayStamp = 42;
  {
    Mission m{};
    m.id = 2;
    m.type = MissionType::Courier;
    m.factionId = 7;
    m.fromSystem = 10;
    m.fromStation = 11;
    m.toSystem = 12;
    m.toStation = 13;
    m.reward = 250.0;
    m.deadlineDay = 45.0;
    s.missionOffers.push_back(m);
  }

  // Reputation + bounties.
  s.reputation.push_back({7, 12.5});
  s.bounties.push_back({7, 500.0});
  s.bountyVouchers.push_back({7, 200.0});

  // Background traffic day-stamps (per visited system).
  s.trafficStamps.push_back({10, 41});

  // Station economy overrides (cached economies persisted).
  {
    StationEconomyOverride ov{};
    ov.stationId = 9999;
    ov.state.lastUpdateDay = 10.0;
    ov.state.lastSampleDay = 9.5;
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      ov.state.inventory[i] = (double)i * 3.0;
    }
    ov.state.history[0].push_back(econ::PricePoint{1.0, 2.0});
    s.stationOverrides.push_back(ov);
  }

  const std::string path = "savegame_test.sav";

  if (!saveToFile(s, path)) {
    std::cerr << "[test_savegame] saveToFile failed\n";
    return 1;
  }

  SaveGame l{};
  if (!loadFromFile(path, l)) {
    std::cerr << "[test_savegame] loadFromFile failed\n";
    std::filesystem::remove(path);
    return 1;
  }

  // Basic scalar fields.
  if (l.seed != s.seed) {
    std::cerr << "[test_savegame] seed mismatch\n";
    ++fails;
  }
  if (!nearly(l.timeDays, s.timeDays)) {
    std::cerr << "[test_savegame] timeDays mismatch\n";
    ++fails;
  }
  if (l.currentSystem != s.currentSystem || l.dockedStation != s.dockedStation) {
    std::cerr << "[test_savegame] system/station mismatch\n";
    ++fails;
  }
  if (!nearly(l.credits, s.credits)) {
    std::cerr << "[test_savegame] credits mismatch\n";
    ++fails;
  }

  if (l.passengerSeats != s.passengerSeats) {
    std::cerr << "[test_savegame] passengerSeats mismatch\n";
    ++fails;
  }

  if (l.smuggleHoldMk != s.smuggleHoldMk) {
    std::cerr << "[test_savegame] smuggleHoldMk mismatch\n";
    ++fails;
  }

  // Cargo array.
  if (!nearly(l.cargo[static_cast<std::size_t>(econ::CommodityId::Metals)],
              s.cargo[static_cast<std::size_t>(econ::CommodityId::Metals)])) {
    std::cerr << "[test_savegame] cargo[Metals] mismatch\n";
    ++fails;
  }

  // Missions.
  if (l.missions.size() != s.missions.size()) {
    std::cerr << "[test_savegame] missions size mismatch\n";
    ++fails;
  } else {
    const auto& a = l.missions.front();
    const auto& b = s.missions.front();
    if (a.id != b.id || a.type != b.type || a.factionId != b.factionId || a.toStation != b.toStation) {
      std::cerr << "[test_savegame] mission fields mismatch\n";
      ++fails;
    }
    if (a.leg != b.leg || a.scanned != b.scanned) {
      std::cerr << "[test_savegame] mission leg/scanned mismatch\n";
      ++fails;
    }
  }

  // Mission board persistence.
  if (l.missionOffersStationId != s.missionOffersStationId || l.missionOffersDayStamp != s.missionOffersDayStamp) {
    std::cerr << "[test_savegame] mission offers header mismatch\n";
    ++fails;
  }
  if (l.missionOffers.size() != s.missionOffers.size()) {
    std::cerr << "[test_savegame] mission offers size mismatch\n";
    ++fails;
  }

  // Traffic stamps.
  if (l.trafficStamps.size() != s.trafficStamps.size()) {
    std::cerr << "[test_savegame] trafficStamps size mismatch\n";
    ++fails;
  } else if (!l.trafficStamps.empty()) {
    if (l.trafficStamps.front().systemId != s.trafficStamps.front().systemId ||
        l.trafficStamps.front().dayStamp != s.trafficStamps.front().dayStamp) {
      std::cerr << "[test_savegame] trafficStamps entry mismatch\n";
      ++fails;
    }
  }

  // Station overrides.
  if (l.stationOverrides.size() != s.stationOverrides.size()) {
    std::cerr << "[test_savegame] stationOverrides size mismatch\n";
    ++fails;
  } else {
    const auto& a = l.stationOverrides.front();
    const auto& b = s.stationOverrides.front();
    if (a.stationId != b.stationId) {
      std::cerr << "[test_savegame] stationId mismatch\n";
      ++fails;
    }
    if (!nearly(a.state.inventory[0], b.state.inventory[0])) {
      std::cerr << "[test_savegame] override inventory mismatch\n";
      ++fails;
    }
    if (a.state.history[0].empty() || b.state.history[0].empty()) {
      std::cerr << "[test_savegame] override history missing\n";
      ++fails;
    }
  }

  std::filesystem::remove(path);

  if (fails == 0) std::cout << "[test_savegame] pass\n";
  return fails;
}
