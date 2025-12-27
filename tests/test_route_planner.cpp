#include "stellar/econ/Economy.h"
#include "stellar/econ/RoutePlanner.h"

#include <cmath>
#include <iostream>

int test_route_planner() {
  int fails = 0;

  using namespace stellar::econ;

  StationEconomyModel fromM{};
  StationEconomyModel toM{};

  // A controlled, deterministic toy economy model for testing.
  fromM.priceVolatility = 1.0;
  toM.priceVolatility = 1.0;
  fromM.shockVolatility = 0.0;
  toM.shockVolatility = 0.0;

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    fromM.desiredStock[i] = 100.0;
    fromM.capacity[i] = 300.0;
    toM.desiredStock[i] = 100.0;
    toM.capacity[i] = 300.0;
  }

  StationEconomyState fromS{};
  StationEconomyState toS{};

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    fromS.inventory[i] = 200.0;
    toS.inventory[i] = 200.0;
  }

  const std::size_t iFood = static_cast<std::size_t>(CommodityId::Food);

  // Make FOOD look profitable by making source very cheap and destination ~base.
  fromM.desiredStock[iFood] = 10.0;
  fromM.capacity[iFood] = 300.0;
  fromS.inventory[iFood] = 300.0; // cheap

  // But make destination storage full (unitsToSpace == 0) so the route is NOT feasible.
  toM.desiredStock[iFood] = 1.0;
  toM.capacity[iFood] = 1.0;
  toS.inventory[iFood] = 1.0;

  {
    const auto routes = bestRoutes(fromS, fromM, toS, toM, 0.10, 8);
    if (!routes.empty()) {
      std::cerr << "[test_route_planner] expected no routes when destination has no space\n";
      ++fails;
    }
  }

  // Now give destination space and make it very expensive (scarce) so FOOD is profitable AND feasible.
  toM.desiredStock[iFood] = 100.0;
  toM.capacity[iFood] = 300.0;
  toS.inventory[iFood] = 0.0; // scarce => expensive

  {
    const auto routes = bestRoutes(fromS, fromM, toS, toM, 0.10, 8);

    bool foundFood = false;
    for (const auto& r : routes) {
      if (r.commodity == CommodityId::Food) {
        foundFood = true;

        const double expectedStationUnits = std::min(fromS.inventory[iFood], toM.capacity[iFood] - toS.inventory[iFood]);
        if (std::abs(r.unitsPossible - expectedStationUnits) > 1e-6) {
          std::cerr << "[test_route_planner] unitsPossible mismatch "
                    << r.unitsPossible << " vs " << expectedStationUnits << "\n";
          ++fails;
        }
        if (!(r.profitPerUnit > 0.0)) {
          std::cerr << "[test_route_planner] expected positive profitPerUnit\n";
          ++fails;
        }
        if (!(r.profitTotal > 0.0)) {
          std::cerr << "[test_route_planner] expected positive profitTotal\n";
          ++fails;
        }
        break;
      }
    }

    if (!foundFood) {
      std::cerr << "[test_route_planner] expected FOOD route to appear\n";
      ++fails;
    }
  }

  // Cargo + fee aware planner should cap units by cargo and compute net profit.
  {
    const double cargoKg = 5.0; // Food is 1kg/unit
    const double feeFrom = 0.10;
    const double feeTo = 0.05;

    const auto routes = bestRoutesForCargo(fromS, fromM, toS, toM, cargoKg, feeFrom, feeTo, 0.10, 8);

    bool foundFood = false;
    for (const auto& r : routes) {
      if (r.commodity == CommodityId::Food) {
        foundFood = true;

        if (std::abs(r.unitsPossible - 5.0) > 1e-6) {
          std::cerr << "[test_route_planner] cargo unit cap not applied (unitsPossible="
                    << r.unitsPossible << ")\n";
          ++fails;
        }

        const double expectedNetPerUnit = r.sellPrice * (1.0 - feeTo) - r.buyPrice * (1.0 + feeFrom);
        if (std::abs(r.netProfitPerUnit - expectedNetPerUnit) > 1e-6) {
          std::cerr << "[test_route_planner] netProfitPerUnit mismatch "
                    << r.netProfitPerUnit << " vs " << expectedNetPerUnit << "\n";
          ++fails;
        }

        const double expectedNetTotal = expectedNetPerUnit * r.unitsPossible;
        if (std::abs(r.netProfitTotal - expectedNetTotal) > 1e-6) {
          std::cerr << "[test_route_planner] netProfitTotal mismatch "
                    << r.netProfitTotal << " vs " << expectedNetTotal << "\n";
          ++fails;
        }
        break;
      }
    }

    if (!foundFood) {
      std::cerr << "[test_route_planner] expected FOOD route in cargo-aware results\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_route_planner] pass\n";
  return fails;
}
