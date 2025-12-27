#pragma once

#include "stellar/econ/Market.h"

#include <vector>

namespace stellar::econ {

struct RouteOpportunity {
  CommodityId commodity{};

  // Raw (no fees): toBid - fromAsk.
  double profitPerUnit{0.0};

  // Raw prices (no fees):
  //  - buyPrice:  source station ask (station sells to player)
  //  - sellPrice: destination station bid (station buys from player)
  double buyPrice{0.0};
  double sellPrice{0.0};

  // Feasibility / sizing helpers.
  // These are populated by the route planner so UI/gameplay can avoid
  // recommending trades that the stations (or cargo hold) cannot support.
  double unitsFrom{0.0};      // available at source station (inventory units)
  double unitsToSpace{0.0};   // free capacity at destination station (units)
  double unitsPossible{0.0};  // min(unitsFrom, unitsToSpace[, cargo cap])

  // Commodity sizing.
  double unitMassKg{0.0};
  double profitPerKg{0.0};

  // Total raw profit for `unitsPossible` (no fees).
  double profitTotal{0.0};

  // Optional fees (only filled by bestRoutesForCargo).
  double feeFrom{0.0};
  double feeTo{0.0};
  double netProfitPerUnit{0.0};
  double netProfitTotal{0.0};
};

std::vector<RouteOpportunity> bestRoutes(const StationEconomyState& fromState,
                                        const StationEconomyModel& fromModel,
                                        const StationEconomyState& toState,
                                        const StationEconomyModel& toModel,
                                        double bidAskSpread = 0.10,
                                        std::size_t maxResults = 5);

// Cargo-aware route planner.
//
// Computes feasible trade opportunities given:
//  - station inventories/capacities
//  - an optional cargo mass limit (kg)
//  - optional station fee rates (0..1)
//
// Results are sorted by netProfitTotal (trip profit) when fees/cargo are provided.
std::vector<RouteOpportunity> bestRoutesForCargo(const StationEconomyState& fromState,
                                                const StationEconomyModel& fromModel,
                                                const StationEconomyState& toState,
                                                const StationEconomyModel& toModel,
                                                double cargoCapacityKg,
                                                double fromFeeRate = 0.0,
                                                double toFeeRate = 0.0,
                                                double bidAskSpread = 0.10,
                                                std::size_t maxResults = 5);

} // namespace stellar::econ
