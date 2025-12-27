#pragma once

#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <string>
#include <vector>

namespace stellar::sim {

// Optional diagnostics about a planned route.
struct RoutePlanStats {
  int expansions{0};      // nodes popped from the open set
  int visited{0};         // nodes closed/expanded
  int hops{0};            // route.size() - 1
  double distanceLy{0.0}; // sum of straight-line leg distances
  bool reached{false};
};

// A* route planner that minimizes hop count (each jump cost = 1).
// The heuristic is ceil(remainingDistance / maxJumpLy), which is admissible.
//
// NOTE: `nodes` must include stubs for both startId and goalId.
// A common pattern is to pass the vector returned by Universe::queryNearby(...).
//
// Returns a list of SystemId including start and goal, or empty if no route found.
std::vector<SystemId> plotRouteAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        RoutePlanStats* outStats = nullptr,
                                        std::size_t maxExpansions = 250000);

// Helper: total straight-line length of a route in ly.
// Returns 0 for empty/single-node routes.
double routeDistanceLy(const std::vector<SystemStub>& nodes,
                       const std::vector<SystemId>& route);

// Helper: validate that a route is contiguous (each hop <= maxJumpLy) and that all ids exist.
// Useful for tests/tooling.
bool validateRoute(const std::vector<SystemStub>& nodes,
                   const std::vector<SystemId>& route,
                   double maxJumpLy,
                   std::string* outError = nullptr);

} // namespace stellar::sim
