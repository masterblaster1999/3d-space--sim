#pragma once

#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace stellar::sim {

// Optional diagnostics about a planned route.
struct RoutePlanStats {
  int expansions{0};      // nodes popped from the open set
  int visited{0};         // nodes closed/expanded
  int hops{0};            // route.size() - 1
  double distanceLy{0.0}; // sum of straight-line leg distances
  double cost{0.0};       // total A* path cost for the chosen cost model
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

// A* route planner with a simple weighted per-leg cost model:
//
//   legCost = costPerJump + costPerLy * legDistanceLy
//
// This allows different "navigation feels" without changing the graph:
//
//   - Min-distance:  costPerJump = 0,         costPerLy = 1
//   - Fuel-like:     costPerJump = fuelBase,  costPerLy = fuelPerLy
//   - Hybrid:        increase costPerJump to prefer fewer hops
//
// Heuristic: ceil(remainingDistance/maxJumpLy) * costPerJump + remainingDistance * costPerLy
// This is admissible as long as costs are non-negative.
//
// Returns a list of SystemId including start and goal, or empty if no route found.
std::vector<SystemId> plotRouteAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        RoutePlanStats* outStats = nullptr,
                                        std::size_t maxExpansions = 250000);

// A* route planner with an additional travel-risk penalty.
//
// Cost model:
//   baseLeg = costPerJump + costPerLy * d
//   riskLeg = riskWeightPerLy * avgRisk01 * d
//   legCost = baseLeg + riskLeg
//
// where avgRisk01 is 0.5*(risk[a] + risk[b]) for the segment endpoints and
// `risk01PerNode[i]` corresponds to nodes[i].
//
// Notes:
//  - riskWeightPerLy is a *cost per ly* scaling factor; higher values bias toward
//    safer routes (potentially longer / more hops).
//  - The heuristic ignores risk (still admissible) which keeps behavior stable
//    even if the risk model changes.
std::vector<SystemId> plotRouteAStarCostRisk(const std::vector<SystemStub>& nodes,
                                            SystemId startId,
                                            SystemId goalId,
                                            double maxJumpLy,
                                            double costPerJump,
                                            double costPerLy,
                                            double riskWeightPerLy,
                                            std::span<const double> risk01PerNode,
                                            RoutePlanStats* outStats = nullptr,
                                            std::size_t maxExpansions = 250000);


// A* route planner with additional hazard penalty.
// Hazard cost is computed as: hazardWeightPerLy * avgNavDisruption01 * segmentDistanceLy
// where avgNavDisruption01 is sampled along the segment using the galaxy hazard field.
std::vector<SystemId> plotRouteAStarCostHazards(const std::vector<SystemStub>& nodes,
                                               SystemId startId,
                                               SystemId goalId,
                                               double maxJumpLy,
                                               double costPerJump,
                                               double costPerLy,
                                               double hazardWeightPerLy,
                                               core::u64 universeSeed,
                                               double timeDays,
                                               RoutePlanStats* outStats = nullptr,
                                               std::size_t maxExpansions = 250000);


// -----------------------------------------------------------------------------
// Hard route constraints
// -----------------------------------------------------------------------------
//
// These helpers let callers exclude specific systems and/or specific jump links
// from routing. This is useful for:
//   - avoiding hostile or locked-down systems
//   - plotting around dynamic hazards / interdictions
//   - "no-fly" corridors (ban edges)
//
// IMPORTANT:
//  - Constraints are *hard* exclusions (not soft penalties).
//  - `bannedEdges` are treated as *undirected*: {a,b} blocks both a->b and b->a.
//  - The heuristic remains admissible because constraints only remove options.
//
// A jump edge to disallow.
struct RouteEdge {
  SystemId from{0};
  SystemId to{0};

  bool operator==(const RouteEdge& o) const { return from == o.from && to == o.to; }
};

// Constrained variants of the A* route planners.
std::vector<SystemId> plotRouteAStarCostConstrained(const std::vector<SystemStub>& nodes,
                                                   SystemId startId,
                                                   SystemId goalId,
                                                   double maxJumpLy,
                                                   double costPerJump,
                                                   double costPerLy,
                                                   std::span<const SystemId> bannedNodes,
                                                   std::span<const RouteEdge> bannedEdges,
                                                   RoutePlanStats* outStats = nullptr,
                                                   std::size_t maxExpansions = 250000);

std::vector<SystemId> plotRouteAStarCostRiskConstrained(const std::vector<SystemStub>& nodes,
                                                       SystemId startId,
                                                       SystemId goalId,
                                                       double maxJumpLy,
                                                       double costPerJump,
                                                       double costPerLy,
                                                       double riskWeightPerLy,
                                                       std::span<const double> risk01PerNode,
                                                       std::span<const SystemId> bannedNodes,
                                                       std::span<const RouteEdge> bannedEdges,
                                                       RoutePlanStats* outStats = nullptr,
                                                       std::size_t maxExpansions = 250000);

std::vector<SystemId> plotRouteAStarCostHazardsConstrained(const std::vector<SystemStub>& nodes,
                                                          SystemId startId,
                                                          SystemId goalId,
                                                          double maxJumpLy,
                                                          double costPerJump,
                                                          double costPerLy,
                                                          double hazardWeightPerLy,
                                                          core::u64 universeSeed,
                                                          double timeDays,
                                                          std::span<const SystemId> bannedNodes,
                                                          std::span<const RouteEdge> bannedEdges,
                                                          RoutePlanStats* outStats = nullptr,
                                                          std::size_t maxExpansions = 250000);

// Result for K-shortest route planning.
struct KRoute {
  std::vector<SystemId> path;
  int hops{0};
  double distanceLy{0.0};
  double cost{0.0};
};

// K-shortest *loopless* routes between startId and goalId using Yen's algorithm.
//
// - Each route is a list of SystemId including start and goal.
// - Results are ordered by increasing total cost, then deterministically by the
//   path's SystemId sequence.
// - For k=1, this behaves similarly to plotRouteAStarCost (but without returning
//   detailed per-solve expansion diagnostics).
std::vector<KRoute> plotKRoutesAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve = 250000);

// K-shortest routes variant using the risk-augmented leg cost model.
// See plotRouteAStarCostRisk() for the cost definition.
std::vector<KRoute> plotKRoutesAStarCostRisk(const std::vector<SystemStub>& nodes,
                                            SystemId startId,
                                            SystemId goalId,
                                            double maxJumpLy,
                                            double costPerJump,
                                            double costPerLy,
                                            double riskWeightPerLy,
                                            std::span<const double> risk01PerNode,
                                            std::size_t k,
                                            std::size_t maxExpansionsPerSolve = 250000);


// K-shortest routes variant using the hazard-augmented leg cost model.
// See plotRouteAStarCostHazards() for the cost definition.
std::vector<KRoute> plotKRoutesAStarCostHazards(const std::vector<SystemStub>& nodes,
                                                SystemId startId,
                                                SystemId goalId,
                                                double maxJumpLy,
                                                double costPerJump,
                                                double costPerLy,
                                                double hazardWeightPerLy,
                                                core::u64 universeSeed,
                                                double timeDays,
                                                std::size_t k,
                                                std::size_t maxExpansionsPerSolve = 250000);

// Convenience wrapper for hop-minimizing K-shortest paths.
std::vector<KRoute> plotKRoutesAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve = 250000);

// -----------------------------------------------------------------------------
// Diversified K-route selection (MMR-style)
// -----------------------------------------------------------------------------
//
// K-shortest routes are often very similar (differing by only a hop or two).
// For UI/AI purposes it is useful to present alternate routes that are
// meaningfully different.
//
// This helper selects up to k routes from a candidate set by trading off
// cost optimality vs path diversity.
//
// Diversity is measured as:
//   diversity = 1 - max_j routeNodeJaccard01(candidate, selected[j])
//
// `diversityWeight01` in [0,1]:
//   0 -> purely cost-ordered selection (first k unique routes)
//   1 -> greedily maximize diversity vs already-selected routes
//
// Determinism:
//  - ties are broken by (lower cost, then lexicographic SystemId path order).
std::vector<KRoute> selectDiversifiedKRoutesMMR(std::span<const KRoute> candidates,
                                                std::size_t k,
                                                double diversityWeight01 = 0.65,
                                                bool ignoreEndpoints = true);

// Convenience: compute a larger candidate set (kCandidates) with Yen's
// algorithm and then diversify down to k.
//
// If kCandidates==0, it defaults to max(k, k*6).
std::vector<KRoute> plotKRoutesAStarCostDiversified(const std::vector<SystemStub>& nodes,
                                                   SystemId startId,
                                                   SystemId goalId,
                                                   double maxJumpLy,
                                                   double costPerJump,
                                                   double costPerLy,
                                                   std::size_t k,
                                                   double diversityWeight01 = 0.65,
                                                   std::size_t kCandidates = 0,
                                                   bool ignoreEndpoints = true,
                                                   std::size_t maxExpansionsPerSolve = 250000);

// Helper: total straight-line length of a route in ly.
// Returns 0 for empty/single-node routes.
double routeDistanceLy(const std::vector<SystemStub>& nodes,
                       const std::vector<SystemId>& route);

// Helper: total cost of a route using the same cost model as plotRouteAStarCost.
// Returns 0 for empty/single-node routes.
double routeCost(const std::vector<SystemStub>& nodes,
                 const std::vector<SystemId>& route,
                 double costPerJump,
                 double costPerLy);

// Helper: total cost of a route using the risk-augmented model.
// Returns 0 for empty/single-node routes.
double routeCostRisk(const std::vector<SystemStub>& nodes,
                     const std::vector<SystemId>& route,
                     double costPerJump,
                     double costPerLy,
                     double riskWeightPerLy,
                     std::span<const double> risk01PerNode);


// Helper: total cost of a route using the hazard-augmented model.
//
// Hazard cost is computed as:
//   hazardWeightPerLy * avgNavDisruption01 * segmentDistanceLy
// where avgNavDisruption01 is sampled along the segment using the galaxy hazard field.
//
// Returns 0 for empty/single-node routes.
double routeCostHazards(const std::vector<SystemStub>& nodes,
                        const std::vector<SystemId>& route,
                        double costPerJump,
                        double costPerLy,
                        double hazardWeightPerLy,
                        core::u64 universeSeed,
                        double timeDays);


// Helper: node-set similarity (Jaccard) between two routes.
//
// If ignoreEndpoints is true, the first and last SystemId of each route are
// excluded from the comparison. This makes it easier to compare alternate routes
// to the same destination without start/goal dominating the similarity.
// Routes with no internal nodes under the chosen endpoint rule yield an empty
// node set; two empty node sets compare as identical (returns 1).
//
// Returns a value in [0,1].
//  - 1: identical node sets (under the chosen endpoint rule)
//  - 0: no shared nodes
double routeNodeJaccard01(const std::vector<SystemId>& a,
                          const std::vector<SystemId>& b,
                          bool ignoreEndpoints = true);

// Helper: validate that a route is contiguous (each hop <= maxJumpLy) and that all ids exist.
// Useful for tests/tooling.
bool validateRoute(const std::vector<SystemStub>& nodes,
                   const std::vector<SystemId>& route,
                   double maxJumpLy,
                   std::string* outError = nullptr);

} // namespace stellar::sim
