#pragma once

#include "stellar/sim/NavRoute.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Route resilience / detour analysis (replacement paths, edge failures)
// -----------------------------------------------------------------------------
//
// Given a base route from start->goal, this helper computes a detour route for
// each jump edge of the base route, assuming that edge is unavailable.
//
// This is useful for:
//  - highlighting choke points in a navigation corridor
//  - estimating "how bad" an interdiction/closure would be
//  - AI planning that prefers routes with easy detours
//
// Implementation note:
//  - This is a simple O(L) * A* approach (L = base route length in hops).
//    It is intended for gameplay tooling and small/medium graphs.

struct RouteEdgeDetour {
  RouteEdge blocked{};

  bool ok{false};
  std::vector<SystemId> path;

  int hops{0};
  double distanceLy{0.0};
  double cost{0.0};

  // cost / baseCost when baseCost > 0, otherwise +inf if cost>0 or 1 if both 0.
  double costFactor{0.0};
};

struct RouteResilienceReport {
  bool ok{false};

  std::vector<SystemId> basePath;
  RoutePlanStats baseStats{};

  int baseHops{0};
  double baseDistanceLy{0.0};
  double baseCost{0.0};

  // One entry per edge in the base route (same order as basePath edges).
  std::vector<RouteEdgeDetour> edgeDetours;
};

// Analyze detours for each base-route edge using the basic cost model
// (costPerJump + costPerLy * d).
//
// Global constraints (bannedNodes/bannedEdges) are applied to both the base route
// and each detour solve.
RouteResilienceReport analyzeRouteEdgeDetoursCost(const std::vector<SystemStub>& nodes,
                                                 SystemId startId,
                                                 SystemId goalId,
                                                 double maxJumpLy,
                                                 double costPerJump,
                                                 double costPerLy,
                                                 std::span<const SystemId> bannedNodes = {},
                                                 std::span<const RouteEdge> bannedEdges = {},
                                                 std::size_t maxExpansionsPerSolve = 250000);

inline RouteResilienceReport analyzeRouteEdgeDetoursCost(const std::vector<SystemStub>& nodes,
                                                        SystemId startId,
                                                        SystemId goalId,
                                                        double maxJumpLy,
                                                        double costPerJump,
                                                        double costPerLy,
                                                        std::span<const SystemId> bannedNodes,
                                                        std::span<const RouteEdge> bannedEdges,
                                                        std::size_t maxExpansionsPerSolve) {
  RouteResilienceReport rep{};

  RoutePlanStats baseStats{};
  rep.basePath = plotRouteAStarCostConstrained(nodes,
                                              startId,
                                              goalId,
                                              maxJumpLy,
                                              costPerJump,
                                              costPerLy,
                                              bannedNodes,
                                              bannedEdges,
                                              &baseStats,
                                              maxExpansionsPerSolve);
  rep.baseStats = baseStats;
  rep.ok = !rep.basePath.empty();

  if (!rep.ok) {
    return rep;
  }

  rep.baseHops = (int)rep.basePath.size() - 1;
  rep.baseDistanceLy = routeDistanceLy(nodes, rep.basePath);
  rep.baseCost = routeCost(nodes, rep.basePath, costPerJump, costPerLy);

  const double baseCost = rep.baseCost;
  const std::size_t edgeCount = rep.basePath.size() > 0 ? rep.basePath.size() - 1 : 0;
  rep.edgeDetours.reserve(edgeCount);

  std::vector<RouteEdge> combinedEdges;
  combinedEdges.reserve(bannedEdges.size() + 1);

  for (std::size_t i = 0; i + 1 < rep.basePath.size(); ++i) {
    const SystemId a = rep.basePath[i];
    const SystemId b = rep.basePath[i + 1];

    RouteEdgeDetour d{};
    d.blocked = RouteEdge{a, b};

    combinedEdges.clear();
    for (const RouteEdge e : bannedEdges) {
      combinedEdges.push_back(e);
    }
    combinedEdges.push_back(d.blocked);

    RoutePlanStats detStats{};
    d.path = plotRouteAStarCostConstrained(nodes,
                                          startId,
                                          goalId,
                                          maxJumpLy,
                                          costPerJump,
                                          costPerLy,
                                          bannedNodes,
                                          std::span<const RouteEdge>(combinedEdges.data(), combinedEdges.size()),
                                          &detStats,
                                          maxExpansionsPerSolve);

    d.ok = !d.path.empty();

    // Sanity: ensure the detour route does not contain the blocked edge.
    if (d.ok) {
      bool usesBlocked = false;
      for (std::size_t j = 0; j + 1 < d.path.size(); ++j) {
        const SystemId u = d.path[j];
        const SystemId v = d.path[j + 1];
        if ((u == a && v == b) || (u == b && v == a)) {
          usesBlocked = true;
          break;
        }
      }
      if (usesBlocked) {
        d.ok = false;
        d.path.clear();
      }
    }

    if (d.ok) {
      d.hops = (int)d.path.size() - 1;
      d.distanceLy = routeDistanceLy(nodes, d.path);
      d.cost = routeCost(nodes, d.path, costPerJump, costPerLy);

      if (std::isfinite(baseCost) && baseCost > 0.0) {
        d.costFactor = d.cost / baseCost;
      } else if (d.cost <= 0.0) {
        d.costFactor = 1.0;
      } else {
        d.costFactor = std::numeric_limits<double>::infinity();
      }
    } else {
      d.hops = 0;
      d.distanceLy = 0.0;
      d.cost = std::numeric_limits<double>::infinity();
      d.costFactor = std::numeric_limits<double>::infinity();
    }

    rep.edgeDetours.push_back(std::move(d));
  }

  return rep;
}

} // namespace stellar::sim
