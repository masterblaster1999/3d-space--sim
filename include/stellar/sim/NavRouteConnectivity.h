#pragma once

#include "stellar/sim/NavRoute.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <span>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// NavRouteConnectivity — edge-disjoint alternate routes + min-cut report
// -----------------------------------------------------------------------------
//
// Motivation:
//   For AI and UI, it's often useful to know if two systems are connected by
//   multiple *independent* corridors. If a blockade, disaster, or interdiction
//   removes a single jump edge, can we still reach the destination?
//
//   The classic graph metric here is the maximum number of *edge-disjoint*
//   paths between start and goal, which equals the value of a unit-capacity
//   max-flow, and also equals the size of a minimum edge cut (max-flow/min-cut).
//
// Design goals:
//   - deterministic (stable results for the same inputs)
//   - zero third-party dependencies
//   - uses the same geometric adjacency rule as NavRoute: an undirected edge
//     exists when distance <= maxJumpLy
//
// Notes:
//   - This is primarily a tooling/analysis helper; it is not meant to replace
//     the A* route planners.

struct EdgeDisjointRoutesResult {
  bool reachable{false};

  // Maximum number of edge-disjoint routes found (capped by maxRoutes in the
  // solver). When reachable==true and maxRoutes was not the limiting factor,
  // this equals minCutEdges.size().
  int maxRoutes{0};

  // Up to maxRoutes edge-disjoint routes (each includes start and goal).
  // Routes are sorted deterministically by their SystemId sequence.
  std::vector<std::vector<SystemId>> routes;

  // One minimum edge cut separating start and goal (undirected edges). This
  // is computed from the final residual graph.
  std::vector<RouteEdge> minCutEdges;

  // SystemIds reachable from start in the final residual graph (source side of
  // the min-cut). Useful for visualization.
  std::vector<SystemId> sourceSide;
};

// Compute up to `maxRoutes` edge-disjoint routes between startId and goalId.
//
// `bannedNodes` / `bannedEdges` behave like the constrained A* functions:
//   - bannedEdges are treated as undirected.
//
// Determinism:
//   - adjacency traversal order is sorted by neighbor SystemId
//   - BFS tie-breaks follow that ordering
inline EdgeDisjointRoutesResult computeEdgeDisjointRoutes(const std::vector<SystemStub>& nodes,
                                                          SystemId startId,
                                                          SystemId goalId,
                                                          double maxJumpLy,
                                                          std::size_t maxRoutes = 16,
                                                          std::span<const SystemId> bannedNodes = {},
                                                          std::span<const RouteEdge> bannedEdges = {}) {
  EdgeDisjointRoutesResult out{};

  if (nodes.empty()) return out;
  if (!std::isfinite(maxJumpLy) || maxJumpLy <= 0.0) return out;

  // ---- Index map -----------------------------------------------------------
  std::unordered_map<SystemId, std::size_t> idToIdx;
  idToIdx.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const auto id = nodes[i].id;
    if (id != 0) idToIdx[id] = i;
  }

  const auto itS = idToIdx.find(startId);
  const auto itT = idToIdx.find(goalId);
  if (itS == idToIdx.end() || itT == idToIdx.end()) return out;

  const std::size_t s = itS->second;
  const std::size_t t = itT->second;

  if (s == t) {
    out.reachable = true;
    out.maxRoutes = 0;
    out.routes.push_back({startId});
    out.sourceSide.push_back(startId);
    return out;
  }

  // ---- Constraint sets -----------------------------------------------------
  std::unordered_set<SystemId> bannedNodeSet;
  bannedNodeSet.reserve(bannedNodes.size());
  for (SystemId id : bannedNodes) {
    if (id != 0) bannedNodeSet.insert(id);
  }
  if (bannedNodeSet.count(startId) || bannedNodeSet.count(goalId)) {
    return out;
  }

  struct EdgeKey {
    SystemId a{0};
    SystemId b{0};
    bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
  };
  struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& k) const {
      std::size_t h = 0;
      // A small, decent hash combine (boost-style).
      auto hc = [&](std::size_t v) {
        h ^= (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
      };
      hc(std::hash<SystemId>()(k.a));
      hc(std::hash<SystemId>()(k.b));
      return h;
    }
  };

  auto canonicalEdgeKey = [](SystemId a, SystemId b) -> EdgeKey {
    if (a == 0 || b == 0) return EdgeKey{};
    if (a < b) return EdgeKey{a, b};
    return EdgeKey{b, a};
  };

  std::unordered_set<EdgeKey, EdgeKeyHash> bannedEdgeSet;
  bannedEdgeSet.reserve(bannedEdges.size());
  for (const RouteEdge& e : bannedEdges) {
    const EdgeKey k = canonicalEdgeKey(e.from, e.to);
    if (k.a != 0) bannedEdgeSet.insert(k);
  }

  auto isNodeAllowed = [&](std::size_t i) -> bool {
    const SystemId id = nodes[i].id;
    if (id == 0) return false;
    return bannedNodeSet.count(id) == 0;
  };

  auto isEdgeAllowed = [&](SystemId a, SystemId b) -> bool {
    const EdgeKey k = canonicalEdgeKey(a, b);
    if (k.a == 0) return false;
    return bannedEdgeSet.count(k) == 0;
  };

  auto systemDistanceLy = [&](std::size_t ia, std::size_t ib) -> double {
    return (nodes[ia].posLy - nodes[ib].posLy).length();
  };

  // ---- Build geometric adjacency (grid accelerated) ------------------------
  // We intentionally mirror NavRoute.cpp's cell scheme for consistency.
  struct CellCoord {
    long long x{0};
    long long y{0};
    long long z{0};
    bool operator==(const CellCoord& o) const { return x == o.x && y == o.y && z == o.z; }
  };
  struct CellHash {
    std::size_t operator()(const CellCoord& c) const {
      std::size_t h = 0;
      auto hc = [&](std::size_t v) {
        h ^= (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
      };
      hc(std::hash<long long>()(c.x));
      hc(std::hash<long long>()(c.y));
      hc(std::hash<long long>()(c.z));
      return h;
    }
  };

  auto cellFor = [&](const math::Vec3d& p, double cellSize) -> CellCoord {
    return CellCoord{(long long)std::floor(p.x / cellSize),
                     (long long)std::floor(p.y / cellSize),
                     (long long)std::floor(p.z / cellSize)};
  };

  auto neighborCells = [&](const CellCoord& c) -> std::array<CellCoord, 27> {
    std::array<CellCoord, 27> outCells{};
    std::size_t n = 0;
    for (long long dz = -1; dz <= 1; ++dz) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dx = -1; dx <= 1; ++dx) {
          outCells[n++] = CellCoord{c.x + dx, c.y + dy, c.z + dz};
        }
      }
    }
    return outCells;
  };

  const double cellSize = maxJumpLy;
  std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash> grid;
  grid.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (nodes[i].id == 0) continue;
    grid[cellFor(nodes[i].posLy, cellSize)].push_back(i);
  }

  struct UEdge {
    std::size_t u{0};
    std::size_t v{0};
    std::size_t other(std::size_t x) const { return (x == u) ? v : u; }
  };

  std::vector<UEdge> edges;
  std::vector<std::vector<int>> adj(nodes.size());

  // Conservative reservation: typical neighbor counts are small.
  edges.reserve(nodes.size() * 8);

  const double maxJump = maxJumpLy;
  const double eps = std::max(1.0e-12, 1.0e-12 * maxJump);

  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (!isNodeAllowed(i)) continue;

    const CellCoord c = cellFor(nodes[i].posLy, cellSize);
    const auto cells = neighborCells(c);
    for (const CellCoord& nc : cells) {
      const auto it = grid.find(nc);
      if (it == grid.end()) continue;

      const auto& bucket = it->second;
      for (std::size_t j : bucket) {
        if (j <= i) continue; // avoid duplicates
        if (!isNodeAllowed(j)) continue;
        if (!isEdgeAllowed(nodes[i].id, nodes[j].id)) continue;

        const double d = systemDistanceLy(i, j);
        if (!(d <= maxJump + eps)) continue;

        const int eIdx = (int)edges.size();
        edges.push_back(UEdge{i, j});
        adj[i].push_back(eIdx);
        adj[j].push_back(eIdx);
      }
    }
  }

  // Deterministic adjacency ordering: sort by neighbor SystemId.
  for (std::size_t i = 0; i < adj.size(); ++i) {
    auto& v = adj[i];
    std::sort(v.begin(), v.end(), [&](int ea, int eb) {
      const std::size_t na = edges[(std::size_t)ea].other(i);
      const std::size_t nb = edges[(std::size_t)eb].other(i);
      const SystemId ida = nodes[na].id;
      const SystemId idb = nodes[nb].id;
      if (ida != idb) return ida < idb;
      return na < nb;
    });
  }

  // Quick reachability check (ignoring capacities).
  auto reachableGraph = [&]() -> bool {
    std::vector<char> vis(nodes.size(), 0);
    std::deque<std::size_t> q;
    q.push_back(s);
    vis[s] = 1;
    while (!q.empty()) {
      const std::size_t x = q.front();
      q.pop_front();
      if (x == t) return true;
      for (int ei : adj[x]) {
        const std::size_t y = edges[(std::size_t)ei].other(x);
        if (!vis[y]) {
          vis[y] = 1;
          q.push_back(y);
        }
      }
    }
    return false;
  };

  if (!reachableGraph()) {
    // Unreachable in the underlying graph.
    return out;
  }
  out.reachable = true;

  // ---- Unit-capacity max-flow on undirected edges (BFS augmentations) ------
  // flow[e] is in {-1,0,+1} relative to the stored orientation u->v.
  //  +1: one unit from u->v
  //  -1: one unit from v->u
  std::vector<int> flow(edges.size(), 0);

  auto residualCanTraverse = [&](std::size_t edgeIdx, int dir) -> bool {
    // dir == +1 means u->v, dir == -1 means v->u.
    const int f = flow[edgeIdx];
    if (dir > 0) return f != +1;
    return f != -1;
  };

  auto augmentOnce = [&]() -> bool {
    std::vector<int> prevNode(nodes.size(), -1);
    std::vector<int> prevEdge(nodes.size(), -1);
    std::vector<int> prevDir(nodes.size(), 0);

    std::deque<std::size_t> q;
    q.push_back(s);
    prevNode[s] = (int)s;

    while (!q.empty() && prevNode[t] == -1) {
      const std::size_t x = q.front();
      q.pop_front();

      for (int ei : adj[x]) {
        const auto& e = edges[(std::size_t)ei];
        const std::size_t y = e.other(x);
        if (prevNode[y] != -1) continue;

        const int dir = (x == e.u) ? +1 : -1;
        if (!residualCanTraverse((std::size_t)ei, dir)) continue;

        prevNode[y] = (int)x;
        prevEdge[y] = ei;
        prevDir[y] = dir;
        q.push_back(y);
        if (y == t) break;
      }
    }

    if (prevNode[t] == -1) return false;

    // Augment by 1 along the discovered path.
    std::size_t cur = t;
    while (cur != s) {
      const int ei = prevEdge[cur];
      const int dir = prevDir[cur];
      if (ei < 0 || dir == 0) break;
      flow[(std::size_t)ei] += dir;
      flow[(std::size_t)ei] = std::clamp(flow[(std::size_t)ei], -1, +1);
      cur = (std::size_t)prevNode[cur];
    }
    return true;
  };

  int maxFlow = 0;
  const std::size_t flowCap = (maxRoutes == 0) ? 0 : maxRoutes;
  while (flowCap == 0 ? false : (std::size_t)maxFlow < flowCap) {
    if (!augmentOnce()) break;
    ++maxFlow;
  }
  out.maxRoutes = maxFlow;

  // ---- Extract explicit routes from the integral flow ----------------------
  // We repeatedly find an s->t path following edges with non-zero flow.
  std::vector<int> flowRem = flow;
  for (int k = 0; k < maxFlow; ++k) {
    std::vector<int> prevNode(nodes.size(), -1);
    std::vector<int> prevEdge(nodes.size(), -1);

    std::deque<std::size_t> q;
    q.push_back(s);
    prevNode[s] = (int)s;

    while (!q.empty() && prevNode[t] == -1) {
      const std::size_t x = q.front();
      q.pop_front();

      for (int ei : adj[x]) {
        const auto& e = edges[(std::size_t)ei];
        const std::size_t y = e.other(x);
        if (prevNode[y] != -1) continue;

        const int dir = (x == e.u) ? +1 : -1;
        const int f = flowRem[(std::size_t)ei];
        const bool hasFlow = (dir > 0) ? (f == +1) : (f == -1);
        if (!hasFlow) continue;

        prevNode[y] = (int)x;
        prevEdge[y] = ei;
        q.push_back(y);
        if (y == t) break;
      }
    }

    if (prevNode[t] == -1) break; // Should not happen, but keep it safe.

    std::vector<SystemId> path;
    std::size_t cur = t;
    while (cur != s) {
      path.push_back(nodes[cur].id);
      const int ei = prevEdge[cur];
      if (ei >= 0) {
        flowRem[(std::size_t)ei] = 0; // consume this edge for subsequent paths
      }
      cur = (std::size_t)prevNode[cur];
    }
    path.push_back(nodes[s].id);
    std::reverse(path.begin(), path.end());
    out.routes.push_back(std::move(path));
  }

  std::sort(out.routes.begin(), out.routes.end());

  // ---- Min-cut from residual reachability ---------------------------------
  std::vector<char> reach(nodes.size(), 0);
  {
    std::deque<std::size_t> q;
    q.push_back(s);
    reach[s] = 1;
    while (!q.empty()) {
      const std::size_t x = q.front();
      q.pop_front();
      for (int ei : adj[x]) {
        const auto& e = edges[(std::size_t)ei];
        const std::size_t y = e.other(x);
        if (reach[y]) continue;

        const int dir = (x == e.u) ? +1 : -1;
        if (!residualCanTraverse((std::size_t)ei, dir)) continue;

        reach[y] = 1;
        q.push_back(y);
      }
    }
  }

  out.sourceSide.clear();
  out.sourceSide.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (reach[i] && nodes[i].id != 0) out.sourceSide.push_back(nodes[i].id);
  }
  std::sort(out.sourceSide.begin(), out.sourceSide.end());

  out.minCutEdges.clear();
  for (std::size_t ei = 0; ei < edges.size(); ++ei) {
    const auto& e = edges[ei];
    const bool a = reach[e.u] != 0;
    const bool b = reach[e.v] != 0;
    if (a == b) continue;

    const SystemId ida = nodes[e.u].id;
    const SystemId idb = nodes[e.v].id;
    if (ida == 0 || idb == 0) continue;
    RouteEdge re;
    if (ida < idb) {
      re.from = ida;
      re.to = idb;
    } else {
      re.from = idb;
      re.to = ida;
    }
    out.minCutEdges.push_back(re);
  }

  std::sort(out.minCutEdges.begin(), out.minCutEdges.end(), [](const RouteEdge& a, const RouteEdge& b) {
    if (a.from != b.from) return a.from < b.from;
    return a.to < b.to;
  });

  // Keep maxRoutes consistent with extracted routes if path extraction was
  // unable to realize all flow units due to unexpected cycles.
  if ((int)out.routes.size() < out.maxRoutes) out.maxRoutes = (int)out.routes.size();

  return out;
}

} // namespace stellar::sim
