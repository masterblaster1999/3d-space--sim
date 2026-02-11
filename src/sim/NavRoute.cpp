#include "stellar/sim/NavRoute.h"
#include "stellar/proc/GalaxyHazards.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace stellar::sim {

namespace {

static double systemDistanceLy(const SystemStub& a, const SystemStub& b) {
  return (a.posLy - b.posLy).length();
}

struct CellCoord {
  long long x{0};
  long long y{0};
  long long z{0};

  bool operator==(const CellCoord& o) const { return x == o.x && y == o.y && z == o.z; }
};

static std::size_t hashCombine(std::size_t h, std::size_t v) {
  // A small, decent hash combine (boost-style).
  return h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
}

struct CellHash {
  std::size_t operator()(const CellCoord& c) const {
    std::size_t h = 0;
    h = hashCombine(h, std::hash<long long>()(c.x));
    h = hashCombine(h, std::hash<long long>()(c.y));
    h = hashCombine(h, std::hash<long long>()(c.z));
    return h;
  }
};

static CellCoord cellFor(const math::Vec3d& p, double cellSize) {
  // cellSize is assumed > 0.
  return CellCoord{
    (long long)std::floor(p.x / cellSize),
    (long long)std::floor(p.y / cellSize),
    (long long)std::floor(p.z / cellSize),
  };
}

static std::unordered_map<SystemId, std::size_t> buildIndex(const std::vector<SystemStub>& nodes) {
  std::unordered_map<SystemId, std::size_t> idx;
  idx.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const auto id = nodes[i].id;
    if (id != 0) idx[id] = i;
  }
  return idx;
}

static void setStats(RoutePlanStats* out, const RoutePlanStats& s) {
  if (out) *out = s;
}

static std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash>
buildGrid(const std::vector<SystemStub>& nodes, double cellSize) {
  std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash> grid;
  grid.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    grid[cellFor(nodes[i].posLy, cellSize)].push_back(i);
  }
  return grid;
}

// --- Back-compat helpers ----------------------------------------------------
//
// A previous WIP iteration of the hazard-aware A* solver used these helper
// names. Keep them as thin wrappers so the algorithm reads clearly.

static std::unordered_map<SystemId, std::size_t> buildIndexMap(const std::vector<SystemStub>& nodes) {
  return buildIndex(nodes);
}

static CellCoord cellForPos(const math::Vec3d& p, double cellSize) {
  return cellFor(p, cellSize);
}

static std::array<CellCoord, 27> neighborsForCell(const CellCoord& c) {
  std::array<CellCoord, 27> out{};
  std::size_t n = 0;
  for (long long dz = -1; dz <= 1; ++dz) {
    for (long long dy = -1; dy <= 1; ++dy) {
      for (long long dx = -1; dx <= 1; ++dx) {
        out[n++] = CellCoord{c.x + dx, c.y + dy, c.z + dz};
      }
    }
  }
  return out;
}

static std::vector<SystemId> rebuildPath(const std::vector<SystemStub>& nodes,
                                        const std::vector<std::size_t>& cameFrom,
                                        std::size_t startIdx,
                                        std::size_t goalIdx,
                                        const std::unordered_map<SystemId, std::size_t>& /*idxMap*/) {
  std::vector<SystemId> path;
  if (nodes.empty()) return path;
  if (startIdx >= nodes.size() || goalIdx >= nodes.size()) return path;
  if (startIdx == goalIdx) {
    path.push_back(nodes[startIdx].id);
    return path;
  }
  if (cameFrom.size() != nodes.size()) return path;

  const std::size_t kNil = (std::size_t)-1;
  std::size_t cur = goalIdx;
  while (cur != kNil && cur != startIdx) {
    path.push_back(nodes[cur].id);
    cur = cameFrom[cur];
  }
  if (cur != startIdx) {
    path.clear();
    return path;
  }
  path.push_back(nodes[startIdx].id);
  std::reverse(path.begin(), path.end());
  return path;
}

struct Edge {
  SystemId from{0};
  SystemId to{0};
  bool operator==(const Edge& o) const { return from == o.from && to == o.to; }
};

struct EdgeHash {
  std::size_t operator()(const Edge& e) const {
    std::size_t h = 0;
    h = hashCombine(h, std::hash<SystemId>()(e.from));
    h = hashCombine(h, std::hash<SystemId>()(e.to));
    return h;
  }
};

static bool pathHasPrefix(const std::vector<SystemId>& path, const std::vector<SystemId>& prefix) {
  if (prefix.size() > path.size()) return false;
  return std::equal(prefix.begin(), prefix.end(), path.begin());
}

static bool pathLexLess(const std::vector<SystemId>& a, const std::vector<SystemId>& b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

static double computeDistanceLy(const std::unordered_map<SystemId, std::size_t>& idx,
                               const std::vector<SystemStub>& nodes,
                               const std::vector<SystemId>& route) {
  if (route.size() < 2) return 0.0;

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    sum += systemDistanceLy(nodes[itA->second], nodes[itB->second]);
  }
  return sum;
}

static double computeCost(const std::unordered_map<SystemId, std::size_t>& idx,
                          const std::vector<SystemStub>& nodes,
                          const std::vector<SystemId>& route,
                          double costPerJump,
                          double costPerLy) {
  if (route.size() < 2) return 0.0;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    const double d = systemDistanceLy(nodes[itA->second], nodes[itB->second]);
    sum += costPerJump + costPerLy * d;
  }
  return sum;
}

static double computeCostRisk(const std::unordered_map<SystemId, std::size_t>& idx,
                              const std::vector<SystemStub>& nodes,
                              const std::vector<SystemId>& route,
                              double costPerJump,
                              double costPerLy,
                              double riskWeightPerLy,
                              std::span<const double> risk01PerNode) {
  if (route.size() < 2) return 0.0;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (riskWeightPerLy < 0.0) riskWeightPerLy = 0.0;

  const bool hasRisk = (risk01PerNode.size() == nodes.size());
  auto riskAt = [&](std::size_t i) -> double {
    if (!hasRisk) return 0.0;
    const double r = risk01PerNode[i];
    if (r < 0.0) return 0.0;
    if (r > 1.0) return 1.0;
    return r;
  };

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    const std::size_t ia = itA->second;
    const std::size_t ib = itB->second;
    const double d = systemDistanceLy(nodes[ia], nodes[ib]);
    const double avgRisk01 = 0.5 * (riskAt(ia) + riskAt(ib));
    sum += costPerJump + costPerLy * d + riskWeightPerLy * avgRisk01 * d;
  }
  return sum;
}


static double computeCostHazards(const std::unordered_map<SystemId, std::size_t>& idx,
                                 const std::vector<SystemStub>& nodes,
                                 const std::vector<SystemId>& route,
                                 double costPerJump,
                                 double costPerLy,
                                 double hazardWeightPerLy,
                                 core::u64 universeSeed,
                                 double timeDays) {
  if (route.size() < 2) return 0.0;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (hazardWeightPerLy < 0.0) hazardWeightPerLy = 0.0;

  // Galaxy hazard sampling parameters: driven by in-game time (drift).
  proc::GalaxyHazardsParams hp{};
  hp.timeDays = timeDays;

  struct EdgeKey {
    std::size_t a{0};
    std::size_t b{0};
    bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
  };
  struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& k) const {
      std::size_t h = std::hash<std::size_t>{}(k.a);
      h = hashCombine(h, std::hash<std::size_t>{}(k.b));
      return h;
    }
  };
  auto edgeKey = [&](std::size_t a, std::size_t b) -> EdgeKey {
    if (a < b) return EdgeKey{a, b};
    return EdgeKey{b, a};
  };

  std::unordered_map<EdgeKey, double, EdgeKeyHash> edgeHazCache;
  edgeHazCache.reserve(route.size() * 2);

  auto navDisrupt01 = [&](std::size_t ia, std::size_t ib) -> double {
    const EdgeKey ek = edgeKey(ia, ib);
    if (const auto it = edgeHazCache.find(ek); it != edgeHazCache.end()) {
      return it->second;
    }
    // Use a small number of midpoint samples; keep consistent with the A* solver.
    const double v = proc::sampleGalaxyNavDisruptionAvgOnSegment(universeSeed,
                                                                 nodes[ia].posLy,
                                                                 nodes[ib].posLy,
                                                                 hp,
                                                                 /*samples=*/5);
    const double clamped = std::clamp(v, 0.0, 1.0);
    edgeHazCache.emplace(ek, clamped);
    return clamped;
  };

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    const std::size_t ia = itA->second;
    const std::size_t ib = itB->second;
    const double d = systemDistanceLy(nodes[ia], nodes[ib]);
    const double nd = navDisrupt01(ia, ib);
    sum += costPerJump + costPerLy * d + hazardWeightPerLy * nd * d;
  }
  return sum;
}

struct AStarSolveResult {
  std::vector<SystemId> path;
  RoutePlanStats stats;
};

static AStarSolveResult aStarCostSolve(const std::vector<SystemStub>& nodes,
                                      const std::unordered_map<SystemId, std::size_t>& idx,
                                      const std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash>& grid,
                                      SystemId startId,
                                      SystemId goalId,
                                      double maxJumpLy,
                                      double costPerJump,
                                      double costPerLy,
                                      const std::unordered_set<SystemId>* bannedNodes,
                                      const std::unordered_set<Edge, EdgeHash>* bannedEdges,
                                      std::size_t maxExpansions) {
  AStarSolveResult out{};
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    out.stats = stats;
    return out;
  }
  if (maxJumpLy <= 0.0) {
    out.stats = stats;
    return out;
  }
  if (nodes.empty()) {
    out.stats = stats;
    return out;
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  if (bannedNodes) {
    if (bannedNodes->count(startId) != 0) {
      out.stats = stats;
      return out;
    }
    if (bannedNodes->count(goalId) != 0) {
      out.stats = stats;
      return out;
    }
  }

  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    out.stats = stats;
    return out;
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    out.stats = stats;
    out.path = {startId};
    return out;
  }

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    const double minHops = std::ceil(d / maxJumpLy);
    return minHops * costPerJump + d * costPerLy;
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = computeDistanceLy(idx, nodes, path);
      stats.cost = gScore[goal];

      out.path = std::move(path);
      out.stats = stats;
      return out;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const SystemId nid = nodes[j].id;
            if (bannedNodes && bannedNodes->count(nid) != 0) continue;
            if (bannedEdges && bannedEdges->count(Edge{nodes[cur.i].id, nid}) != 0) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double legCost = costPerJump + costPerLy * d;
            const double tentative = gScore[cur.i] + legCost;

            if (tentative + 1e-12 < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const double f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  out.stats = stats;
  return out;
}

static AStarSolveResult aStarCostSolveRisk(const std::vector<SystemStub>& nodes,
                                          const std::unordered_map<SystemId, std::size_t>& idx,
                                          const std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash>& grid,
                                          SystemId startId,
                                          SystemId goalId,
                                          double maxJumpLy,
                                          double costPerJump,
                                          double costPerLy,
                                          double riskWeightPerLy,
                                          std::span<const double> risk01PerNode,
                                          const std::unordered_set<SystemId>* bannedNodes,
                                          const std::unordered_set<Edge, EdgeHash>* bannedEdges,
                                          std::size_t maxExpansions) {
  AStarSolveResult out{};
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    out.stats = stats;
    return out;
  }
  if (maxJumpLy <= 0.0) {
    out.stats = stats;
    return out;
  }
  if (nodes.empty()) {
    out.stats = stats;
    return out;
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (riskWeightPerLy < 0.0) riskWeightPerLy = 0.0;

  if (bannedNodes) {
    if (bannedNodes->count(startId) != 0) {
      out.stats = stats;
      return out;
    }
    if (bannedNodes->count(goalId) != 0) {
      out.stats = stats;
      return out;
    }
  }

  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    out.stats = stats;
    return out;
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    out.stats = stats;
    out.path = {startId};
    return out;
  }

  const bool hasRisk = (risk01PerNode.size() == nodes.size());
  auto riskAt = [&](std::size_t i) -> double {
    if (!hasRisk) return 0.0;
    const double r = risk01PerNode[i];
    if (r < 0.0) return 0.0;
    if (r > 1.0) return 1.0;
    return r;
  };

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  // Admissible heuristic: ignore risk.
  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    const double minHops = std::ceil(d / maxJumpLy);
    return minHops * costPerJump + d * costPerLy;
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = computeDistanceLy(idx, nodes, path);
      stats.cost = gScore[goal];

      out.path = std::move(path);
      out.stats = stats;
      return out;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const SystemId nid = nodes[j].id;
            if (bannedNodes && bannedNodes->count(nid) != 0) continue;
            if (bannedEdges && bannedEdges->count(Edge{nodes[cur.i].id, nid}) != 0) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double avgRisk01 = 0.5 * (riskAt(cur.i) + riskAt(j));
            const double legCost = costPerJump + costPerLy * d + riskWeightPerLy * avgRisk01 * d;
            const double tentative = gScore[cur.i] + legCost;

            if (tentative + 1e-12 < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const double f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  out.stats = stats;
  return out;
}

} // namespace


static AStarSolveResult aStarCostSolveHazards(const std::vector<SystemStub>& nodes,
                                             const std::unordered_map<SystemId, std::size_t>& idx,
                                             const std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash>& grid,
                                             SystemId startId,
                                             SystemId goalId,
                                             double maxJumpLy,
                                             double costPerJump,
                                             double costPerLy,
                                             double hazardWeightPerLy,
                                             core::u64 universeSeed,
                                             double timeDays,
                                             const std::unordered_set<SystemId>* bannedNodes,
                                             const std::unordered_set<Edge, EdgeHash>* bannedEdges,
                                             std::size_t maxExpansions) {
  AStarSolveResult out{};
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    out.stats = stats;
    return out;
  }
  if (maxJumpLy <= 0.0) {
    out.stats = stats;
    return out;
  }
  if (nodes.empty()) {
    out.stats = stats;
    return out;
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (hazardWeightPerLy < 0.0) hazardWeightPerLy = 0.0;

  if (bannedNodes) {
    if (bannedNodes->count(startId) != 0) {
      out.stats = stats;
      return out;
    }
    if (bannedNodes->count(goalId) != 0) {
      out.stats = stats;
      return out;
    }
  }

  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    out.stats = stats;
    return out;
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    out.stats = stats;
    out.path = {startId};
    return out;
  }

  // Galaxy hazard sampling parameters: driven by in-game time (drift).
  proc::GalaxyHazardsParams hp{};
  hp.timeDays = timeDays;

  struct EdgeKey {
    std::size_t a{0};
    std::size_t b{0};
    bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
  };
  struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& k) const {
      std::size_t h = std::hash<std::size_t>{}(k.a);
      h = hashCombine(h, std::hash<std::size_t>{}(k.b));
      return h;
    }
  };
  auto edgeKey = [&](std::size_t a, std::size_t b) -> EdgeKey {
    if (a < b) return EdgeKey{a, b};
    return EdgeKey{b, a};
  };

  std::unordered_map<EdgeKey, double, EdgeKeyHash> edgeHazCache;
  edgeHazCache.reserve(4096);

  auto navDisrupt01 = [&](std::size_t ia, std::size_t ib) -> double {
    const EdgeKey ek = edgeKey(ia, ib);
    if (const auto it = edgeHazCache.find(ek); it != edgeHazCache.end()) {
      return it->second;
    }

    // Midpoint sampling to approximate the average hazard along the segment.
    // Keep this small; it's called in the inner loop of A*.
    const double v = proc::sampleGalaxyNavDisruptionAvgOnSegment(universeSeed,
                                                                 nodes[ia].posLy,
                                                                 nodes[ib].posLy,
                                                                 hp,
                                                                 /*samples=*/5);
    const double clamped = std::clamp(v, 0.0, 1.0);
    edgeHazCache.emplace(ek, clamped);
    return clamped;
  };

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  // Admissible heuristic: ignore hazards.
  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    const double minHops = std::ceil(d / maxJumpLy);
    return minHops * costPerJump + d * costPerLy;
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = computeDistanceLy(idx, nodes, path);
      stats.cost = gScore[goal];

      out.path = std::move(path);
      out.stats = stats;
      return out;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const SystemId nid = nodes[j].id;
            if (bannedNodes && bannedNodes->count(nid) != 0) continue;
            if (bannedEdges && bannedEdges->count(Edge{nodes[cur.i].id, nid}) != 0) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double nd = navDisrupt01(cur.i, j);
            const double legCost = costPerJump + costPerLy * d + hazardWeightPerLy * nd * d;
            const double tentative = gScore[cur.i] + legCost;

            if (tentative + 1e-12 < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const double f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  out.stats = stats;
  return out;
}

std::vector<SystemId> plotRouteAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        RoutePlanStats* outStats,
                                        std::size_t maxExpansions) {
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0) {
    setStats(outStats, stats);
    return {};
  }
  if (nodes.empty()) {
    setStats(outStats, stats);
    return {};
  }

  const auto idx = buildIndex(nodes);
  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    setStats(outStats, stats);
    return {startId};
  }

  // Spatial hash grid for neighbor queries.
  // Cell size is maxJumpLy; any reachable neighbor must be in the same or adjacent cell.
  const auto grid = buildGrid(nodes, maxJumpLy);

  std::vector<int> cameFrom(N, -1);
  std::vector<int> gScore(N, std::numeric_limits<int>::max());
  std::vector<char> closed(N, 0);

  auto heuristic = [&](std::size_t i) -> int {
    if (i == goal) return 0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    return (int)std::ceil(d / maxJumpLy);
  };

  struct QN {
    int f{0};
    int g{0};
    std::size_t i{0};
  };

  // Prefer lower f, then lower g, then lower index (stable/deterministic tie-break).
  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0;
  open.push(QN{heuristic(start), 0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = routeDistanceLy(nodes, path);
      stats.cost = (double)stats.hops;

      setStats(outStats, stats);
      return path;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    // Neighbors: any node within maxJumpLy (search the surrounding 3x3x3 cell block).
    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const int tentative = gScore[cur.i] + 1;
            if (tentative < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const int f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  // Not found (or hit expansion cap).
  setStats(outStats, stats);
  return {};
}

std::vector<SystemId> plotRouteAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        RoutePlanStats* outStats,
                                        std::size_t maxExpansions) {
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0) {
    setStats(outStats, stats);
    return {};
  }
  if (nodes.empty()) {
    setStats(outStats, stats);
    return {};
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  // Degenerate: no optimization signal. Fall back to the classic hop planner.
  if (costPerJump <= 0.0 && costPerLy <= 0.0) {
    return plotRouteAStarHops(nodes, startId, goalId, maxJumpLy, outStats, maxExpansions);
  }

  const auto idx = buildIndex(nodes);
  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    setStats(outStats, stats);
    return {startId};
  }

  // Spatial hash grid for neighbor queries.
  const auto grid = buildGrid(nodes, maxJumpLy);

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    const double minHops = std::ceil(d / maxJumpLy);
    return minHops * costPerJump + d * costPerLy;
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = routeDistanceLy(nodes, path);
      stats.cost = gScore[goal];

      setStats(outStats, stats);
      return path;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double legCost = costPerJump + costPerLy * d;
            const double tentative = gScore[cur.i] + legCost;

            if (tentative + 1e-12 < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const double f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  setStats(outStats, stats);
  return {};
}

std::vector<SystemId> plotRouteAStarCostRisk(const std::vector<SystemStub>& nodes,
                                            SystemId startId,
                                            SystemId goalId,
                                            double maxJumpLy,
                                            double costPerJump,
                                            double costPerLy,
                                            double riskWeightPerLy,
                                            std::span<const double> risk01PerNode,
                                            RoutePlanStats* outStats,
                                            std::size_t maxExpansions) {
  RoutePlanStats stats{};

  if (riskWeightPerLy <= 0.0) {
    return plotRouteAStarCost(nodes, startId, goalId, maxJumpLy, costPerJump, costPerLy, outStats, maxExpansions);
  }

  if (startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0) {
    setStats(outStats, stats);
    return {};
  }
  if (nodes.empty()) {
    setStats(outStats, stats);
    return {};
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (riskWeightPerLy < 0.0) riskWeightPerLy = 0.0;

  // Degenerate: no optimization signal. Fall back to hop planner.
  if (costPerJump <= 0.0 && costPerLy <= 0.0 && riskWeightPerLy <= 0.0) {
    return plotRouteAStarHops(nodes, startId, goalId, maxJumpLy, outStats, maxExpansions);
  }

  const auto idx = buildIndex(nodes);
  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    setStats(outStats, stats);
    return {startId};
  }

  const bool hasRisk = (risk01PerNode.size() == nodes.size());
  auto riskAt = [&](std::size_t i) -> double {
    if (!hasRisk) return 0.0;
    const double r = risk01PerNode[i];
    if (r < 0.0) return 0.0;
    if (r > 1.0) return 1.0;
    return r;
  };

  // Spatial hash grid for neighbor queries.
  const auto grid = buildGrid(nodes, maxJumpLy);

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  // Admissible heuristic: ignore risk.
  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    const double minHops = std::ceil(d / maxJumpLy);
    return minHops * costPerJump + d * costPerLy;
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = routeDistanceLy(nodes, path);
      stats.cost = gScore[goal];

      setStats(outStats, stats);
      return path;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double avgRisk01 = 0.5 * (riskAt(cur.i) + riskAt(j));
            const double legCost = costPerJump + costPerLy * d + riskWeightPerLy * avgRisk01 * d;
            const double tentative = gScore[cur.i] + legCost;

            if (tentative + 1e-12 < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const double f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  setStats(outStats, stats);
  return {};
}

std::vector<SystemId> plotRouteAStarCostHazards(const std::vector<SystemStub>& nodes,
                                               SystemId startId,
                                               SystemId goalId,
                                               double maxJumpLy,
                                               double costPerJump,
                                               double costPerLy,
                                               double hazardWeightPerLy,
                                               core::u64 universeSeed,
                                               double timeDays,
                                               RoutePlanStats* outStats,
                                               std::size_t maxExpansions) {
  // If hazard weighting is disabled, fall back to the plain cost model.
  if (hazardWeightPerLy <= 0.0) {
    return plotRouteAStarCost(nodes, startId, goalId, maxJumpLy, costPerJump, costPerLy, outStats, maxExpansions);
  }

  RoutePlanStats stats{};

  if (nodes.empty() || startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0 || maxExpansions == 0) {
    setStats(outStats, stats);
    return {};
  }

  if (hazardWeightPerLy < 0.0) hazardWeightPerLy = 0.0;

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const auto grid = buildGrid(nodes, maxJumpLy);

  const auto res = aStarCostSolveHazards(nodes, idx, grid,
                                        startId, goalId,
                                        maxJumpLy,
                                        costPerJump, costPerLy, hazardWeightPerLy,
                                        universeSeed, timeDays,
                                        nullptr, nullptr,
                                        maxExpansions);

  setStats(outStats, res.stats);
  return res.path;
}



std::vector<SystemId> plotRouteAStarCostConstrained(const std::vector<SystemStub>& nodes,
                                                   SystemId startId,
                                                   SystemId goalId,
                                                   double maxJumpLy,
                                                   double costPerJump,
                                                   double costPerLy,
                                                   std::span<const SystemId> bannedNodes,
                                                   std::span<const RouteEdge> bannedEdges,
                                                   RoutePlanStats* outStats,
                                                   std::size_t maxExpansions) {
  RoutePlanStats stats{};

  if (nodes.empty() || startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0 || maxExpansions == 0) {
    setStats(outStats, stats);
    return {};
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  // Degenerate: treat as hop-minimizing (while still respecting constraints).
  if (costPerJump <= 0.0 && costPerLy <= 0.0) {
    costPerJump = 1.0;
    costPerLy = 0.0;
  }

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const auto grid = buildGrid(nodes, maxJumpLy);

  std::unordered_set<SystemId> bannedNodeSet;
  std::unordered_set<Edge, EdgeHash> bannedEdgeSet;

  const std::unordered_set<SystemId>* bannedNodePtr = nullptr;
  const std::unordered_set<Edge, EdgeHash>* bannedEdgePtr = nullptr;

  if (!bannedNodes.empty()) {
    bannedNodeSet.reserve(bannedNodes.size());
    for (const SystemId id : bannedNodes) {
      if (id != 0) bannedNodeSet.insert(id);
    }
    bannedNodePtr = &bannedNodeSet;
  }

  if (!bannedEdges.empty()) {
    bannedEdgeSet.reserve(bannedEdges.size() * 2);
    for (const RouteEdge e : bannedEdges) {
      if (e.from == 0 || e.to == 0 || e.from == e.to) continue;
      // Treat as undirected by inserting both directions.
      bannedEdgeSet.insert(Edge{e.from, e.to});
      bannedEdgeSet.insert(Edge{e.to, e.from});
    }
    bannedEdgePtr = &bannedEdgeSet;
  }

  const auto res = aStarCostSolve(nodes, idx, grid,
                                 startId, goalId,
                                 maxJumpLy,
                                 costPerJump, costPerLy,
                                 bannedNodePtr, bannedEdgePtr,
                                 maxExpansions);

  setStats(outStats, res.stats);
  return res.path;
}

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
                                                       RoutePlanStats* outStats,
                                                       std::size_t maxExpansions) {
  // If risk weighting is disabled, fall back to the plain constrained cost model.
  if (riskWeightPerLy <= 0.0) {
    return plotRouteAStarCostConstrained(nodes, startId, goalId, maxJumpLy,
                                         costPerJump, costPerLy,
                                         bannedNodes, bannedEdges,
                                         outStats, maxExpansions);
  }

  RoutePlanStats stats{};

  if (nodes.empty() || startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0 || maxExpansions == 0) {
    setStats(outStats, stats);
    return {};
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (riskWeightPerLy < 0.0) riskWeightPerLy = 0.0;

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const auto grid = buildGrid(nodes, maxJumpLy);

  std::unordered_set<SystemId> bannedNodeSet;
  std::unordered_set<Edge, EdgeHash> bannedEdgeSet;

  const std::unordered_set<SystemId>* bannedNodePtr = nullptr;
  const std::unordered_set<Edge, EdgeHash>* bannedEdgePtr = nullptr;

  if (!bannedNodes.empty()) {
    bannedNodeSet.reserve(bannedNodes.size());
    for (const SystemId id : bannedNodes) {
      if (id != 0) bannedNodeSet.insert(id);
    }
    bannedNodePtr = &bannedNodeSet;
  }

  if (!bannedEdges.empty()) {
    bannedEdgeSet.reserve(bannedEdges.size() * 2);
    for (const RouteEdge e : bannedEdges) {
      if (e.from == 0 || e.to == 0 || e.from == e.to) continue;
      bannedEdgeSet.insert(Edge{e.from, e.to});
      bannedEdgeSet.insert(Edge{e.to, e.from});
    }
    bannedEdgePtr = &bannedEdgeSet;
  }

  const auto res = aStarCostSolveRisk(nodes, idx, grid,
                                     startId, goalId,
                                     maxJumpLy,
                                     costPerJump, costPerLy, riskWeightPerLy,
                                     risk01PerNode,
                                     bannedNodePtr, bannedEdgePtr,
                                     maxExpansions);

  setStats(outStats, res.stats);
  return res.path;
}

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
                                                          RoutePlanStats* outStats,
                                                          std::size_t maxExpansions) {
  // If hazard weighting is disabled, fall back to the plain constrained cost model.
  if (hazardWeightPerLy <= 0.0) {
    return plotRouteAStarCostConstrained(nodes, startId, goalId, maxJumpLy,
                                         costPerJump, costPerLy,
                                         bannedNodes, bannedEdges,
                                         outStats, maxExpansions);
  }

  RoutePlanStats stats{};

  if (nodes.empty() || startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0 || maxExpansions == 0) {
    setStats(outStats, stats);
    return {};
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (hazardWeightPerLy < 0.0) hazardWeightPerLy = 0.0;

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const auto grid = buildGrid(nodes, maxJumpLy);

  std::unordered_set<SystemId> bannedNodeSet;
  std::unordered_set<Edge, EdgeHash> bannedEdgeSet;

  const std::unordered_set<SystemId>* bannedNodePtr = nullptr;
  const std::unordered_set<Edge, EdgeHash>* bannedEdgePtr = nullptr;

  if (!bannedNodes.empty()) {
    bannedNodeSet.reserve(bannedNodes.size());
    for (const SystemId id : bannedNodes) {
      if (id != 0) bannedNodeSet.insert(id);
    }
    bannedNodePtr = &bannedNodeSet;
  }

  if (!bannedEdges.empty()) {
    bannedEdgeSet.reserve(bannedEdges.size() * 2);
    for (const RouteEdge e : bannedEdges) {
      if (e.from == 0 || e.to == 0 || e.from == e.to) continue;
      bannedEdgeSet.insert(Edge{e.from, e.to});
      bannedEdgeSet.insert(Edge{e.to, e.from});
    }
    bannedEdgePtr = &bannedEdgeSet;
  }

  const auto res = aStarCostSolveHazards(nodes, idx, grid,
                                        startId, goalId,
                                        maxJumpLy,
                                        costPerJump, costPerLy, hazardWeightPerLy,
                                        universeSeed, timeDays,
                                        bannedNodePtr, bannedEdgePtr,
                                        maxExpansions);

  setStats(outStats, res.stats);
  return res.path;
}


std::vector<KRoute> plotKRoutesAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve) {
  std::vector<KRoute> out;
  if (k == 0) return out;

  if (startId == 0 || goalId == 0) return out;
  if (maxJumpLy <= 0.0) return out;
  if (nodes.empty()) return out;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  // Degenerate: treat as hop-minimizing.
  if (costPerJump <= 0.0 && costPerLy <= 0.0) {
    costPerJump = 1.0;
    costPerLy = 0.0;
  }

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) return out;

  const auto grid = buildGrid(nodes, maxJumpLy);

  // First shortest path.
  {
    const auto base = aStarCostSolve(nodes, idx, grid, startId, goalId, maxJumpLy,
                                    costPerJump, costPerLy,
                                    nullptr, nullptr,
                                    maxExpansionsPerSolve);
    if (base.path.empty()) return out;

    out.push_back(KRoute{base.path, base.stats.hops, base.stats.distanceLy, base.stats.cost});
  }

  struct Candidate {
    std::vector<SystemId> path;
    int hops{0};
    double distanceLy{0.0};
    double cost{0.0};
  };

  const auto containsPath = [](const std::vector<KRoute>& routes, const std::vector<SystemId>& p) {
    for (const auto& r : routes) {
      if (r.path == p) return true;
    }
    return false;
  };

  const auto containsPathCand = [](const std::vector<Candidate>& cands, const std::vector<SystemId>& p) {
    for (const auto& c : cands) {
      if (c.path == p) return true;
    }
    return false;
  };

  const auto betterCand = [](const Candidate& a, const Candidate& b) {
    const double da = a.cost;
    const double db = b.cost;
    if (std::abs(da - db) > 1e-12) return da < db;
    return pathLexLess(a.path, b.path);
  };

  std::vector<Candidate> candidates;

  // Yen's algorithm:
  //  - out[0] is the shortest.
  //  - candidates holds the next-best deviations.
  for (std::size_t kth = 1; kth < k; ++kth) {
    const auto& prev = out[kth - 1].path;
    if (prev.size() < 2) break;

    for (std::size_t i = 0; i + 1 < prev.size(); ++i) {
      const SystemId spurNode = prev[i];

      // Root path includes spur node.
      std::vector<SystemId> root(prev.begin(), prev.begin() + (long long)i + 1);

      // Ban nodes in the root path *except* the spur node to enforce loopless paths.
      std::unordered_set<SystemId> bannedNodes;
      bannedNodes.reserve(i);
      for (std::size_t r = 0; r < i; ++r) bannedNodes.insert(root[r]);

      // Ban edges that would recreate any previously found shortest path with this root.
      std::unordered_set<Edge, EdgeHash> bannedEdges;
      for (const auto& found : out) {
        const auto& p = found.path;
        if (p.size() > i + 1 && pathHasPrefix(p, root)) {
          bannedEdges.insert(Edge{p[i], p[i + 1]});
        }
      }

      const auto spur = aStarCostSolve(nodes, idx, grid,
                                      spurNode, goalId,
                                      maxJumpLy,
                                      costPerJump, costPerLy,
                                      &bannedNodes,
                                      &bannedEdges,
                                      maxExpansionsPerSolve);

      if (spur.path.empty()) continue;

      // Combine root + spur (skip spurNode duplicate).
      std::vector<SystemId> total = root;
      if (spur.path.size() > 1) {
        total.insert(total.end(), spur.path.begin() + 1, spur.path.end());
      }

      if (total.size() < 2) continue;
      if (containsPath(out, total)) continue;
      if (containsPathCand(candidates, total)) continue;

      const double dist = computeDistanceLy(idx, nodes, total);
      const double cst  = computeCost(idx, nodes, total, costPerJump, costPerLy);

      const int hops = (int)total.size() - 1;
      candidates.push_back(Candidate{std::move(total), hops, dist, cst});
    }

    if (candidates.empty()) break;

    // Pick the best candidate.
    std::size_t bestIdx = 0;
    for (std::size_t i = 1; i < candidates.size(); ++i) {
      if (betterCand(candidates[i], candidates[bestIdx])) bestIdx = i;
    }

    Candidate best = std::move(candidates[bestIdx]);
    candidates.erase(candidates.begin() + (long long)bestIdx);

    out.push_back(KRoute{std::move(best.path), best.hops, best.distanceLy, best.cost});
  }

  return out;
}

std::vector<KRoute> plotKRoutesAStarCostRisk(const std::vector<SystemStub>& nodes,
                                            SystemId startId,
                                            SystemId goalId,
                                            double maxJumpLy,
                                            double costPerJump,
                                            double costPerLy,
                                            double riskWeightPerLy,
                                            std::span<const double> risk01PerNode,
                                            std::size_t k,
                                            std::size_t maxExpansionsPerSolve) {
  if (riskWeightPerLy <= 0.0) {
    return plotKRoutesAStarCost(nodes, startId, goalId, maxJumpLy, costPerJump, costPerLy, k, maxExpansionsPerSolve);
  }

  std::vector<KRoute> out;
  if (k == 0) return out;

  if (startId == 0 || goalId == 0) return out;
  if (maxJumpLy <= 0.0) return out;
  if (nodes.empty()) return out;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (riskWeightPerLy < 0.0) riskWeightPerLy = 0.0;

  // Degenerate: if nothing contributes to cost, treat as hop-minimizing.
  if (costPerJump <= 0.0 && costPerLy <= 0.0 && riskWeightPerLy <= 0.0) {
    costPerJump = 1.0;
    costPerLy = 0.0;
  }

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) return out;

  const auto grid = buildGrid(nodes, maxJumpLy);

  // First shortest path.
  {
    const auto base = aStarCostSolveRisk(nodes, idx, grid, startId, goalId, maxJumpLy,
                                        costPerJump, costPerLy, riskWeightPerLy, risk01PerNode,
                                        nullptr, nullptr,
                                        maxExpansionsPerSolve);
    if (base.path.empty()) return out;
    out.push_back(KRoute{base.path, base.stats.hops, base.stats.distanceLy, base.stats.cost});
  }

  struct Candidate {
    std::vector<SystemId> path;
    int hops{0};
    double distanceLy{0.0};
    double cost{0.0};
  };

  const auto containsPath = [](const std::vector<KRoute>& routes, const std::vector<SystemId>& p) {
    for (const auto& r : routes) {
      if (r.path == p) return true;
    }
    return false;
  };

  const auto containsPathCand = [](const std::vector<Candidate>& cands, const std::vector<SystemId>& p) {
    for (const auto& c : cands) {
      if (c.path == p) return true;
    }
    return false;
  };

  const auto betterCand = [](const Candidate& a, const Candidate& b) {
    const double da = a.cost;
    const double db = b.cost;
    if (std::abs(da - db) > 1e-12) return da < db;
    return pathLexLess(a.path, b.path);
  };

  std::vector<Candidate> candidates;

  // Yen's algorithm:
  //  - out[0] is the shortest.
  //  - candidates holds the next-best deviations.
  for (std::size_t kth = 1; kth < k; ++kth) {
    const auto& prev = out[kth - 1].path;
    if (prev.size() < 2) break;

    for (std::size_t i = 0; i + 1 < prev.size(); ++i) {
      const SystemId spurNode = prev[i];

      // Root path includes spur node.
      std::vector<SystemId> root(prev.begin(), prev.begin() + (long long)i + 1);

      // Ban nodes in the root path *except* the spur node to enforce loopless paths.
      std::unordered_set<SystemId> bannedNodes;
      bannedNodes.reserve(i);
      for (std::size_t r = 0; r < i; ++r) bannedNodes.insert(root[r]);

      // Ban edges that would recreate any previously found shortest path with this root.
      std::unordered_set<Edge, EdgeHash> bannedEdges;
      for (const auto& found : out) {
        const auto& p = found.path;
        if (p.size() > i + 1 && pathHasPrefix(p, root)) {
          bannedEdges.insert(Edge{p[i], p[i + 1]});
        }
      }

      const auto spur = aStarCostSolveRisk(nodes, idx, grid,
                                          spurNode, goalId,
                                          maxJumpLy,
                                          costPerJump, costPerLy, riskWeightPerLy, risk01PerNode,
                                          &bannedNodes,
                                          &bannedEdges,
                                          maxExpansionsPerSolve);

      if (spur.path.empty()) continue;

      // Combine root + spur (skip spurNode duplicate).
      std::vector<SystemId> total = root;
      if (spur.path.size() > 1) {
        total.insert(total.end(), spur.path.begin() + 1, spur.path.end());
      }

      if (total.size() < 2) continue;
      if (containsPath(out, total)) continue;
      if (containsPathCand(candidates, total)) continue;

      const double dist = computeDistanceLy(idx, nodes, total);
      const double cst  = computeCostRisk(idx, nodes, total, costPerJump, costPerLy, riskWeightPerLy, risk01PerNode);

      const int hops = (int)total.size() - 1;
      candidates.push_back(Candidate{std::move(total), hops, dist, cst});
    }

    if (candidates.empty()) break;

    // Pick the best candidate.
    std::size_t bestIdx = 0;
    for (std::size_t i = 1; i < candidates.size(); ++i) {
      if (betterCand(candidates[i], candidates[bestIdx])) bestIdx = i;
    }

    Candidate best = std::move(candidates[bestIdx]);
    candidates.erase(candidates.begin() + (long long)bestIdx);

    out.push_back(KRoute{std::move(best.path), best.hops, best.distanceLy, best.cost});
  }

  return out;
}


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
                                                std::size_t maxExpansionsPerSolve) {
  if (hazardWeightPerLy <= 0.0) {
    return plotKRoutesAStarCost(nodes, startId, goalId, maxJumpLy, costPerJump, costPerLy, k, maxExpansionsPerSolve);
  }

  std::vector<KRoute> out;
  if (k == 0) return out;

  if (startId == 0 || goalId == 0) return out;
  if (maxJumpLy <= 0.0) return out;
  if (nodes.empty()) return out;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (hazardWeightPerLy < 0.0) hazardWeightPerLy = 0.0;

  // Degenerate: if nothing contributes to cost, treat as hop-minimizing.
  if (costPerJump <= 0.0 && costPerLy <= 0.0 && hazardWeightPerLy <= 0.0) {
    costPerJump = 1.0;
    costPerLy = 0.0;
  }

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) return out;

  const auto grid = buildGrid(nodes, maxJumpLy);

  // First shortest path.
  {
    const auto base = aStarCostSolveHazards(nodes, idx, grid, startId, goalId, maxJumpLy,
                                           costPerJump, costPerLy, hazardWeightPerLy,
                                           universeSeed, timeDays,
                                           nullptr, nullptr,
                                           maxExpansionsPerSolve);
    if (base.path.empty()) return out;
    out.push_back(KRoute{base.path, base.stats.hops, base.stats.distanceLy, base.stats.cost});
  }

  struct Candidate {
    std::vector<SystemId> path;
    int hops{0};
    double distanceLy{0.0};
    double cost{0.0};
  };

  const auto containsPath = [](const std::vector<KRoute>& routes, const std::vector<SystemId>& p) {
    for (const auto& r : routes) {
      if (r.path == p) return true;
    }
    return false;
  };

  const auto containsPathCand = [](const std::vector<Candidate>& cands, const std::vector<SystemId>& p) {
    for (const auto& c : cands) {
      if (c.path == p) return true;
    }
    return false;
  };

  const auto betterCand = [](const Candidate& a, const Candidate& b) {
    const double da = a.cost;
    const double db = b.cost;
    if (std::abs(da - db) > 1e-12) return da < db;
    return pathLexLess(a.path, b.path);
  };

  std::vector<Candidate> candidates;

  // Yen's algorithm:
  //  - out[0] is the shortest.
  //  - candidates holds the next-best deviations.
  for (std::size_t kth = 1; kth < k; ++kth) {
    const auto& prev = out[kth - 1].path;
    if (prev.size() < 2) break;

    for (std::size_t i = 0; i + 1 < prev.size(); ++i) {
      const SystemId spurNode = prev[i];

      // Root path includes spur node.
      std::vector<SystemId> root(prev.begin(), prev.begin() + (long long)i + 1);

      // Ban nodes in the root path *except* the spur node to enforce loopless paths.
      std::unordered_set<SystemId> bannedNodes;
      bannedNodes.reserve(i);
      for (std::size_t r = 0; r < i; ++r) bannedNodes.insert(root[r]);

      // Ban edges that would recreate any previously found shortest path with this root.
      std::unordered_set<Edge, EdgeHash> bannedEdges;
      for (const auto& found : out) {
        const auto& p = found.path;
        if (p.size() > i + 1 && pathHasPrefix(p, root)) {
          bannedEdges.insert(Edge{p[i], p[i + 1]});
        }
      }

      const auto spur = aStarCostSolveHazards(nodes, idx, grid,
                                             spurNode, goalId,
                                             maxJumpLy,
                                             costPerJump, costPerLy, hazardWeightPerLy,
                                             universeSeed, timeDays,
                                             &bannedNodes,
                                             &bannedEdges,
                                             maxExpansionsPerSolve);

      if (spur.path.empty()) continue;

      // Combine root + spur (skip spurNode duplicate).
      std::vector<SystemId> total = root;
      if (spur.path.size() > 1) {
        total.insert(total.end(), spur.path.begin() + 1, spur.path.end());
      }

      if (total.size() < 2) continue;
      if (containsPath(out, total)) continue;
      if (containsPathCand(candidates, total)) continue;

      const double dist = computeDistanceLy(idx, nodes, total);
      const double cst  = computeCostHazards(idx, nodes, total,
                                             costPerJump, costPerLy, hazardWeightPerLy,
                                             universeSeed, timeDays);

      const int hops = (int)total.size() - 1;
      candidates.push_back(Candidate{std::move(total), hops, dist, cst});
    }

    if (candidates.empty()) break;

    // Pick the best candidate.
    std::size_t bestIdx = 0;
    for (std::size_t i = 1; i < candidates.size(); ++i) {
      if (betterCand(candidates[i], candidates[bestIdx])) bestIdx = i;
    }

    Candidate best = std::move(candidates[bestIdx]);
    candidates.erase(candidates.begin() + (long long)bestIdx);

    out.push_back(KRoute{std::move(best.path), best.hops, best.distanceLy, best.cost});
  }

  return out;
}

std::vector<KRoute> plotKRoutesAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve) {
  return plotKRoutesAStarCost(nodes, startId, goalId, maxJumpLy, 1.0, 0.0, k, maxExpansionsPerSolve);
}

std::vector<KRoute> selectDiversifiedKRoutesMMR(std::span<const KRoute> candidates,
                                                std::size_t k,
                                                double diversityWeight01,
                                                bool ignoreEndpoints) {
  std::vector<KRoute> out;
  if (k == 0) return out;
  if (candidates.empty()) return out;

  double alpha = diversityWeight01;
  if (!std::isfinite(alpha)) alpha = 0.0;
  alpha = std::clamp(alpha, 0.0, 1.0);

  // Deduplicate by path (keep the lowest-cost entry).
  std::vector<KRoute> uniq;
  uniq.reserve(candidates.size());

  for (const KRoute& c : candidates) {
    if (c.path.empty()) continue;

    bool merged = false;
    for (auto& u : uniq) {
      if (u.path == c.path) {
        // Prefer lower cost; if costs tie, prefer lexicographically smaller path.
        const double du = u.cost;
        const double dc = c.cost;
        if ((dc < du - 1e-12) || (std::abs(dc - du) <= 1e-12 && pathLexLess(c.path, u.path))) {
          u = c;
        }
        merged = true;
        break;
      }
    }
    if (!merged) {
      uniq.push_back(c);
    }
  }

  if (uniq.empty()) return out;

  // Stable total order: (cost, then path lex).
  std::sort(uniq.begin(), uniq.end(), [&](const KRoute& a, const KRoute& b) {
    const double da = a.cost;
    const double db = b.cost;
    if (std::abs(da - db) > 1e-12) return da < db;
    return pathLexLess(a.path, b.path);
  });

  if (k > uniq.size()) k = uniq.size();
  out.reserve(k);

  // Compute cost normalization range.
  double minCost = uniq.front().cost;
  double maxCost = uniq.front().cost;
  for (const auto& r : uniq) {
    if (!std::isfinite(r.cost)) continue;
    minCost = std::min(minCost, r.cost);
    maxCost = std::max(maxCost, r.cost);
  }
  if (!std::isfinite(minCost)) minCost = 0.0;
  if (!std::isfinite(maxCost)) maxCost = minCost;

  const auto quality01 = [&](double cost) {
    if (!std::isfinite(cost)) return 0.0;
    if (maxCost <= minCost + 1e-12) return 1.0;
    double t = (cost - minCost) / (maxCost - minCost);
    t = std::clamp(t, 0.0, 1.0);
    return 1.0 - t;
  };

  // Greedy MMR-like selection.
  std::vector<bool> used(uniq.size(), false);

  // Always take the best-cost route first.
  out.push_back(uniq[0]);
  used[0] = true;

  for (std::size_t step = 1; step < k; ++step) {
    std::size_t bestIdx = (std::size_t)-1;
    double bestScore = -1e300;

    for (std::size_t i = 0; i < uniq.size(); ++i) {
      if (used[i]) continue;

      const KRoute& cand = uniq[i];
      const double q = quality01(cand.cost);

      double maxSim = 0.0;
      for (const auto& sel : out) {
        maxSim = std::max(maxSim, routeNodeJaccard01(cand.path, sel.path, ignoreEndpoints));
      }
      maxSim = std::clamp(maxSim, 0.0, 1.0);
      const double diversity01 = 1.0 - maxSim;
      const double score = (1.0 - alpha) * q + alpha * diversity01;

      bool take = false;
      if (bestIdx == (std::size_t)-1) {
        take = true;
      } else if (score > bestScore + 1e-12) {
        take = true;
      } else if (std::abs(score - bestScore) <= 1e-12) {
        // Deterministic tie-break.
        const KRoute& best = uniq[bestIdx];
        if (cand.cost < best.cost - 1e-12) {
          take = true;
        } else if (std::abs(cand.cost - best.cost) <= 1e-12 && pathLexLess(cand.path, best.path)) {
          take = true;
        }
      }

      if (take) {
        bestIdx = i;
        bestScore = score;
      }
    }

    if (bestIdx == (std::size_t)-1) break;
    used[bestIdx] = true;
    out.push_back(uniq[bestIdx]);
  }

  return out;
}

std::vector<KRoute> plotKRoutesAStarCostDiversified(const std::vector<SystemStub>& nodes,
                                                   SystemId startId,
                                                   SystemId goalId,
                                                   double maxJumpLy,
                                                   double costPerJump,
                                                   double costPerLy,
                                                   std::size_t k,
                                                   double diversityWeight01,
                                                   std::size_t kCandidates,
                                                   bool ignoreEndpoints,
                                                   std::size_t maxExpansionsPerSolve) {
  if (k == 0) return {};

  if (kCandidates == 0) {
    kCandidates = std::max<std::size_t>(k, k * 6);
  }
  if (kCandidates < k) kCandidates = k;

  const auto candidates = plotKRoutesAStarCost(nodes,
                                               startId,
                                               goalId,
                                               maxJumpLy,
                                               costPerJump,
                                               costPerLy,
                                               kCandidates,
                                               maxExpansionsPerSolve);

  return selectDiversifiedKRoutesMMR(candidates, k, diversityWeight01, ignoreEndpoints);
}

double routeDistanceLy(const std::vector<SystemStub>& nodes,
                       const std::vector<SystemId>& route) {
  if (route.size() < 2) return 0.0;

  const auto idx = buildIndex(nodes);

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    sum += systemDistanceLy(nodes[itA->second], nodes[itB->second]);
  }
  return sum;
}

double routeCost(const std::vector<SystemStub>& nodes,
                 const std::vector<SystemId>& route,
                 double costPerJump,
                 double costPerLy) {
  if (route.size() < 2) return 0.0;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  const auto idx = buildIndex(nodes);

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    const double d = systemDistanceLy(nodes[itA->second], nodes[itB->second]);
    sum += costPerJump + costPerLy * d;
  }
  return sum;
}

double routeCostRisk(const std::vector<SystemStub>& nodes,
                     const std::vector<SystemId>& route,
                     double costPerJump,
                     double costPerLy,
                     double riskWeightPerLy,
                     std::span<const double> risk01PerNode) {
  if (route.size() < 2) return 0.0;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;
  if (riskWeightPerLy < 0.0) riskWeightPerLy = 0.0;

  const bool hasRisk = (risk01PerNode.size() == nodes.size());
  auto riskAt = [&](std::size_t i) -> double {
    if (!hasRisk) return 0.0;
    const double r = risk01PerNode[i];
    if (r < 0.0) return 0.0;
    if (r > 1.0) return 1.0;
    return r;
  };

  const auto idx = buildIndex(nodes);

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;

    const std::size_t ia = itA->second;
    const std::size_t ib = itB->second;
    const double d = systemDistanceLy(nodes[ia], nodes[ib]);
    const double avgRisk01 = 0.5 * (riskAt(ia) + riskAt(ib));
    sum += costPerJump + costPerLy * d + riskWeightPerLy * avgRisk01 * d;
  }
  return sum;
}


double routeCostHazards(const std::vector<SystemStub>& nodes,
                        const std::vector<SystemId>& route,
                        double costPerJump,
                        double costPerLy,
                        double hazardWeightPerLy,
                        core::u64 universeSeed,
                        double timeDays) {
  if (route.size() < 2) return 0.0;

  const auto idx = buildIndex(nodes);
  return computeCostHazards(idx, nodes, route,
                            costPerJump,
                            costPerLy,
                            hazardWeightPerLy,
                            universeSeed,
                            timeDays);
}

double routeNodeJaccard01(const std::vector<SystemId>& a,
                          const std::vector<SystemId>& b,
                          bool ignoreEndpoints) {
  auto buildSet = [&](const std::vector<SystemId>& r) {
    std::unordered_set<SystemId> s;
    if (r.empty()) return s;

    std::size_t begin = 0;
    std::size_t end = r.size();
    if (ignoreEndpoints) {
      if (r.size() <= 2) return s;
      begin = 1;
      end = r.size() - 1;
    }

    for (std::size_t i = begin; i < end; ++i) {
      s.insert(r[i]);
    }
    return s;
  };

  const auto sa = buildSet(a);
  const auto sb = buildSet(b);

  if (sa.empty() && sb.empty()) return 1.0;

  std::size_t inter = 0;
  if (!sa.empty() && !sb.empty()) {
    // Iterate the smaller set for a modest speedup.
    const auto* small = &sa;
    const auto* large = &sb;
    if (sb.size() < sa.size()) {
      small = &sb;
      large = &sa;
    }
    for (const auto id : *small) {
      if (large->find(id) != large->end()) {
        ++inter;
      }
    }
  }

  const std::size_t uni = sa.size() + sb.size() - inter;
  if (uni == 0) return 1.0;

  const double j = (double)inter / (double)uni;
  if (j < 0.0) return 0.0;
  if (j > 1.0) return 1.0;
  return j;
}


bool validateRoute(const std::vector<SystemStub>& nodes,
                   const std::vector<SystemId>& route,
                   double maxJumpLy,
                   std::string* outError) {
  if (route.empty()) {
    if (outError) *outError = "route is empty";
    return false;
  }
  if (maxJumpLy <= 0.0) {
    if (outError) *outError = "maxJumpLy must be > 0";
    return false;
  }

  const auto idx = buildIndex(nodes);

  for (const auto id : route) {
    if (idx.find(id) == idx.end()) {
      if (outError) *outError = "route references unknown system id: " + std::to_string((unsigned long long)id);
      return false;
    }
  }

  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto a = idx.find(route[i]);
    const auto b = idx.find(route[i + 1]);
    if (a == idx.end() || b == idx.end()) continue;
    const double d = systemDistanceLy(nodes[a->second], nodes[b->second]);
    if (d > maxJumpLy + 1e-9) {
      if (outError) {
        *outError = "jump " + std::to_string(i) + " exceeds maxJumpLy (" +
                    std::to_string(d) + " > " + std::to_string(maxJumpLy) + ")";
      }
      return false;
    }
  }

  return true;
}

} // namespace stellar::sim
