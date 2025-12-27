#include "stellar/sim/NavRoute.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>

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

} // namespace

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
    setStats(outStats, stats);
    return {startId};
  }

  // Spatial hash grid for neighbor queries.
  // Cell size is maxJumpLy; any reachable neighbor must be in the same or adjacent cell.
  const double cellSize = maxJumpLy;
  std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash> grid;
  grid.reserve(nodes.size());

  for (std::size_t i = 0; i < nodes.size(); ++i) {
    grid[cellFor(nodes[i].posLy, cellSize)].push_back(i);
  }

  std::vector<int> cameFrom(N, -1);
  std::vector<int> gScore(N, std::numeric_limits<int>::max());
  std::vector<char> closed(N, 0);

  auto heuristic = [&](std::size_t i) -> int {
    if (i == goal) return 0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    return (int)std::ceil(d / maxJumpLy);
  };

  struct QN {
    int f;
    int g;
    std::size_t i;
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const { return a.f > b.f; }
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

      setStats(outStats, stats);
      return path;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, cellSize);

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
