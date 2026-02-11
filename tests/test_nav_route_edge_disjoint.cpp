#include "stellar/sim/NavRouteConnectivity.h"

#include "test_harness.h"

#include <algorithm>
#include <set>
#include <vector>

using stellar::sim::EdgeDisjointRoutesResult;
using stellar::sim::RouteEdge;
using stellar::sim::SystemId;
using stellar::sim::SystemStub;

namespace {

static SystemStub stub(SystemId id, double x, double y, double z = 0.0) {
  SystemStub s;
  s.id = id;
  s.posLy = stellar::math::Vec3d{x, y, z};
  return s;
}

static std::pair<SystemId, SystemId> canon(SystemId a, SystemId b) {
  if (a < b) return {a, b};
  return {b, a};
}

} // namespace

int test_nav_route_edge_disjoint() {
  int failures = 0;

  // Diamond-like graph with 3 independent corridors:
  //   1-2-5
  //   1-3-5
  //   1-4-5
  // plus extra cross edges among {2,3,4} that should not change the max
  // edge-disjoint count (still limited by deg(1)=3).
  std::vector<SystemStub> nodes;
  nodes.push_back(stub(1, 0.0, 0.0));
  nodes.push_back(stub(2, 1.0, 0.0));
  nodes.push_back(stub(3, 1.0, 1.0));
  nodes.push_back(stub(4, 1.0, -1.0));
  nodes.push_back(stub(5, 2.0, 0.0));

  const double maxJumpLy = 1.5;

  {
    const EdgeDisjointRoutesResult r = stellar::sim::computeEdgeDisjointRoutes(nodes, 1, 5, maxJumpLy, 8);
    CHECK(r.reachable);
    CHECK_EQ(r.maxRoutes, 3);
    CHECK_EQ(r.routes.size(), (std::size_t)3);

    // Deterministic lex ordering.
    CHECK(r.routes[0] == (std::vector<SystemId>{1, 2, 5}));
    CHECK(r.routes[1] == (std::vector<SystemId>{1, 3, 5}));
    CHECK(r.routes[2] == (std::vector<SystemId>{1, 4, 5}));

    // Verify edge-disjointness in the undirected sense.
    std::set<std::pair<SystemId, SystemId>> used;
    for (const auto& path : r.routes) {
      CHECK(path.size() >= 2);
      for (std::size_t i = 0; i + 1 < path.size(); ++i) {
        const auto e = canon(path[i], path[i + 1]);
        CHECK(used.insert(e).second);
      }
    }

    // The minimum cut should be the 3 saturated edges out of the start.
    CHECK_EQ(r.minCutEdges.size(), (std::size_t)3);
    CHECK(r.minCutEdges[0].from == 1 && r.minCutEdges[0].to == 2);
    CHECK(r.minCutEdges[1].from == 1 && r.minCutEdges[1].to == 3);
    CHECK(r.minCutEdges[2].from == 1 && r.minCutEdges[2].to == 4);
  }

  {
    // Ban one corridor (edge 1-2) and confirm the max disjoint routes drop.
    const std::vector<RouteEdge> bannedEdges = {RouteEdge{1, 2}};
    const EdgeDisjointRoutesResult r = stellar::sim::computeEdgeDisjointRoutes(nodes, 1, 5, maxJumpLy, 8, {}, bannedEdges);
    CHECK(r.reachable);
    CHECK_EQ(r.maxRoutes, 2);
    CHECK_EQ(r.routes.size(), (std::size_t)2);
    CHECK(r.routes[0] == (std::vector<SystemId>{1, 3, 5}));
    CHECK(r.routes[1] == (std::vector<SystemId>{1, 4, 5}));
    CHECK_EQ(r.minCutEdges.size(), (std::size_t)2);
    CHECK(r.minCutEdges[0].from == 1 && r.minCutEdges[0].to == 3);
    CHECK(r.minCutEdges[1].from == 1 && r.minCutEdges[1].to == 4);
  }

  {
    // If the start is banned, the result should be unreachable.
    const std::vector<SystemId> bannedNodes = {1};
    const EdgeDisjointRoutesResult r = stellar::sim::computeEdgeDisjointRoutes(nodes, 1, 5, maxJumpLy, 8, bannedNodes, {});
    CHECK(!r.reachable);
    CHECK_EQ(r.maxRoutes, 0);
    CHECK(r.routes.empty());
  }

  return failures;
}
