#include "stellar/sim/NavRouteResilience.h"

#include "test_harness.h"

#include <cmath>
#include <string>
#include <vector>

using stellar::math::Vec3d;
using stellar::sim::RouteEdgeDetour;
using stellar::sim::RouteResilienceReport;
using stellar::sim::SystemId;
using stellar::sim::SystemStub;

namespace {

SystemStub stub(SystemId id, double x) {
  SystemStub s;
  s.id = id;
  s.posLy = Vec3d{x, 0.0, 0.0};
  return s;
}

bool routeUsesUndirectedEdge(const std::vector<SystemId>& route, SystemId a, SystemId b) {
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const SystemId u = route[i];
    const SystemId v = route[i + 1];
    if ((u == a && v == b) || (u == b && v == a)) return true;
  }
  return false;
}

} // namespace

int test_nav_route_resilience() {
  int failures = 0;

  // Graph layout (1D for easy jump-range control):
  //
  //   1 -- 2 -- 3 -- 4 -- 5   (range-limited, base path)
  //         \  /
  //          6      (local detour only for edge 2-3)
  //
  std::vector<SystemStub> nodes;
  nodes.push_back(stub(1, 0.0));
  nodes.push_back(stub(2, 1.0));
  nodes.push_back(stub(3, 2.0));
  nodes.push_back(stub(4, 3.0));
  nodes.push_back(stub(5, 4.0));
  nodes.push_back(stub(6, 1.5)); // within range of both 2 and 3

  const double maxJumpLy = 1.1;
  const double costPerJump = 1.0;
  const double costPerLy = 0.0;

  const RouteResilienceReport rep =
      stellar::sim::analyzeRouteEdgeDetoursCost(nodes, 1, 5, maxJumpLy, costPerJump, costPerLy);

  CHECK(rep.ok);
  CHECK_EQ(rep.basePath.size(), (std::size_t)5);

  CHECK_EQ(rep.basePath[0], (SystemId)1);
  CHECK_EQ(rep.basePath[1], (SystemId)2);
  CHECK_EQ(rep.basePath[2], (SystemId)3);
  CHECK_EQ(rep.basePath[3], (SystemId)4);
  CHECK_EQ(rep.basePath[4], (SystemId)5);

  CHECK_EQ(rep.edgeDetours.size(), (std::size_t)4);

  // Edge 1-2 is critical (no alternate neighbor within range).
  CHECK(!rep.edgeDetours[0].ok);

  // Edge 2-3 has a detour via 6.
  CHECK(rep.edgeDetours[1].ok);
  CHECK(routeUsesUndirectedEdge(rep.edgeDetours[1].path, 2, 3) == false);

  // Expect the exact deterministic detour sequence.
  CHECK_EQ(rep.edgeDetours[1].path.size(), (std::size_t)6);
  CHECK_EQ(rep.edgeDetours[1].path[0], (SystemId)1);
  CHECK_EQ(rep.edgeDetours[1].path[1], (SystemId)2);
  CHECK_EQ(rep.edgeDetours[1].path[2], (SystemId)6);
  CHECK_EQ(rep.edgeDetours[1].path[3], (SystemId)3);
  CHECK_EQ(rep.edgeDetours[1].path[4], (SystemId)4);
  CHECK_EQ(rep.edgeDetours[1].path[5], (SystemId)5);

  // Remaining edges are critical in this toy graph.
  CHECK(!rep.edgeDetours[2].ok);
  CHECK(!rep.edgeDetours[3].ok);

  // Base cost is 4 hops (costPerJump=1). Detour adds one hop (via 6).
  CHECK(std::abs(rep.baseCost - 4.0) < 1e-9);
  CHECK(std::abs(rep.edgeDetours[1].cost - 5.0) < 1e-9);
  CHECK(std::abs(rep.edgeDetours[1].costFactor - 1.25) < 1e-9);

  // Contiguity sanity checks (detour route should still be valid).
  std::string err;
  CHECK(stellar::sim::validateRoute(nodes, rep.basePath, maxJumpLy, &err));
  CHECK(stellar::sim::validateRoute(nodes, rep.edgeDetours[1].path, maxJumpLy, &err));

  return failures;
}
