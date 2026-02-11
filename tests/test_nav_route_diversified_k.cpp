#include "stellar/sim/NavRoute.h"

#include "test_harness.h"

#include <algorithm>
#include <vector>

using stellar::sim::KRoute;
using stellar::sim::SystemId;

namespace {

bool samePath(const std::vector<SystemId>& a, const std::vector<SystemId>& b) {
  return a == b;
}

} // namespace

int test_nav_route_diversified_k() {
  int failures = 0;

  // Synthetic candidate routes (paths + costs). This isolates the selection
  // logic from the underlying A* / geometry graph.
  std::vector<KRoute> cands;
  cands.push_back(KRoute{{1, 2, 3, 10}, 3, 0.0, 10.0});
  cands.push_back(KRoute{{1, 2, 4, 10}, 3, 0.0, 10.1});
  cands.push_back(KRoute{{1, 5, 6, 10}, 3, 0.0, 11.0});
  cands.push_back(KRoute{{1, 2, 3, 7, 10}, 4, 0.0, 10.5});

  // With diversityWeight=0, selection should be cost-ordered.
  {
    const auto picked = stellar::sim::selectDiversifiedKRoutesMMR(cands, 2, 0.0, true);
    CHECK_EQ(picked.size(), (std::size_t)2);
    CHECK(samePath(picked[0].path, {1, 2, 3, 10}));
    CHECK(samePath(picked[1].path, {1, 2, 4, 10}));
  }

  // With diversityWeight=1, selection should prefer the most different route.
  {
    const auto picked = stellar::sim::selectDiversifiedKRoutesMMR(cands, 2, 1.0, true);
    CHECK_EQ(picked.size(), (std::size_t)2);
    CHECK(samePath(picked[0].path, {1, 2, 3, 10}));
    CHECK(samePath(picked[1].path, {1, 5, 6, 10}));
  }

  // The greedy choice should remain deterministic even if the input order changes.
  {
    auto c2 = cands;
    std::reverse(c2.begin(), c2.end());

    const auto picked = stellar::sim::selectDiversifiedKRoutesMMR(c2, 3, 1.0, true);
    CHECK_EQ(picked.size(), (std::size_t)3);
    CHECK(samePath(picked[0].path, {1, 2, 3, 10}));
    CHECK(samePath(picked[1].path, {1, 5, 6, 10}));
    CHECK(samePath(picked[2].path, {1, 2, 4, 10}));
  }

  return failures;
}
