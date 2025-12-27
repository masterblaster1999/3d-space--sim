#include "stellar/sim/Universe.h"

#include <algorithm>
#include <iostream>

int test_proc() {
  int fails = 0;

  const stellar::core::u64 seed = 12345;
  stellar::sim::Universe a(seed);
  stellar::sim::Universe b(seed);

  const stellar::math::Vec3d pos{0,0,0};
  const auto la = a.queryNearby(pos, 60.0, 16);
  const auto lb = b.queryNearby(pos, 60.0, 16);

  if (la.empty() || lb.empty()) {
    std::cerr << "[test_proc] expected at least one system stub\n";
    return 1;
  }

  // Determinism: nearby query should match for the same seed.
  if (la.size() != lb.size()) {
    std::cerr << "[test_proc] size mismatch " << la.size() << " vs " << lb.size() << "\n";
    ++fails;
  }

  const std::size_t n = std::min<std::size_t>({la.size(), lb.size(), 8});
  for (std::size_t i = 0; i < n; ++i) {
    if (la[i].id != lb[i].id || la[i].name != lb[i].name) {
      std::cerr << "[test_proc] stub mismatch at i=" << i << "\n";
      ++fails;
      break;
    }
  }

  // Determinism: full system generation should match for the same seed + id.
  {
    const auto& sa = a.getSystem(la.front().id, &la.front());
    const auto& sb = b.getSystem(lb.front().id, &lb.front());

    if (sa.stub.id != sb.stub.id || sa.stub.name != sb.stub.name) {
      std::cerr << "[test_proc] system stub mismatch\n";
      ++fails;
    }
    if (sa.planets.size() != sb.planets.size()) {
      std::cerr << "[test_proc] planet count mismatch\n";
      ++fails;
    }
    if (sa.stations.size() != sb.stations.size()) {
      std::cerr << "[test_proc] station count mismatch\n";
      ++fails;
    }
    if (!sa.planets.empty() && !sb.planets.empty() && sa.planets.front().name != sb.planets.front().name) {
      std::cerr << "[test_proc] first planet name mismatch\n";
      ++fails;
    }
    if (!sa.stations.empty() && !sb.stations.empty() && sa.stations.front().name != sb.stations.front().name) {
      std::cerr << "[test_proc] first station name mismatch\n";
      ++fails;
    }
  }

  // Different seed should very likely produce a different first nearby system.
  {
    stellar::sim::Universe c(seed + 1);
    const auto lc = c.queryNearby(pos, 60.0, 16);
    if (!lc.empty()) {
      if (lc.front().id == la.front().id && lc.front().name == la.front().name) {
        std::cerr << "[test_proc] different seed produced identical first stub (unexpected)\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_proc] pass\n";
  return fails;
}
