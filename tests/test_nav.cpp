#include "stellar/sim/NavRoute.h"

#include <iostream>
#include <string>
#include <vector>

int test_nav() {
  int fails = 0;

  using namespace stellar;

  auto makeStub = [](sim::SystemId id, double x) {
    sim::SystemStub s{};
    s.id = id;
    s.seed = id * 1337ULL;
    s.name = "S" + std::to_string((unsigned long long)id);
    s.posLy = math::Vec3d{x, 0.0, 0.0};
    s.primaryClass = sim::StarClass::G;
    s.planetCount = 1;
    s.stationCount = 1;
    s.factionId = 0;
    return s;
  };

  // A simple linear chain of systems 8 ly apart.
  std::vector<sim::SystemStub> nodes;
  nodes.push_back(makeStub(1, 0.0));
  nodes.push_back(makeStub(2, 8.0));
  nodes.push_back(makeStub(3, 16.0));
  nodes.push_back(makeStub(4, 24.0));
  nodes.push_back(makeStub(5, 32.0));
  nodes.push_back(makeStub(6, 40.0));
  nodes.push_back(makeStub(99, 1000.0)); // distractor far away

  sim::RoutePlanStats stats{};
  const auto route = sim::plotRouteAStarHops(nodes, 1, 6, 10.0, &stats);

  if (route.empty()) {
    std::cerr << "[test_nav] expected a route, got empty\n";
    ++fails;
  } else {
    if (route.front() != 1 || route.back() != 6) {
      std::cerr << "[test_nav] route endpoints mismatch\n";
      ++fails;
    }
    if (route.size() != 6) {
      std::cerr << "[test_nav] expected 6 nodes in route, got " << route.size() << "\n";
      ++fails;
    }
    std::string err;
    if (!sim::validateRoute(nodes, route, 10.0, &err)) {
      std::cerr << "[test_nav] route validation failed: " << err << "\n";
      ++fails;
    }
    if (!stats.reached || stats.hops != 5) {
      std::cerr << "[test_nav] stats mismatch reached=" << stats.reached << " hops=" << stats.hops << "\n";
      ++fails;
    }

    // Determinism check (same input => same output).
    sim::RoutePlanStats stats2{};
    const auto route2 = sim::plotRouteAStarHops(nodes, 1, 6, 10.0, &stats2);
    if (route2 != route) {
      std::cerr << "[test_nav] non-deterministic result\n";
      ++fails;
    }
  }

  // Unreachable: gap larger than jump range.
  std::vector<sim::SystemStub> gap;
  gap.push_back(makeStub(1, 0.0));
  gap.push_back(makeStub(2, 8.0));
  gap.push_back(makeStub(4, 24.0)); // gap (8 -> 24) = 16 > 10
  gap.push_back(makeStub(6, 40.0));

  const auto noRoute = sim::plotRouteAStarHops(gap, 1, 6, 10.0, nullptr);
  if (!noRoute.empty()) {
    std::cerr << "[test_nav] expected no route in gapped graph, got size=" << noRoute.size() << "\n";
    ++fails;
  }

  // start==goal: should trivially succeed (if present).
  const auto self = sim::plotRouteAStarHops(nodes, 3, 3, 10.0, nullptr);
  if (self.size() != 1 || self.front() != 3) {
    std::cerr << "[test_nav] expected self route [3], got size=" << self.size() << "\n";
    ++fails;
  }

  if (fails == 0) std::cout << "[test_nav] pass\n";
  return fails;
}
