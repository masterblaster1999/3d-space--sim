#include "stellar/sim/Universe.h"

#include <iostream>

int test_proc() {
  int fails = 0;

  stellar::sim::UniverseConfig cfg;
  cfg.seed = 12345;
  cfg.systemCount = 10;

  stellar::sim::Universe a(cfg);
  stellar::sim::Universe b(cfg);

  a.generate();
  b.generate();

  if (a.size() != b.size()) {
    std::cerr << "[test_proc] size mismatch\n";
    ++fails;
  }

  // Determinism check: first system star name + planet count should match.
  if (a.size() > 0 && b.size() > 0) {
    const auto& sa = a.system(0);
    const auto& sb = b.system(0);

    if (sa.primary.name != sb.primary.name) {
      std::cerr << "[test_proc] star name mismatch: " << sa.primary.name << " vs " << sb.primary.name << "\n";
      ++fails;
    }
    if (sa.planets.size() != sb.planets.size()) {
      std::cerr << "[test_proc] planet count mismatch\n";
      ++fails;
    }
    if (sa.id != sb.id) {
      std::cerr << "[test_proc] id mismatch\n";
      ++fails;
    }
  }

  // Different seed should produce different output (very likely).
  stellar::sim::UniverseConfig cfg2 = cfg;
  cfg2.seed = 54321;
  stellar::sim::Universe c(cfg2);
  c.generate();

  if (a.size() > 0 && c.size() > 0) {
    if (a.system(0).primary.name == c.system(0).primary.name && a.system(0).id == c.system(0).id) {
      std::cerr << "[test_proc] different seed produced identical first system (unexpected)\n";
      ++fails;
    }
  }

  return fails;
}
