#include "stellar/sim/Orbit.h"

#include <cmath>
#include <iostream>

namespace {
bool approx(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}
} // namespace

int test_orbit() {
  using stellar::sim::Orbit;

  int fails = 0;

  Orbit o;
  o.semiMajorAxisAU = 2.0;
  o.eccentricity = 0.0;
  o.inclinationDeg = 0.0;
  o.longitudeAscendingNodeDeg = 0.0;
  o.argumentPeriapsisDeg = 0.0;
  o.meanAnomalyAtEpochDeg = 0.0;
  o.orbitalPeriodDays = 100.0;

  auto p0 = o.positionAU(0.0);
  auto p1 = o.positionAU(25.0);  // quarter period
  auto p2 = o.positionAU(50.0);  // half period

  if (!approx(p0.x, 2.0) || !approx(p0.y, 0.0) || !approx(p0.z, 0.0)) {
    std::cerr << "[test_orbit] expected p0=(2,0,0) got " << p0 << "\n";
    ++fails;
  }

  if (!approx(p1.x, 0.0, 1e-6) || !approx(p1.y, 2.0, 1e-6) || !approx(p1.z, 0.0, 1e-6)) {
    std::cerr << "[test_orbit] expected p1≈(0,2,0) got " << p1 << "\n";
    ++fails;
  }

  if (!approx(p2.x, -2.0, 1e-6) || !approx(p2.y, 0.0, 1e-6) || !approx(p2.z, 0.0, 1e-6)) {
    std::cerr << "[test_orbit] expected p2≈(-2,0,0) got " << p2 << "\n";
    ++fails;
  }

  return fails;
}
