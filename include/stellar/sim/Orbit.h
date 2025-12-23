#pragma once

#include "stellar/math/Vec3.h"

namespace stellar::sim {

// Simple Keplerian orbit (two-body) for procedural solar system simulation.
// Units:
// - semiMajorAxisAU: AU
// - orbitalPeriodDays: days
// - angles: degrees
struct Orbit {
  double semiMajorAxisAU = 1.0;
  double eccentricity = 0.0;

  double inclinationDeg = 0.0;              // i
  double longitudeAscendingNodeDeg = 0.0;   // Ω
  double argumentPeriapsisDeg = 0.0;        // ω
  double meanAnomalyAtEpochDeg = 0.0;       // M0 at t=0

  double orbitalPeriodDays = 365.25;

  // Returns position in AU relative to the primary body.
  stellar::math::Vec3d positionAU(double daysSinceEpoch) const;
};

} // namespace stellar::sim
