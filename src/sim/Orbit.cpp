#include "stellar/sim/Orbit.h"

#include "stellar/math/Math.h"

#include <cmath>

namespace stellar::sim {
namespace {
  double wrapTwoPi(double x) {
    x = std::fmod(x, stellar::math::twoPi);
    if (x < 0.0) x += stellar::math::twoPi;
    return x;
  }

  // Solve Kepler's equation: E - e sin(E) = M
  // Returns E (eccentric anomaly) in radians.
  double solveKepler(double M, double e) {
    // Good initial guess:
    double E = (e < 0.8) ? M : stellar::math::pi;

    // Newton-Raphson iterations
    for (int i = 0; i < 12; ++i) {
      const double f = E - e * std::sin(E) - M;
      const double fp = 1.0 - e * std::cos(E);
      const double dE = f / fp;
      E -= dE;
      if (std::abs(dE) < 1e-12) break;
    }
    return E;
  }
} // namespace

stellar::math::Vec3d Orbit::positionAU(double daysSinceEpoch) const {
  const double a = semiMajorAxisAU;
  const double e = stellar::math::clamp(eccentricity, 0.0, 0.999999);

  const double P = (orbitalPeriodDays <= 0.0) ? 1.0 : orbitalPeriodDays;
  const double n = stellar::math::twoPi / P; // mean motion rad/day

  const double M0 = stellar::math::degToRad(meanAnomalyAtEpochDeg);
  const double M = wrapTwoPi(M0 + n * daysSinceEpoch);

  const double E = solveKepler(M, e);

  // True anomaly ν
  const double sinE2 = std::sin(E * 0.5);
  const double cosE2 = std::cos(E * 0.5);
  const double sqrt1pe = std::sqrt(1.0 + e);
  const double sqrt1me = std::sqrt(1.0 - e);
  const double nu = 2.0 * std::atan2(sqrt1pe * sinE2, sqrt1me * cosE2);

  // Radius
  const double r = a * (1.0 - e * std::cos(E));

  // Angles
  const double i = stellar::math::degToRad(inclinationDeg);
  const double Omega = stellar::math::degToRad(longitudeAscendingNodeDeg);
  const double omega = stellar::math::degToRad(argumentPeriapsisDeg);

  const double cosO = std::cos(Omega);
  const double sinO = std::sin(Omega);
  const double cosi = std::cos(i);
  const double sini = std::sin(i);

  const double arg = omega + nu;
  const double cosA = std::cos(arg);
  const double sinA = std::sin(arg);

  // Inertial coordinates (standard orbital element transform)
  const double x = r * (cosO * cosA - sinO * sinA * cosi);
  const double y = r * (sinO * cosA + cosO * sinA * cosi);
  const double z = r * (sinA * sini);

  return {x, y, z};
}

} // namespace stellar::sim
