#include "stellar/sim/BearingTrack.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

static bool finitePos(double x) {
  return std::isfinite(x) && x > 0.0;
}

static math::Vec3d col(const math::Mat3d& m, int i) {
  // Mat3d is column-major: m[col * 3 + row]
  const int base = i * 3;
  return {m.m[base + 0], m.m[base + 1], m.m[base + 2]};
}

static bool axesOrthonormal(const math::Mat3d& axes, double eps = 1e-6) {
  const math::Vec3d a0 = col(axes, 0);
  const math::Vec3d a1 = col(axes, 1);
  const math::Vec3d a2 = col(axes, 2);

  if (!approx(a0.length(), 1.0, eps)) return false;
  if (!approx(a1.length(), 1.0, eps)) return false;
  if (!approx(a2.length(), 1.0, eps)) return false;

  if (std::abs(math::dot(a0, a1)) > eps) return false;
  if (std::abs(math::dot(a0, a2)) > eps) return false;
  if (std::abs(math::dot(a1, a2)) > eps) return false;

  return true;
}

int test_geometry_bearing_uncertainty_ellipsoid() {
  int fails = 0;

  sim::BearingTrackParams tp{};
  tp.observationHalfLifeSec = 0.0;
  tp.velHalfLifeSec = 0.0;
  tp.solveRegularization = 0.0;
  tp.determinantEps = 1.0e-14;
  tp.minEffectiveWeight = 1.8;

  sim::BearingTrackUncertaintyEllipsoidParams ep{};
  ep.sigmaScaleKm = 1.0;
  ep.minInfoEigenvalue = 1.0e-12;
  ep.maxAxisKm = 1.0e12;

  const math::Vec3d target{0.0, 0.0, 100000.0};

  auto makeTrack = [&](double w) {
    sim::BearingTrack3d tr{};

    const math::Vec3d o1{0.0, 0.0, 0.0};
    const math::Vec3d o2{20000.0, 0.0, 0.0};

    const math::Vec3d d1 = (target - o1).normalized();
    const math::Vec3d d2 = (target - o2).normalized();

    sim::updateBearingTrack(tr, 0.0, true, o1, d1, w, tp);
    sim::updateBearingTrack(tr, 0.0, true, o2, d2, w, tp);

    return tr;
  };

  // --- Solvable track: ellipsoid should be valid and orthonormal ---
  {
    const auto tr = makeTrack(1.0);
    if (!tr.initialized) {
      std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] expected initialized track with two bearings.\n";
      ++fails;
    }

    const auto ell = sim::bearingTrackUncertaintyEllipsoid(tr, ep);
    if (!ell.valid) {
      std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] expected valid ellipsoid.\n";
      ++fails;
    } else {
      if (!axesOrthonormal(ell.axes, 1e-5)) {
        std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] axes not orthonormal.\n";
        ++fails;
      }

      if (!(finitePos(ell.infoEigenvalues.x) && finitePos(ell.infoEigenvalues.y) && finitePos(ell.infoEigenvalues.z))) {
        std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] expected positive finite info eigenvalues.\n";
        ++fails;
      }

      if (!(std::isfinite(ell.conditionNumber) && ell.conditionNumber >= 1.0)) {
        std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] bad condition number. got=" << ell.conditionNumber << "\n";
        ++fails;
      }

      if (!(std::isfinite(ell.sigmaAxisKm.x) && std::isfinite(ell.sigmaAxisKm.y) && std::isfinite(ell.sigmaAxisKm.z))) {
        std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] non-finite sigma axes.\n";
        ++fails;
      }

      if (!(ell.sigmaAxisKm.x >= 0.0 && ell.sigmaAxisKm.y >= 0.0 && ell.sigmaAxisKm.z >= 0.0)) {
        std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] negative sigma axis length(s).\n";
        ++fails;
      }
    }
  }

  // --- Weight scaling: larger weights should shrink the sigma axes ---
  {
    const auto tr1 = makeTrack(1.0);
    const auto tr4 = makeTrack(4.0);

    const auto e1 = sim::bearingTrackUncertaintyEllipsoid(tr1, ep);
    const auto e4 = sim::bearingTrackUncertaintyEllipsoid(tr4, ep);

    if (!(e1.valid && e4.valid)) {
      std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] expected valid ellipsoids for scaling test.\n";
      ++fails;
    } else {
      // Expect roughly 1/sqrt(4) = 0.5 scaling on all axes.
      const auto checkRatio = [&](double a, double b, const char* name) {
        if (!(a > 1.0e-15)) return;
        const double r = b / a;
        if (!(r > 0.35 && r < 0.65)) {
          std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] bad axis ratio for " << name
                    << " r=" << r << " (expected ~0.5)\n";
          ++fails;
        }
      };

      checkRatio(e1.sigmaAxisKm.x, e4.sigmaAxisKm.x, "x");
      checkRatio(e1.sigmaAxisKm.y, e4.sigmaAxisKm.y, "y");
      checkRatio(e1.sigmaAxisKm.z, e4.sigmaAxisKm.z, "z");
    }
  }

  // --- Rank-deficient case: one bearing should produce a very elongated axis ---
  {
    sim::BearingTrack3d tr{};
    const math::Vec3d o{0.0, 0.0, 0.0};
    const math::Vec3d d = (target - o).normalized();
    sim::updateBearingTrack(tr, 0.0, true, o, d, 1.0, tp);

    const auto ell = sim::bearingTrackUncertaintyEllipsoid(tr, ep);
    if (!ell.valid) {
      std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] expected valid ellipsoid for single-bearing case.\n";
      ++fails;
    } else {
      // With minInfoEigenvalue=1e-12 and sigma=1, the weakest axis should be ~1e6 km.
      const double maxAxis = std::max({ell.sigmaAxisKm.x, ell.sigmaAxisKm.y, ell.sigmaAxisKm.z});
      if (!(maxAxis > 1.0e4)) {
        std::cerr << "[test_geometry_bearing_uncertainty_ellipsoid] expected elongated axis for rank deficiency. maxAxis="
                  << maxAxis << "\n";
        ++fails;
      }
    }
  }

  if (fails == 0) {
    std::cout << "[test_geometry_bearing_uncertainty_ellipsoid] PASS\n";
  }

  return fails;
}
