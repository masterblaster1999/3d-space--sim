#include "stellar/sim/BearingTrack.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_bearing_track_diagnostics() {
  int fails = 0;

  sim::BearingTrack3d tr{};
  sim::BearingTrackParams p{};
  p.observationHalfLifeSec = 0.0; // deterministic for this test
  p.solveRegularization = 0.0;
  p.determinantEps = 1.0e-14;
  p.minEffectiveWeight = 1.8;

  // Before any measurements.
  {
    const auto d = sim::bearingTrackSolveDiagnostics(tr, p);
    if (d.canSolve) {
      std::cerr << "[test_bearing_track_diagnostics] empty track should not be solvable.\n";
      ++fails;
    }
    if (!approx(d.effectiveWeight, 0.0)) {
      std::cerr << "[test_bearing_track_diagnostics] expected zero weight.\n";
      ++fails;
    }
    if (!(d.progress01 <= 1.0e-9)) {
      std::cerr << "[test_bearing_track_diagnostics] expected near-zero progress for empty track.\n";
      ++fails;
    }
  }

  const math::Vec3d target{100.0, -50.0, 20.0};

  // One bearing: weight increases but determinant stays ~0.
  {
    const math::Vec3d o0{0.0, 0.0, 0.0};
    const math::Vec3d d0 = (target - o0).normalized();
    sim::updateBearingTrack(tr, 0.0, true, o0, d0, 1.0, p);

    const auto d = sim::bearingTrackSolveDiagnostics(tr, p);
    if (d.canSolve) {
      std::cerr << "[test_bearing_track_diagnostics] single bearing should not be solvable.\n";
      ++fails;
    }
    if (!(d.weightProgress01 > 0.0)) {
      std::cerr << "[test_bearing_track_diagnostics] expected weight progress > 0 after a measurement.\n";
      ++fails;
    }
    if (!(d.detAbs <= d.detThresh + 1e-12)) {
      std::cerr << "[test_bearing_track_diagnostics] expected determinant to be below threshold for single bearing.\n";
      ++fails;
    }
  }

  // Second bearing with a baseline should become solvable.
  {
    const math::Vec3d o1{50.0, 20.0, -10.0};
    const math::Vec3d d1 = (target - o1).normalized();
    sim::updateBearingTrack(tr, 1.0, true, o1, d1, 1.0, p);

    const auto d = sim::bearingTrackSolveDiagnostics(tr, p);
    if (!d.canSolve) {
      std::cerr << "[test_bearing_track_diagnostics] expected canSolve after two bearings with baseline.\n";
      ++fails;
    }
    if (!(d.progress01 >= 0.999)) {
      std::cerr << "[test_bearing_track_diagnostics] expected progress near 1 after solvable geometry. progress=" << d.progress01
                << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_bearing_track_diagnostics] PASS\n";
  }

  return fails;
}
