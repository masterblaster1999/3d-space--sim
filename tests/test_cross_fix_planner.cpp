#include "stellar/sim/CrossFixPlanner.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

static bool approxVec3(const math::Vec3d& a, const math::Vec3d& b, double eps = 1e-6) {
  return approx(a.x, b.x, eps) && approx(a.y, b.y, eps) && approx(a.z, b.z, eps);
}

int test_cross_fix_planner() {
  int fails = 0;

  // --- No fix: baseline derived from range guess ---
  {
    sim::CrossFixParams p{};
    p.baselineFromSigma = true;
    p.baselineFromRange = true;
    p.baselineFromRangeFrac = 0.10;
    p.baselineKm = 7777.0;
    p.minBaselineKm = 10.0;
    p.maxBaselineKm = 1.0e9;
    p.seed = 1234;

    const math::Vec3d obs{0.0, 0.0, 0.0};
    const math::Vec3d bearing{0.0, 0.0, 1.0};
    const double rangeGuessKm = 100000.0;

    const auto rec = sim::recommendCrossFixWaypoint((const sim::BearingTrack3d*)nullptr,
                                                    obs,
                                                    bearing,
                                                    rangeGuessKm,
                                                    p);

    if (!rec.valid) {
      std::cerr << "[test_cross_fix_planner] expected valid recommendation for range-guess case.\n";
      ++fails;
    } else {
      const double expectedBaseline = rangeGuessKm * 0.10;
      if (!approx(rec.baselineKm, expectedBaseline, 1e-6)) {
        std::cerr << "[test_cross_fix_planner] bad baseline. got=" << rec.baselineKm
                  << " expected=" << expectedBaseline << "\n";
        ++fails;
      }

      const math::Vec3d bUnit = bearing.normalized();
      const double dotPerp = std::abs(math::dot(rec.preferred.moveDirUnit, bUnit));
      if (dotPerp > 1e-6) {
        std::cerr << "[test_cross_fix_planner] preferred move dir not perpendicular to bearing. dot=" << dotPerp << "\n";
        ++fails;
      }

      const double dist = (rec.preferred.waypointKm - obs).length();
      if (!approx(dist, rec.baselineKm, 1e-6)) {
        std::cerr << "[test_cross_fix_planner] waypoint distance mismatch. got=" << dist
                  << " baseline=" << rec.baselineKm << "\n";
        ++fails;
      }

      // Alternate should be the opposite direction.
      const math::Vec3d sum = rec.preferred.moveDirUnit + rec.alternate.moveDirUnit;
      if (sum.length() > 1e-6) {
        std::cerr << "[test_cross_fix_planner] alternate direction not opposite of preferred.\n";
        ++fails;
      }

      if (rec.preferred.sideSign != -rec.alternate.sideSign) {
        std::cerr << "[test_cross_fix_planner] sideSign mismatch between preferred and alternate.\n";
        ++fails;
      }

      if (!(rec.expectedBearingChangeDeg > 0.0)) {
        std::cerr << "[test_cross_fix_planner] expected non-zero bearing change estimate.\n";
        ++fails;
      }
    }
  }

  // --- With fix: baseline derived from sigma ---
  {
    sim::BearingTrack3d tr{};
    tr.initialized = true;
    tr.posKm = {0.0, 0.0, 100000.0};
    tr.sigmaKm = 5000.0;

    sim::CrossFixParams p{};
    p.baselineFromSigma = true;
    p.baselineFromSigmaMult = 2.0;
    p.baselineFromRange = false;
    p.seed = 42;

    const math::Vec3d obs{0.0, 0.0, 0.0};
    const math::Vec3d bearing = tr.posKm - obs;

    const auto rec = sim::recommendCrossFixWaypoint(&tr, obs, bearing, -1.0, p);

    if (!rec.valid || !rec.trackHasFix) {
      std::cerr << "[test_cross_fix_planner] expected valid recommendation using track fix.\n";
      ++fails;
    } else {
      const double expectedBaseline = tr.sigmaKm * 2.0;
      if (!approx(rec.baselineKm, expectedBaseline, 1e-6)) {
        std::cerr << "[test_cross_fix_planner] bad sigma baseline. got=" << rec.baselineKm
                  << " expected=" << expectedBaseline << "\n";
        ++fails;
      }

      const math::Vec3d rUnit = (tr.posKm - obs).normalized();
      const double dotPerp = std::abs(math::dot(rec.preferred.moveDirUnit, rUnit));
      if (dotPerp > 1e-6) {
        std::cerr << "[test_cross_fix_planner] preferred move dir not perpendicular to LOS. dot=" << dotPerp << "\n";
        ++fails;
      }

      if (!(rec.expectedBearingChangeDeg > 0.0)) {
        std::cerr << "[test_cross_fix_planner] expected non-zero bearing change estimate for fixed target.\n";
        ++fails;
      }

      if (!approxVec3(rec.aimPointKm, tr.posKm, 1e-9)) {
        std::cerr << "[test_cross_fix_planner] expected aim point to equal track position when fixed.\n";
        ++fails;
      }
    }
  }

  // --- Determinism: repeated calls with same seed must match ---
  {
    sim::CrossFixParams p{};
    p.seed = 999;

    const math::Vec3d obs{1.0, 2.0, 3.0};
    const math::Vec3d bearing{0.2, 0.1, 1.0};
    const double rangeGuessKm = 200000.0;

    const auto a = sim::recommendCrossFixWaypoint((const sim::BearingTrack3d*)nullptr,
                                                  obs,
                                                  bearing,
                                                  rangeGuessKm,
                                                  p);
    const auto b = sim::recommendCrossFixWaypoint((const sim::BearingTrack3d*)nullptr,
                                                  obs,
                                                  bearing,
                                                  rangeGuessKm,
                                                  p);

    if (!a.valid || !b.valid) {
      std::cerr << "[test_cross_fix_planner] expected valid recommendations for determinism case.\n";
      ++fails;
    } else {
      if (a.preferred.sideSign != b.preferred.sideSign) {
        std::cerr << "[test_cross_fix_planner] side sign changed between repeated calls.\n";
        ++fails;
      }
      if (!approxVec3(a.preferred.moveDirUnit, b.preferred.moveDirUnit, 1e-9)) {
        std::cerr << "[test_cross_fix_planner] moveDir changed between repeated calls.\n";
        ++fails;
      }
    }
  }



  // --- Prediction: side + baseline optimization by info gain ---
  {
    sim::BearingTrack3d tr{};
    sim::BearingTrackParams tp{};
    tp.observationHalfLifeSec = 0.0;
    tp.solveRegularization = 0.0;
    tp.determinantEps = 1.0e-14;
    tp.minEffectiveWeight = 1.8;

    const math::Vec3d target{0.0, 0.0, 100000.0};

    // Seed the track with one prior bearing from a non-symmetric location.
    const math::Vec3d o0{5000.0, 2000.0, 0.0};
    const math::Vec3d d0 = (target - o0).normalized();
    sim::updateBearingTrack(tr, 0.0, true, o0, d0, 1.0, tp);

    sim::CrossFixParams p{};
    p.baselineFromSigma = false;
    p.baselineFromRange = true;
    p.baselineFromRangeFrac = 0.05;
    p.baselineKm = 1000.0;
    p.minBaselineKm = 1000.0;
    p.maxBaselineKm = 20000.0;
    p.seed = 123;

    p.includePrediction = true;
    p.chooseSideByPrediction = true;
    p.optimizeBaselineByPrediction = true;
    p.predictedMeasurementWeight = 1.0;
    p.baselineSearchSteps = 5;
    p.baselineSearchMinMult = 0.5;
    p.baselineSearchMaxMult = 4.0;
    p.predictedSolveParams = tp;

    const math::Vec3d obs{0.0, 0.0, 0.0};
    const math::Vec3d bearing = target - obs;
    const double rangeGuessKm = bearing.length();

    const auto rec = sim::recommendCrossFixWaypoint(&tr, obs, bearing, rangeGuessKm, p);

    if (!rec.valid) {
      std::cerr << "[test_cross_fix_planner] expected valid recommendation for prediction case.\n";
      ++fails;
    } else {
      if (!rec.hasPrediction) {
        std::cerr << "[test_cross_fix_planner] expected prediction fields to be populated.\n";
        ++fails;
      }
      if (!(rec.preferred.hasPrediction && rec.alternate.hasPrediction)) {
        std::cerr << "[test_cross_fix_planner] expected both candidates to have prediction data.\n";
        ++fails;
      }

      // Preferred should not be worse than alternate under predicted score.
      const double sPref = sim::predictedScore(rec.preferred.predictedSolve);
      const double sAlt = sim::predictedScore(rec.alternate.predictedSolve);
      if (sPref + 1.0e-12 < sAlt) {
        std::cerr << "[test_cross_fix_planner] expected preferred predicted score >= alternate.\n";
        ++fails;
      }

      // In this geometry, bigger baseline should yield better determinant margin.
      // The optimizer should therefore prefer the max clamp.
      if (!approx(rec.baselineKm, p.maxBaselineKm, 1e-6)) {
        std::cerr << "[test_cross_fix_planner] expected optimized baseline to hit max clamp. got=" << rec.baselineKm
                  << " expected=" << p.maxBaselineKm << "\n";
        ++fails;
      }
    }
  }
  if (fails == 0) {
    std::cout << "[test_cross_fix_planner] PASS\n";
  }

  return fails;
}
