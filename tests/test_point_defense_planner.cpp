#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>

using namespace stellar;


int test_point_defense_planner() {
  int fails = 0;

  const math::Vec3d shooterPos{0, 0, 0};
  const math::Vec3d shooterVel{0, 0, 0};

  sim::Missile m{};
  m.posKm = {0, 0, 50};
  m.velKmS = {0, 0, -10};
  m.ttlSimSec = 10.0;
  m.turnRateRadS = 0.0;
  m.hasTarget = false;

  sim::PointDefenseParams p{};
  p.muzzleSpeedKmS = 20.0;
  p.maxPlanSimSec = 6.0;
  p.simStepSec = 0.005;
  p.useFullMissileSimulation = false; // constant-velocity scenario

  // --- Case 1: unlimited slew, already pointed roughly at the intercept ---
  {
    const math::Vec3d aimNow{0, 0, 1};
    p.turretSlewRateRadS = 0.0;

    const auto sol = sim::planPointDefenseShot(m,
                                              shooterPos,
                                              shooterVel,
                                              aimNow,
                                              /*worldTargets=*/nullptr,
                                              /*worldTargetCount=*/0,
                                              /*seed=*/123,
                                              p);

    if (!sol.valid) {
      std::cerr << "FAIL: expected valid point-defense solution (immediate fire)\n";
      ++fails;
    } else {
      if (sol.fireDelaySec > 0.05) {
        std::cerr << "FAIL: expected near-immediate fireDelaySec, got " << sol.fireDelaySec << "\n";
        ++fails;
      }
      if (sol.interceptTimeSec < 1.5 || sol.interceptTimeSec > 2.0) {
        std::cerr << "FAIL: interceptTimeSec out of expected range, got " << sol.interceptTimeSec << "\n";
        ++fails;
      }
      const double aimDot = math::dot(sol.aimDirWorld.normalized(), math::Vec3d{0, 0, 1});
      if (aimDot < 0.999) {
        std::cerr << "FAIL: expected aimDirWorld roughly +Z, dot=" << aimDot << "\n";
        ++fails;
      }
    }
  }

  // --- Case 2: slew-limited turret must wait to rotate before firing ---
  {
    const math::Vec3d aimNow{1, 0, 0};
    p.turretSlewRateRadS = 0.5;

    const auto sol = sim::planPointDefenseShot(m,
                                              shooterPos,
                                              shooterVel,
                                              aimNow,
                                              /*worldTargets=*/nullptr,
                                              /*worldTargetCount=*/0,
                                              /*seed=*/999,
                                              p);

    if (!sol.valid) {
      std::cerr << "FAIL: expected valid point-defense solution (slew-limited)\n";
      ++fails;
    } else {
      // ~90deg / 0.5rad/s = 3.14s required slew before firing.
      if (sol.fireDelaySec < 3.0 || sol.fireDelaySec > 3.35) {
        std::cerr << "FAIL: fireDelaySec not consistent with slew constraint. got=" << sol.fireDelaySec << "\n";
        ++fails;
      }
      if (sol.interceptTimeSec < 3.55 || sol.interceptTimeSec > 3.95) {
        std::cerr << "FAIL: interceptTimeSec not in expected range. got=" << sol.interceptTimeSec << "\n";
        ++fails;
      }
      if (sol.requiredSlewTimeSec < 3.0 || sol.requiredSlewTimeSec > 3.35) {
        std::cerr << "FAIL: requiredSlewTimeSec not in expected range. got=" << sol.requiredSlewTimeSec << "\n";
        ++fails;
      }

      const double aimDot = math::dot(sol.aimDirWorld.normalized(), math::Vec3d{0, 0, 1});
      if (aimDot < 0.999) {
        std::cerr << "FAIL: expected aimDirWorld roughly +Z after slew, dot=" << aimDot << "\n";
        ++fails;
      }
    }
  }

  // --- Case 3: turret too slow -> no feasible intercept before impact ---
  {
    const math::Vec3d aimNow{1, 0, 0};
    p.turretSlewRateRadS = 0.1; // requires ~15.7s to rotate 90deg

    const auto sol = sim::planPointDefenseShot(m,
                                              shooterPos,
                                              shooterVel,
                                              aimNow,
                                              /*worldTargets=*/nullptr,
                                              /*worldTargetCount=*/0,
                                              /*seed=*/555,
                                              p);

    if (sol.valid) {
      std::cerr << "FAIL: expected invalid solution (turret too slow)\n";
      ++fails;
    }
  }

  return fails;
}
