#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>

using namespace stellar;

int test_point_defense_engagement() {
  int fails = 0;

  const sim::CombatTargetKind targetKind = sim::CombatTargetKind::Player;
  const core::u64 targetId = 42;

  const math::Vec3d targetPos{0, 0, 0};
  const math::Vec3d targetVel{0, 0, 0};

  sim::Missile missiles[2]{};

  // Nearest inbound missile.
  missiles[0].posKm = {0, 0, 50};
  missiles[0].velKmS = {0, 0, -10};
  missiles[0].ttlSimSec = 30.0;
  missiles[0].turnRateRadS = 0.0;
  missiles[0].hasTarget = true;
  missiles[0].targetKind = targetKind;
  missiles[0].targetId = targetId;

  // Second missile farther out.
  missiles[1].posKm = {0, 0, 100};
  missiles[1].velKmS = {0, 0, -10};
  missiles[1].ttlSimSec = 30.0;
  missiles[1].turnRateRadS = 0.0;
  missiles[1].hasTarget = true;
  missiles[1].targetKind = targetKind;
  missiles[1].targetId = targetId;

  sim::SphereTarget world[1]{};
  world[0].kind = targetKind;
  world[0].id = targetId;
  world[0].centerKm = targetPos;
  world[0].velKmS = targetVel;
  world[0].radiusKm = 1.0;

  sim::PointDefenseEngagementParams ep{};
  ep.pointDefense.muzzleSpeedKmS = 20.0;
  ep.pointDefense.maxPlanSimSec = 6.0;
  ep.pointDefense.simStepSec = 0.005;
  ep.pointDefense.useFullMissileSimulation = false; // constant-velocity scenario
  ep.maxShots = 2;
  ep.shotCooldownSec = 0.10;
  ep.maxPlanSimSec = 6.0;

  const math::Vec3d turretAimNow{0, 0, 1};

  // --- Case 1: schedule two shots ---
  {
    const auto plan = sim::planPointDefenseEngagement(missiles,
                                                      2,
                                                      targetKind,
                                                      targetId,
                                                      targetPos,
                                                      targetVel,
                                                      turretAimNow,
                                                      world,
                                                      1,
                                                      /*seed=*/1234,
                                                      ep);

    if (!plan.valid) {
      std::cerr << "FAIL: expected valid point defense engagement plan\n";
      ++fails;
    } else {
      if (plan.shots.size() != 2) {
        std::cerr << "FAIL: expected 2 scheduled shots, got " << plan.shots.size() << "\n";
        ++fails;
      } else {
        // Should engage the nearest inbound missile first.
        if (plan.shots[0].threat.missileIndex != 0 || plan.shots[1].threat.missileIndex != 1) {
          std::cerr << "FAIL: expected engagement order [0,1], got ["
                    << plan.shots[0].threat.missileIndex << ", "
                    << plan.shots[1].threat.missileIndex << "]\n";
          ++fails;
        }

        const double f0 = plan.shots[0].solution.fireDelaySec;
        const double f1 = plan.shots[1].solution.fireDelaySec;
        if (f0 < -1e-6 || f1 < -1e-6 || f1 + 1e-9 < f0) {
          std::cerr << "FAIL: fire schedule is not monotonic. f0=" << f0 << " f1=" << f1 << "\n";
          ++fails;
        }

        const double t0 = plan.shots[0].solution.interceptTimeSec;
        const double t1 = plan.shots[1].solution.interceptTimeSec;
        if (t0 <= 0.0 || t1 <= 0.0 || t1 + 1e-9 < t0) {
          std::cerr << "FAIL: intercept schedule is not monotonic. t0=" << t0 << " t1=" << t1 << "\n";
          ++fails;
        }

        // Sanity: first intercept should happen around 1.67s (50km closing @ 30km/s).
        if (t0 < 1.5 || t0 > 2.0) {
          std::cerr << "FAIL: unexpected intercept time for shot0. got=" << t0 << "\n";
          ++fails;
        }
      }
    }
  }

  // --- Case 2: cooldown too large -> only one shot fits in horizon ---
  {
    sim::PointDefenseEngagementParams slow = ep;
    slow.shotCooldownSec = 10.0;

    const auto plan = sim::planPointDefenseEngagement(missiles,
                                                      2,
                                                      targetKind,
                                                      targetId,
                                                      targetPos,
                                                      targetVel,
                                                      turretAimNow,
                                                      world,
                                                      1,
                                                      /*seed=*/5678,
                                                      slow);

    if (!plan.valid) {
      std::cerr << "FAIL: expected valid plan even with large cooldown\n";
      ++fails;
    } else {
      if (plan.shots.size() != 1) {
        std::cerr << "FAIL: expected 1 scheduled shot with large cooldown, got " << plan.shots.size() << "\n";
        ++fails;
      }
    }
  }

  return fails;
}
