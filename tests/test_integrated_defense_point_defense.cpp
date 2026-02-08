#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

int test_integrated_defense_point_defense() {
  int fails = 0;

  const sim::CombatTargetKind kind = sim::CombatTargetKind::Ship;
  const core::u64 targetId = 42;

  const math::Vec3d targetPos{0, 0, 0};
  const math::Vec3d targetVel{0, 0, 0};

  sim::Missile m{};
  m.hasTarget = true;
  m.targetKind = kind;
  m.targetId = targetId;
  m.seeker = sim::MissileSeekerType::Heat;
  m.posKm = {0, 0, 50};
  m.velKmS = {0, 0, -10};
  m.ttlSimSec = 20.0;
  m.turnRateRadS = 0.0;

  std::vector<sim::Missile> missiles{m};

  sim::IntegratedDefenseWithPointDefenseParams p{};
  p.enablePointDefense = true;
  p.maxTtiForPointDefenseSec = 10.0;
  p.interceptBeforeImpactMarginSec = 0.01;
  p.clampPointDefenseHorizonToTti = true;

  p.pointDefense.muzzleSpeedKmS = 20.0;
  p.pointDefense.maxPlanSimSec = 6.0;
  p.pointDefense.simStepSec = 0.005;
  p.pointDefense.useFullMissileSimulation = false;  // constant-velocity scenario
  p.pointDefense.mountConeCos = -1.0;
  p.pointDefense.turretSlewRateRadS = 0.0;
  p.pointDefense.requireLineOfSight = false;

  // Case 1: Point defense is feasible -> expect a valid PD shot suggestion.
  {
    const math::Vec3d aimNow{0, 0, 1};

    const auto plan = sim::planIntegratedDefenseWithPointDefense(missiles.data(),
                                                                 missiles.size(),
                                                                 kind,
                                                                 targetId,
                                                                 targetPos,
                                                                 targetVel,
                                                                 aimNow,
                                                                 /*worldTargets=*/nullptr,
                                                                 /*worldTargetCount=*/0,
                                                                 /*seed=*/123,
                                                                 p);

    if (!plan.valid || !plan.defense.valid || !plan.defense.threat.inbound) {
      std::cerr << "FAIL: expected valid integrated defense plan (PD case)\n";
      ++fails;
    } else {
      // Threat time-to-impact in this setup should be ~5s.
      if (!approx(plan.defense.threat.ttiSec, 5.0, 1e-6)) {
        std::cerr << "FAIL: unexpected ttiSec. got=" << plan.defense.threat.ttiSec << " expected=5\n";
        ++fails;
      }

      if (!plan.pointDefense.valid) {
        std::cerr << "FAIL: expected valid point-defense solution\n";
        ++fails;
      } else {
        if (plan.pointDefense.fireDelaySec > 0.05) {
          std::cerr << "FAIL: expected near-immediate fireDelaySec. got=" << plan.pointDefense.fireDelaySec << "\n";
          ++fails;
        }

        // Solve 20*t = 50 - 10*t => t = 50/30 = 1.666...
        if (plan.pointDefense.interceptTimeSec < 1.5 || plan.pointDefense.interceptTimeSec > 1.9) {
          std::cerr << "FAIL: interceptTimeSec out of expected range. got=" << plan.pointDefense.interceptTimeSec
                    << " expected~1.67\n";
          ++fails;
        }

        const double aimDot = math::dot(plan.pointDefense.aimDirWorld.normalized(), math::Vec3d{0, 0, 1});
        if (aimDot < 0.999) {
          std::cerr << "FAIL: expected aimDirWorld roughly +Z. dot=" << aimDot << "\n";
          ++fails;
        }
      }

      // Defense maneuver is still recommended (PD is additive).
      if (plan.defense.maneuver == sim::DefenseManeuverKind::None) {
        std::cerr << "FAIL: expected a non-None defensive maneuver suggestion\n";
        ++fails;
      }

      if (!plan.defense.countermeasures.valid) {
        std::cerr << "FAIL: expected countermeasure response plan to be valid\n";
        ++fails;
      }
    }
  }

  // Case 2: Turret too slow -> integrated plan remains valid but PD should be invalid.
  {
    sim::IntegratedDefenseWithPointDefenseParams p2 = p;
    p2.pointDefense.turretSlewRateRadS = 0.1;  // ~15.7s for 90deg turn

    const math::Vec3d aimNow{1, 0, 0};

    const auto plan = sim::planIntegratedDefenseWithPointDefense(missiles.data(),
                                                                 missiles.size(),
                                                                 kind,
                                                                 targetId,
                                                                 targetPos,
                                                                 targetVel,
                                                                 aimNow,
                                                                 /*worldTargets=*/nullptr,
                                                                 /*worldTargetCount=*/0,
                                                                 /*seed=*/999,
                                                                 p2);

    if (!plan.valid || !plan.defense.valid || !plan.defense.threat.inbound) {
      std::cerr << "FAIL: expected valid integrated defense plan (slow turret case)\n";
      ++fails;
    } else {
      if (plan.pointDefense.valid) {
        std::cerr << "FAIL: expected invalid point-defense solution when turret is too slow\n";
        ++fails;
      }
    }
  }

  // Case 3: Time-to-impact gate disables PD planning.
  {
    sim::IntegratedDefenseWithPointDefenseParams p3 = p;
    p3.maxTtiForPointDefenseSec = 1.0;  // tti ~ 5s, so skip

    const math::Vec3d aimNow{0, 0, 1};

    const auto plan = sim::planIntegratedDefenseWithPointDefense(missiles.data(),
                                                                 missiles.size(),
                                                                 kind,
                                                                 targetId,
                                                                 targetPos,
                                                                 targetVel,
                                                                 aimNow,
                                                                 /*worldTargets=*/nullptr,
                                                                 /*worldTargetCount=*/0,
                                                                 /*seed=*/1,
                                                                 p3);

    if (!plan.valid || !plan.defense.valid || !plan.defense.threat.inbound) {
      std::cerr << "FAIL: expected valid integrated defense plan (tti gate case)\n";
      ++fails;
    } else {
      if (plan.pointDefense.valid) {
        std::cerr << "FAIL: expected point-defense planning to be skipped by maxTti gate\n";
        ++fails;
      }
    }
  }

  return fails;
}
