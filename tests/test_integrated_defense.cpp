#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

int test_integrated_defense() {
  int fails = 0;

  const sim::CombatTargetKind kind = sim::CombatTargetKind::Ship;
  const core::u64 targetId = 1;

  const math::Vec3d targetPos{10, 0, 0};
  const math::Vec3d targetVel{0, 0, 0};

  sim::SphereTarget a{};
  a.kind = sim::CombatTargetKind::Asteroid;
  a.id = 111;
  a.centerKm = {5, 3, 0};
  a.radiusKm = 2.0;

  const sim::SphereTarget worldTargets[] = {a};

  sim::IntegratedDefenseParams p{};
  p.enableTerrainMasking = true;
  p.minTtiForTerrainSec = 1.0;
  p.terrain.maxCoverTravelKm = 100.0;
  p.terrain.lookaheadSec = 0.0;
  p.terrain.coverPadKm = 0.01;
  p.terrain.ignoreIfAlreadyOccluded = true;

  // Case 1: LOS-dependent inbound missile -> terrain masking should be selected.
  {
    sim::Missile m{};
    m.hasTarget = true;
    m.targetKind = kind;
    m.targetId = targetId;
    m.seeker = sim::MissileSeekerType::Radar;
    m.requireLineOfSight = true;
    m.posKm = {0, 0, 0};
    m.velKmS = {1, 0, 0};
    m.ttlSimSec = 100.0;

    std::vector<sim::Missile> missiles{m};

    const auto plan = sim::planIntegratedDefense(missiles.data(),
                                                 missiles.size(),
                                                 kind,
                                                 targetId,
                                                 targetPos,
                                                 targetVel,
                                                 worldTargets,
                                                 1,
                                                 /*seed=*/7,
                                                 p);

    if (!plan.valid || !plan.threat.inbound) {
      std::cerr << "FAIL: expected valid integrated defense plan (LOS case)\n";
      fails++;
    } else {
      if (plan.maneuver != sim::DefenseManeuverKind::TerrainMask) {
        std::cerr << "FAIL: expected TerrainMask maneuver when LOS-dependent\n";
        fails++;
      }
      if (!plan.terrain.valid) {
        std::cerr << "FAIL: expected valid terrain sub-plan\n";
        fails++;
      }
      if (plan.maneuverDirWorld.lengthSq() < 1e-12) {
        std::cerr << "FAIL: expected non-zero maneuverDirWorld\n";
        fails++;
      } else {
        const double d = math::dot(plan.maneuverDirWorld, plan.terrain.dirWorld);
        if (d < 0.999) {
          std::cerr << "FAIL: expected maneuverDirWorld to match terrain.dirWorld (dot=" << d << ")\n";
          fails++;
        }
      }

      if (!plan.countermeasures.valid) {
        std::cerr << "FAIL: expected countermeasure plan to be valid (LOS case)\n";
        fails++;
      }
      if (plan.countermeasures.type != sim::CountermeasureType::Chaff) {
        std::cerr << "FAIL: expected chaff for radar seeker\n";
        fails++;
      }

      if (!(plan.urgency01 >= 0.0 && plan.urgency01 <= 1.0) || !std::isfinite(plan.urgency01)) {
        std::cerr << "FAIL: urgency01 out of range\n";
        fails++;
      }
    }
  }

  // Case 2: Same geometry but no LOS requirement -> fall back to classic evasion jink.
  {
    sim::Missile m{};
    m.hasTarget = true;
    m.targetKind = kind;
    m.targetId = targetId;
    m.seeker = sim::MissileSeekerType::Radar;
    m.requireLineOfSight = false;
    m.datalinkRequireLineOfSight = false;
    m.posKm = {0, 0, 0};
    m.velKmS = {1, 0, 0};
    m.ttlSimSec = 100.0;

    std::vector<sim::Missile> missiles{m};

    const auto plan = sim::planIntegratedDefense(missiles.data(),
                                                 missiles.size(),
                                                 kind,
                                                 targetId,
                                                 targetPos,
                                                 targetVel,
                                                 worldTargets,
                                                 1,
                                                 /*seed=*/7,
                                                 p);

    if (!plan.valid || !plan.threat.inbound) {
      std::cerr << "FAIL: expected valid integrated defense plan (jink case)\n";
      fails++;
    } else {
      if (plan.maneuver != sim::DefenseManeuverKind::EvasionJink) {
        std::cerr << "FAIL: expected EvasionJink maneuver when terrain is not useful\n";
        fails++;
      }
      if (!plan.evasion.valid) {
        std::cerr << "FAIL: expected valid evasion sub-plan\n";
        fails++;
      }

      const double d = math::dot(plan.maneuverDirWorld, plan.evasion.dirWorld);
      if (d < 0.999) {
        std::cerr << "FAIL: expected maneuverDirWorld to match evasion.dirWorld (dot=" << d << ")\n";
        fails++;
      }

      // Countermeasures should still be planned.
      if (!plan.countermeasures.valid) {
        std::cerr << "FAIL: expected countermeasure plan to be valid (jink case)\n";
        fails++;
      }
      if (!approx(plan.countermeasures.firstReleaseDelaySec, 3.95, 1e-3)) {
        std::cerr << "FAIL: unexpected firstReleaseDelaySec. got=" << plan.countermeasures.firstReleaseDelaySec
                  << " expected~3.95\n";
        fails++;
      }
    }
  }

  return fails;
}
