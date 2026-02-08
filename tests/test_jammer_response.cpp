#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_jammer_response() {
  int fails = 0;

  const math::Vec3d targetPos{0, 0, 10};
  const math::Vec3d targetVel{0, 0, 0};

  // Case 1: Radar missile without HOJ -> jammer should be recommended when inbound and within window.
  {
    sim::Missile m{};
    m.seeker = sim::MissileSeekerType::Radar;
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 5};
    m.ttlSimSec = 100.0;
    m.homeOnJam = false;

    const auto plan = sim::planJammerResponse(m, targetPos, targetVel);
    if (!plan.valid) {
      std::cerr << "FAIL: expected valid jammer plan for radar threat\n";
      fails++;
    } else {
      if (!plan.jammerOn) {
        std::cerr << "FAIL: expected jammerOn for non-HOJ radar threat\n";
        fails++;
      }
      if (!approx(plan.jammerPower, 1.0, 1e-12)) {
        std::cerr << "FAIL: expected jammerPower=1.0, got " << plan.jammerPower << "\n";
        fails++;
      }
      if (!(plan.receivedJamming01 > 0.9)) {
        std::cerr << "FAIL: expected high receivedJamming01, got " << plan.receivedJamming01 << "\n";
        fails++;
      }
    }
  }

  // Case 2: HOJ-capable missile near doppler notch -> recommend jammer OFF to deny HOJ.
  {
    sim::Missile m{};
    m.seeker = sim::MissileSeekerType::Radar;
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 0.2};
    m.ttlSimSec = 100.0;

    m.homeOnJam = true;
    m.homeOnJamMinJammerPower = 0.5;
    m.homeOnJamMinJamming01 = 0.2;
    m.radarDopplerNotchKmS = 1.0;

    const math::Vec3d closeTarget{0, 0, 1};

    const auto plan = sim::planJammerResponse(m, closeTarget, targetVel);
    if (!plan.valid) {
      std::cerr << "FAIL: expected valid jammer plan (HOJ notch case)\n";
      fails++;
    } else {
      if (!plan.homeOnJamCapable) {
        std::cerr << "FAIL: expected homeOnJamCapable true\n";
        fails++;
      }
      if (plan.jammerOn) {
        std::cerr << "FAIL: expected jammerOff to deny HOJ when near notch\n";
        fails++;
      }
      if (plan.reason != sim::JammerResponseReason::DenyHomeOnJamNotch) {
        std::cerr << "FAIL: expected DenyHomeOnJamNotch reason\n";
        fails++;
      }
    }
  }

  // Case 3: HOJ-capable missile extremely close -> recommend jammer OFF (deny close-range HOJ cue).
  {
    sim::Missile m{};
    m.seeker = sim::MissileSeekerType::Radar;
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 1.0};
    m.ttlSimSec = 100.0;

    m.homeOnJam = true;
    m.homeOnJamMinJammerPower = 0.5;
    m.homeOnJamMinJamming01 = 0.2;
    m.radarDopplerNotchKmS = 0.05;

    sim::JammerResponseParams jp{};
    jp.denyHomeOnJamTtiSec = 2.0;

    const math::Vec3d closeTarget{0, 0, 1};

    const auto plan = sim::planJammerResponse(m, closeTarget, targetVel, jp);
    if (!plan.valid) {
      std::cerr << "FAIL: expected valid jammer plan (HOJ close case)\n";
      fails++;
    } else {
      if (!plan.homeOnJamCapable) {
        std::cerr << "FAIL: expected homeOnJamCapable true (HOJ close case)\n";
        fails++;
      }
      if (plan.jammerOn) {
        std::cerr << "FAIL: expected jammerOff to deny HOJ at close range\n";
        fails++;
      }
      if (plan.reason != sim::JammerResponseReason::DenyHomeOnJamClose) {
        std::cerr << "FAIL: expected DenyHomeOnJamClose reason\n";
        fails++;
      }
    }
  }

  // Case 4: Too-far gating -> jammer not recommended if tti is beyond maxTtiToJamSec.
  {
    sim::Missile m{};
    m.seeker = sim::MissileSeekerType::Radar;
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 1.0};
    m.ttlSimSec = 10000.0;

    const math::Vec3d farTarget{0, 0, 5000};

    sim::JammerResponseParams jp{};
    jp.maxTtiToJamSec = 60.0;

    const auto plan = sim::planJammerResponse(m, farTarget, targetVel, jp);
    if (!plan.valid) {
      std::cerr << "FAIL: expected valid jammer plan (far case)\n";
      fails++;
    } else {
      if (plan.jammerOn) {
        std::cerr << "FAIL: expected jammerOff when too far\n";
        fails++;
      }
      if (plan.reason != sim::JammerResponseReason::TooFar) {
        std::cerr << "FAIL: expected TooFar reason\n";
        fails++;
      }
    }
  }

  // Case 5: Integrated defense terrain masking + HOJ -> override recommends jammer OFF.
  {
    const sim::CombatTargetKind kind = sim::CombatTargetKind::Ship;
    const core::u64 targetId = 1;

    const math::Vec3d tPos{10, 0, 0};
    const math::Vec3d tVel{0, 0, 0};

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

    sim::Missile m{};
    m.hasTarget = true;
    m.targetKind = kind;
    m.targetId = targetId;
    m.seeker = sim::MissileSeekerType::Radar;
    m.requireLineOfSight = true;
    m.posKm = {0, 0, 0};
    m.velKmS = {1, 0, 0};
    m.ttlSimSec = 100.0;

    m.homeOnJam = true;
    m.homeOnJamMinJammerPower = 0.5;
    m.homeOnJamMinJamming01 = 0.2;
    m.radarDopplerNotchKmS = 0.2;

    std::vector<sim::Missile> missiles{m};

    const auto plan = sim::planIntegratedDefense(missiles.data(),
                                                 missiles.size(),
                                                 kind,
                                                 targetId,
                                                 tPos,
                                                 tVel,
                                                 worldTargets,
                                                 1,
                                                 /*seed=*/7,
                                                 p);

    if (!plan.valid || !plan.threat.inbound) {
      std::cerr << "FAIL: expected valid integrated defense plan (HOJ terrain case)\n";
      fails++;
    } else {
      if (plan.maneuver != sim::DefenseManeuverKind::TerrainMask) {
        std::cerr << "FAIL: expected TerrainMask maneuver (HOJ terrain case)\n";
        fails++;
      }

      if (!plan.jammer.valid) {
        std::cerr << "FAIL: expected jammer plan valid in integrated defense\n";
        fails++;
      } else {
        if (plan.jammer.jammerOn) {
          std::cerr << "FAIL: expected jammerOff when terrain masking against HOJ-capable missile\n";
          fails++;
        }
        if (plan.jammer.reason != sim::JammerResponseReason::DenyHomeOnJamTerrainMask) {
          std::cerr << "FAIL: expected DenyHomeOnJamTerrainMask reason\n";
          fails++;
        }
      }
    }
  }

  return fails;
}
