#include "stellar/sim/Combat.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

int test_terminal_weave() {
  int fails = 0;

  // Simple stationary target far away so the missile doesn't impact during the test.
  sim::SphereTarget tgt{};
  tgt.kind = sim::CombatTargetKind::Ship;
  tgt.index = 0;
  tgt.id = 1;
  tgt.centerKm = {1000, 0, 0};
  tgt.velKmS = {0, 0, 0};
  tgt.radiusKm = 1.0;

  sim::Missile m{};
  m.prevKm = {0, 0, 0};
  m.posKm = {0, 0, 0};
  m.velKmS = {10, 0, 0};
  m.ttlSimSec = 10.0;
  m.turnRateRadS = 100.0;

  // Ensure the seeker is active so terminal weaving is enabled.
  m.seeker = sim::MissileSeekerType::Radar;
  m.seekerActivationSimSec = 0.0;
  m.seekerFovCos = -1.0;

  m.hasTarget = true;
  m.targetKind = sim::CombatTargetKind::Ship;
  m.targetId = 1;

  // Enable weave immediately with a known phase.
  m.terminalWeaveAmplitudeRad = 0.25;
  m.terminalWeaveFrequencyHz = 1.0;
  m.terminalWeaveStartTtiSec = 0.0;
  m.terminalWeaveMaxRad = 0.0;

  m.hasTerminalWeavePhase = true;
  m.terminalWeavePhaseRad = 0.0;

  std::vector<sim::Missile> missiles{m};
  std::vector<sim::MissileDetonation> dets;
  std::vector<sim::MissileHit> hits;

  // Step once: phase == 0 -> lateral basis points along -Z in our chosen basis,
  // so the missile should pick up a negative Z component.
  {
    sim::SphereTarget targets[1]{tgt};
    sim::stepMissiles(missiles, 0.05, targets, 1, dets, hits);

    if (missiles.empty()) {
      std::cerr << "FAIL: missile destroyed unexpectedly (first step)\n";
      return 1;
    }

    const double z = missiles[0].velKmS.z;
    if (!(z < -1e-6)) {
      std::cerr << "FAIL: expected negative Z velocity from weave at phase=0, got z=" << z << "\n";
      ++fails;
    }
  }

  // Advance to just after t=0.25s; the step starting at t=0.25 uses phase=pi/2,
  // biasing toward +Y, so Y velocity should become positive.
  double t = 0.05;
  while (t < 0.30 - 1e-12) {
    sim::SphereTarget targets[1]{tgt};
    sim::stepMissiles(missiles, 0.05, targets, 1, dets, hits);
    t += 0.05;

    if (missiles.empty()) {
      std::cerr << "FAIL: missile destroyed unexpectedly while advancing\n";
      return 1;
    }
  }

  const double yEnd = missiles[0].velKmS.y;
  if (!(yEnd > 1e-6)) {
    std::cerr << "FAIL: expected positive Y velocity after 0.25s weave phase, got y=" << yEnd << "\n";
    ++fails;
  }

  // Determinism check: repeat with identical initial conditions and compare end state.
  {
    sim::Missile m2 = m;
    std::vector<sim::Missile> ms2{m2};
    std::vector<sim::MissileDetonation> dets2;
    std::vector<sim::MissileHit> hits2;

    double t2 = 0.0;
    while (t2 < 0.30 - 1e-12) {
      sim::SphereTarget targets[1]{tgt};
      sim::stepMissiles(ms2, 0.05, targets, 1, dets2, hits2);
      t2 += 0.05;
      if (ms2.empty()) break;
    }

    if (ms2.empty()) {
      std::cerr << "FAIL: second run missile destroyed unexpectedly\n";
      ++fails;
    } else {
      const math::Vec3d dv = ms2[0].velKmS - missiles[0].velKmS;
      if (dv.lengthSq() > 1e-12) {
        std::cerr << "FAIL: weave simulation not deterministic; dv=" << dv.x << "," << dv.y << "," << dv.z << "\n";
        ++fails;
      }
    }
  }

  return fails;
}
