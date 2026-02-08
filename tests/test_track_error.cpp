#include "stellar/sim/Combat.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

int test_track_error() {
  int fails = 0;

  // Stationary target far away so the missile does not impact during the test.
  sim::SphereTarget tgt{};
  tgt.kind = sim::CombatTargetKind::Ship;
  tgt.index = 0;
  tgt.id = 1;
  tgt.centerKm = {1000, 0, 0};
  tgt.velKmS = {0, 0, 0};
  tgt.radiusKm = 1.0;

  // --- Case 1: poor track quality -> track error should bias the commanded direction. ---
  sim::Missile m{};
  m.prevKm = {0, 0, 0};
  m.posKm = {0, 0, 0};
  m.velKmS = {10, 0, 0};
  m.ttlSimSec = 10.0;
  m.turnRateRadS = 100.0;

  // Ensure the seeker is active.
  m.seeker = sim::MissileSeekerType::Radar;
  m.seekerActivationSimSec = 0.0;
  m.seekerFovCos = -1.0;

  m.hasTarget = true;
  m.targetKind = sim::CombatTargetKind::Ship;
  m.targetId = 1;

  // Force track quality to remain very low.
  m.enableTrackQuality = true;
  m.trackQuality = 0.0;
  m.trackQualityRiseHalfLifeSimSec = 1e9;
  m.trackQualityFallHalfLifeSimSec = 1e9;

  // Enable deterministic track error with a known phase.
  m.trackErrorMaxRad = 0.25;
  m.trackErrorFrequencyHz = 1.0;
  m.hasTrackErrorPhase = true;
  m.trackErrorPhaseRad = 0.0;

  std::vector<sim::Missile> missiles{m};
  std::vector<sim::MissileDetonation> dets;
  std::vector<sim::MissileHit> hits;

  {
    sim::SphereTarget targets[1]{tgt};
    sim::stepMissiles(missiles, 0.05, targets, 1, dets, hits);

    if (missiles.empty()) {
      std::cerr << "FAIL: missile destroyed unexpectedly (track error case)\\n";
      return 1;
    }

    // With forward along +X, our basis uses right == -Z. Phase==0 biases toward right, so Z < 0.
    const double z = missiles[0].velKmS.z;
    if (!(z < -1e-6)) {
      std::cerr << "FAIL: expected negative Z velocity from track error at phase=0, got z=" << z << "\\n";
      ++fails;
    }
  }

  // --- Case 2: perfect track quality -> error amplitude should be zero. ---
  {
    sim::Missile m2 = m;
    m2.trackQuality = 1.0;

    std::vector<sim::Missile> ms2{m2};
    std::vector<sim::MissileDetonation> dets2;
    std::vector<sim::MissileHit> hits2;

    sim::SphereTarget targets[1]{tgt};
    sim::stepMissiles(ms2, 0.05, targets, 1, dets2, hits2);

    if (ms2.empty()) {
      std::cerr << "FAIL: missile destroyed unexpectedly (perfect quality case)\\n";
      ++fails;
    } else {
      const double z2 = ms2[0].velKmS.z;
      if (!(std::fabs(z2) < 1e-9)) {
        std::cerr << "FAIL: expected ~0 Z velocity when trackQuality==1, got z=" << z2 << "\\n";
        ++fails;
      }
    }
  }

  // Determinism check: repeat Case 1 and compare end state.
  {
    sim::Missile m3 = m;
    std::vector<sim::Missile> ms3{m3};
    std::vector<sim::MissileDetonation> dets3;
    std::vector<sim::MissileHit> hits3;

    sim::SphereTarget targets[1]{tgt};
    sim::stepMissiles(ms3, 0.05, targets, 1, dets3, hits3);

    if (ms3.empty()) {
      std::cerr << "FAIL: determinism run missile destroyed unexpectedly\\n";
      ++fails;
    } else if (!missiles.empty()) {
      const math::Vec3d dv = ms3[0].velKmS - missiles[0].velKmS;
      if (dv.lengthSq() > 1e-12) {
        std::cerr << "FAIL: track error simulation not deterministic; dv=" << dv.x << "," << dv.y << "," << dv.z
                  << "\\n";
        ++fails;
      }
    }
  }

  return fails;
}
