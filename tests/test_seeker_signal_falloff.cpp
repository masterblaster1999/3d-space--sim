#include "stellar/sim/Combat.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

int test_seeker_signal_falloff() {
  int fails = 0;

  auto makeMissile = []() {
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {10, 0, 0};
    m.ttlSimSec = 10.0;
    m.turnRateRadS = 100.0;

    // Avoid accidental close-range collisions in the tests.
    m.radiusKm = 0.0;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = -1.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 1;

    // Track quality should instantly match measurement quality for test predictability.
    m.enableTrackQuality = true;
    m.trackQuality = 1.0;
    m.trackQualityRiseHalfLifeSimSec = 0.0;
    m.trackQualityFallHalfLifeSimSec = 0.0;
    m.trackQualityResistFloor = 0.25;

    // Enable the new range-based signal falloff model.
    m.seekerSignalHalfRangeKm = 200.0;
    m.seekerSignalMinDistKm = 1.0;

    // Ensure that "equal score" decoys only win when track quality is degraded.
    m.decoyResistance = 1.2;
    m.decoyCommitSimSec = 1.0;

    return m;
  };

  auto makeTarget = [](double xKm) {
    sim::SphereTarget t{};
    t.kind = sim::CombatTargetKind::Ship;
    t.index = 0;
    t.id = 1;
    t.centerKm = {xKm, 0, 0};
    t.velKmS = {0, 0, 0};
    t.radiusKm = 1.0;
    t.radarSignature = 1.0;
    return t;
  };

  auto makeDecoy = [](double xKm) {
    sim::SphereTarget d{};
    d.kind = sim::CombatTargetKind::Decoy;
    d.index = 1;
    d.id = 42;
    d.centerKm = {xKm, 0, 0};
    d.velKmS = {0, 0, 0};
    d.radiusKm = 1.0;
    d.decoyRadar = 1.0;
    return d;
  };

  // --- Case 1: long range -> weak received signal -> low track quality -> decoy should win. ---
  {
    sim::Missile m = makeMissile();

    // Far enough to avoid collisions and ensure signal falloff is severe.
    sim::SphereTarget tgt = makeTarget(1000.0);
    sim::SphereTarget decoy = makeDecoy(1000.0);

    std::vector<sim::Missile> missiles{m};
    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    sim::SphereTarget targets[2]{tgt, decoy};
    sim::stepMissiles(missiles, 0.05, targets, 2, dets, hits);

    if (missiles.empty()) {
      std::cerr << "FAIL: missile destroyed unexpectedly (long-range signal falloff)\n";
      return 1;
    }

    const double q = missiles[0].trackQuality;
    if (!(q < 0.10 && q > 0.01)) {
      std::cerr << "FAIL: expected low trackQuality at long range, got " << q << "\n";
      ++fails;
    }

    if (missiles[0].committedDecoyId != 42) {
      std::cerr << "FAIL: expected decoy commitment at long range; committedDecoyId="
                << missiles[0].committedDecoyId << "\n";
      ++fails;
    }
  }

  // --- Case 2: short range -> strong received signal -> high track quality -> decoy should NOT win. ---
  {
    sim::Missile m = makeMissile();

    sim::SphereTarget tgt = makeTarget(50.0);
    sim::SphereTarget decoy = makeDecoy(50.0);

    std::vector<sim::Missile> missiles{m};
    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    sim::SphereTarget targets[2]{tgt, decoy};
    sim::stepMissiles(missiles, 0.05, targets, 2, dets, hits);

    if (missiles.empty()) {
      std::cerr << "FAIL: missile destroyed unexpectedly (short-range signal falloff)\n";
      ++fails;
    } else {
      const double q = missiles[0].trackQuality;
      if (!(q > 0.80)) {
        std::cerr << "FAIL: expected higher trackQuality at shorter range, got " << q << "\n";
        ++fails;
      }

      if (missiles[0].committedDecoyId != 0) {
        std::cerr << "FAIL: expected no decoy commitment at shorter range; committedDecoyId="
                  << missiles[0].committedDecoyId << "\n";
        ++fails;
      }
    }
  }

  // Determinism check: repeat Case 1 and ensure the end state matches.
  {
    sim::Missile m = makeMissile();
    sim::SphereTarget tgt = makeTarget(1000.0);
    sim::SphereTarget decoy = makeDecoy(1000.0);

    std::vector<sim::Missile> a{m};
    std::vector<sim::Missile> b{m};
    std::vector<sim::MissileDetonation> detA;
    std::vector<sim::MissileHit> hitA;
    std::vector<sim::MissileDetonation> detB;
    std::vector<sim::MissileHit> hitB;

    sim::SphereTarget targets[2]{tgt, decoy};
    sim::stepMissiles(a, 0.05, targets, 2, detA, hitA);
    sim::stepMissiles(b, 0.05, targets, 2, detB, hitB);

    if (a.empty() || b.empty()) {
      std::cerr << "FAIL: determinism run missile destroyed unexpectedly\n";
      ++fails;
    } else {
      const double dq = std::fabs(a[0].trackQuality - b[0].trackQuality);
      if (dq > 1e-12) {
        std::cerr << "FAIL: trackQuality not deterministic; dq=" << dq << "\n";
        ++fails;
      }
      if (a[0].committedDecoyId != b[0].committedDecoyId) {
        std::cerr << "FAIL: committedDecoyId not deterministic\n";
        ++fails;
      }
    }
  }

  return fails;
}
