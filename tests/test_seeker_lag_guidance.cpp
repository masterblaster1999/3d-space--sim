#include "stellar/sim/Combat.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_seeker_lag_guidance() {
  int fails = 0;

  // Scenario:
  // - Missile is initially flying +Z.
  // - Locked target is due +X (90 deg off-boresight).
  // - Missile has a very high turn rate, but the seeker gimbal slew is very low.
  // - With seekerLagGuidanceGain disabled, missile snaps to the true lead direction.
  // - With seekerLagGuidanceGain enabled, the commanded guidance direction is
  //   biased toward the lagged seeker LOS, leaving a noticeable +Z component.

  sim::SphereTarget tgt{};
  tgt.kind = sim::CombatTargetKind::Ship;
  tgt.index = 0;
  tgt.id = 42;
  tgt.centerKm = {100, 0, 0};
  tgt.velKmS = {0, 0, 0};
  tgt.radiusKm = 1.0;

  const math::Vec3d toDir = (tgt.centerKm - math::Vec3d{0, 0, 0}).normalized();

  auto makeMissile = [&](double lagGain) {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 5.0;
    m.radiusKm = 0.1;
    m.blastRadiusKm = 0.0;
    m.dmg = 0.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 42;

    // Wide FOV so the target is always inside the seeker cone.
    m.seekerFovCos = -1.0;

    // Seeker is active immediately and updates continuously.
    m.seekerActivationSimSec = 0.0;
    m.seekerUpdatePeriodSimSec = 0.0;

    // Missile can turn essentially instantly.
    m.turnRateRadS = 100.0;

    // Seeker gimbal is slow.
    m.seekerSlewRateRadS = 0.1;

    // Force an initial seeker pointing direction that is aligned to missile nose,
    // not the target (simulating a mid-flight misalignment / high-LOS-rate case).
    m.hasSeekerDir = true;
    m.seekerDirWorld = {0, 0, 1};

    m.seekerLagGuidanceGain = lagGain;
    return m;
  };

  // --- Baseline: no lag influence ---
  {
    std::vector<sim::Missile> ms;
    ms.push_back(makeMissile(/*lagGain=*/0.0));

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (!dets.empty() || !hits.empty()) {
      std::cerr << "[test_seeker_lag_guidance] baseline should not detonate/hit. dets="
                << dets.size() << " hits=" << hits.size() << "\n";
      ++fails;
    }

    const math::Vec3d dir = math::safeNormalized(ms[0].velKmS, math::Vec3d{0, 0, 1}, 1e-12);
    const double dotTo = math::dot(dir, toDir);

    if (!(dotTo > 0.999)) {
      std::cerr << "[test_seeker_lag_guidance] baseline should snap to target dir. dot="
                << dotTo << "\n";
      ++fails;
    }
  }

  // --- Enabled: lag should bias guidance toward seeker LOS ---
  {
    std::vector<sim::Missile> ms;
    ms.push_back(makeMissile(/*lagGain=*/1.0));

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (!dets.empty() || !hits.empty()) {
      std::cerr << "[test_seeker_lag_guidance] lag case should not detonate/hit. dets="
                << dets.size() << " hits=" << hits.size() << "\n";
      ++fails;
    }

    const math::Vec3d dir = math::safeNormalized(ms[0].velKmS, math::Vec3d{0, 0, 1}, 1e-12);
    const double dotTo = math::dot(dir, toDir);

    // With a slow seeker and large initial misalignment, guidance should NOT fully
    // rotate to +X in a single step.
    if (!(dotTo < 0.90)) {
      std::cerr << "[test_seeker_lag_guidance] lag should reduce alignment to target. dot="
                << dotTo << "\n";
      ++fails;
    }

    // Direction should retain a meaningful forward (+Z) component.
    if (!(dir.z > 0.30)) {
      std::cerr << "[test_seeker_lag_guidance] lag should retain +Z component. z="
                << dir.z << "\n";
      ++fails;
    }

    // Sanity: ensure normalization stays intact.
    if (!approx(dir.length(), 1.0, 1e-9)) {
      std::cerr << "[test_seeker_lag_guidance] direction should be unit length. |dir|="
                << dir.length() << "\n";
      ++fails;
    }
  }

  return fails;
}
