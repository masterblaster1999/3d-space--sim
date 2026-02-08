#include "stellar/sim/Combat.h"

#include "test_harness.h"

#include <cmath>
#include <iostream>
#include <vector>

using stellar::math::Vec3d;

static bool approx(double a, double b, double eps) {
  return std::fabs(a - b) <= eps;
}

int test_projectile_missile_interception() {
  using namespace stellar::sim;

  int failures = 0;

  // A missile and a projectile cross paths in 1 second. The projectile should
  // intercept and detonate the missile before it can impact the ship target.
  constexpr double dt = 1.0;

  std::vector<Projectile> projectiles;
  {
    Projectile p{};
    p.posKm = Vec3d{0, 0, 10};
    p.prevKm = p.posKm;
    p.velKmS = Vec3d{0, 0, -10};
    p.ttlSimSec = 2.0;
    p.radiusKm = 0.10;
    p.dmg = 1.0;
    p.fromPlayer = false;
    p.shooterId = 42;
    projectiles.push_back(p);
  }

  std::vector<Missile> missiles;
  {
    Missile m{};
    m.posKm = Vec3d{0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = Vec3d{0, 0, 10};
    m.ttlSimSec = 2.0;
    m.radiusKm = 0.15;
    m.blastRadiusKm = 1.0;
    m.dmg = 5.0;
    m.fromPlayer = true;
    m.shooterId = 7;

    // "Locked" target (ship at z=10). Without interception, this would be a direct hit at t=1.
    m.hasTarget = true;
    m.targetKind = CombatTargetKind::Ship;
    m.targetId = 999;

    missiles.push_back(m);
  }

  SphereTarget ship{};
  ship.kind = CombatTargetKind::Ship;
  ship.index = 0;
  ship.id = 999;
  ship.centerKm = Vec3d{0, 0, 10};
  ship.radiusKm = 0.5;
  ship.velKmS = Vec3d{0, 0, 0};

  const SphereTarget targets[] = {ship};

  // Step the projectile WITHOUT targets so it doesn't immediately collide with the ship.
  // (The interception logic uses the projectile's swept [prevKm->posKm] segment.)
  std::vector<ProjectileHit> projHits;
  stepProjectiles(projectiles, dt, /*targets=*/nullptr, /*targetCount=*/0, projHits);
  CHECK_EQ(projHits.size(), (std::size_t)0);

  std::vector<MissileDetonation> dets;
  std::vector<MissileHit> hits;
  stepMissilesWithProjectileInterception(missiles, projectiles, dt, targets, 1, dets, hits);

  CHECK_EQ(dets.size(), (std::size_t)1);
  if (dets.size() == 1) {
    const MissileDetonation& d = dets[0];

    CHECK(d.intercepted);
    CHECK_EQ(d.interceptorShooterId, (stellar::core::u64)42);
    CHECK_EQ(d.interceptorFromPlayer, false);

    // Relative segment is -10 -> +10, radius is (0.15 + 0.10)=0.25, so tEnter=(10-0.25)/20=0.4875.
    // Missile z = 10 * tEnter = 4.875.
    CHECK(approx(d.pointKm.z, 4.875, 1e-6));
  }

  // Explosion was far from the ship; the ship should take no splash damage.
  CHECK_EQ(hits.size(), (std::size_t)0);

  // Both the missile and projectile should have been consumed.
  CHECK_EQ(missiles.size(), (std::size_t)1);
  CHECK(missiles[0].ttlSimSec <= 0.0);

  CHECK_EQ(projectiles.size(), (std::size_t)1);
  CHECK(projectiles[0].ttlSimSec <= 0.0);

  return failures;
}
