#include "stellar/sim/Combat.h"
#include "stellar/sim/Countermeasures.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_combat() {
  int fails = 0;

  // --- Damage splits shield->hull ---
  {
    double shield = 10.0;
    double hull = 25.0;
    sim::applyDamage(6.0, shield, hull);
    if (!approx(shield, 4.0) || !approx(hull, 25.0)) {
      std::cerr << "[test_combat] damage should drain shield first. got shield="
                << shield << " hull=" << hull << "\n";
      ++fails;
    }
    sim::applyDamage(9.0, shield, hull);
    if (!approx(shield, 0.0) || !approx(hull, 20.0)) {
      std::cerr << "[test_combat] overflow should spill into hull. got shield="
                << shield << " hull=" << hull << "\n";
      ++fails;
    }
  }

  // --- Ray-sphere intersection (entry distance) ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 1};
    const math::Vec3d c{0, 0, 10};
    double t = 0.0;
    const bool hit = sim::raySphereIntersectKm(o, d, c, 1.0, t);
    if (!hit || !approx(t, 9.0, 1e-9)) {
      std::cerr << "[test_combat] raySphereIntersectKm wrong. hit=" << hit << " t=" << t << " expected t=9\n";
      ++fails;
    }
  }

  // --- Ray-sphere should reject spheres entirely behind the origin ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 1};
    const math::Vec3d c{0, 0, -10};
    double t = 123.0;
    const bool hit = sim::raySphereIntersectKm(o, d, c, 1.0, t);
    if (hit) {
      std::cerr << "[test_combat] raySphereIntersectKm should not hit spheres behind the ray. t=" << t << "\n";
      ++fails;
    }
  }

  // --- Ray-sphere should accept non-normalized direction vectors ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 2}; // not normalized
    const math::Vec3d c{0, 0, 10};
    double t = 0.0;
    const bool hit = sim::raySphereIntersectKm(o, d, c, 1.0, t);
    if (!hit || !approx(t, 9.0, 1e-9)) {
      std::cerr << "[test_combat] raySphereIntersectKm should treat dir as a direction. hit=" << hit << " t=" << t
                << " expected t=9\n";
      ++fails;
    }
  }

  // --- Ray-sphere starting inside should return t=0 ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 1};
    const math::Vec3d c{0, 0, 0};
    double t = 1.0;
    const bool hit = sim::raySphereIntersectKm(o, d, c, 1.0, t);
    if (!hit || !approx(t, 0.0, 1e-9)) {
      std::cerr << "[test_combat] raySphereIntersectKm should clamp entry when inside. hit=" << hit << " t=" << t
                << " expected t=0\n";
      ++fails;
    }
  }

  // --- Nearest target selection ---
  {
    sim::SphereTarget targets[2]{};
    targets[0].kind = sim::CombatTargetKind::Ship;
    targets[0].index = 0;
    targets[0].id = 111;
    targets[0].centerKm = {0, 0, 10};
    targets[0].radiusKm = 1.0;

    targets[1].kind = sim::CombatTargetKind::Asteroid;
    targets[1].index = 1;
    targets[1].id = 222;
    targets[1].centerKm = {0, 0, 6};
    targets[1].radiusKm = 1.0;

    const auto hit = sim::raycastNearestSphereKm({0, 0, 0}, {0, 0, 1}, 50.0, targets, 2);
    if (!hit.hit || hit.id != 222) {
      std::cerr << "[test_combat] raycastNearestSphereKm should hit nearer target first. got hit="
                << hit.hit << " id=" << hit.id << "\n";
      ++fails;
    }
  }

  // --- Projectile step + hit emission ---
  {
    std::vector<sim::Projectile> ps;
    sim::Projectile p{};
    p.posKm = {0, 0, 0};
    p.prevKm = p.posKm;
    p.velKmS = {0, 0, 10};
    p.ttlSimSec = 2.0;
    p.radiusKm = 0.1;
    p.dmg = 5.0;
    p.fromPlayer = true;
    p.shooterId = 0;
    ps.push_back(p);

    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 7;
    tgt.id = 999;
    tgt.centerKm = {0, 0, 5};
    tgt.radiusKm = 0.5;

    std::vector<sim::ProjectileHit> hits;
    sim::stepProjectiles(ps, 1.0, &tgt, 1, hits);

    if (hits.size() != 1 || hits[0].targetId != 999 || !approx(hits[0].dmg, 5.0)) {
      std::cerr << "[test_combat] stepProjectiles should emit one hit. hits=" << hits.size()
                << " id=" << (hits.empty() ? 0 : hits[0].targetId) << "\n";
      ++fails;
    }

    if (ps.empty() || ps[0].ttlSimSec > 0.0) {
      std::cerr << "[test_combat] projectile should be consumed on hit. ttl="
                << (ps.empty() ? -1.0 : ps[0].ttlSimSec) << "\n";
      ++fails;
    }
  }

  // --- Projectile step should account for moving targets (swept collision) ---
  {
    std::vector<sim::Projectile> ps;
    sim::Projectile p{};
    p.posKm = {0, 0, 0};
    p.prevKm = p.posKm;
    p.velKmS = {0, 0, 10};
    p.ttlSimSec = 2.0;
    p.radiusKm = 0.1;
    p.dmg = 5.0;
    p.fromPlayer = true;
    p.shooterId = 0;
    ps.push_back(p);

    // Target starts offset in X, but sweeps across the projectile path during dt=1.
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 3;
    tgt.id = 123;
    tgt.centerKm = {-1, 0, 5};
    tgt.velKmS = {2, 0, 0};
    tgt.radiusKm = 0.5;

    std::vector<sim::ProjectileHit> hits;
    sim::stepProjectiles(ps, 1.0, &tgt, 1, hits);

    // The relative motion should produce an entry hit before closest approach.
    if (hits.size() != 1 || hits[0].targetId != 123 || !approx(hits[0].pointKm.z, 4.4116515945854475, 1e-6)) {
      std::cerr << "[test_combat] swept projectile collision failed. hits=" << hits.size()
                << " id=" << (hits.empty() ? 0 : hits[0].targetId)
                << " z=" << (hits.empty() ? 0.0 : hits[0].pointKm.z) << "\n";
      ++fails;
    }
  }

  // --- Missile step + detonation + splash hit ---
  {
    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 2.0;
    m.radiusKm = 0.1;
    m.blastRadiusKm = 2.0;
    m.dmg = 5.0;
    m.fromPlayer = true;
    m.shooterId = 0;
    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 999;
    ms.push_back(m);

    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 7;
    tgt.id = 999;
    tgt.centerKm = {0, 0, 5};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 0.5;

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (dets.size() != 1) {
      std::cerr << "[test_combat] stepMissiles should emit one detonation. got=" << dets.size() << "\n";
      ++fails;
    }
    if (hits.size() != 1 || hits[0].targetId != 999 || !approx(hits[0].dmg, 5.0)) {
      std::cerr << "[test_combat] stepMissiles should splash-hit the target. hits=" << hits.size()
                << " id=" << (hits.empty() ? 0 : hits[0].targetId)
                << " dmg=" << (hits.empty() ? 0.0 : hits[0].dmg) << "\n";
      ++fails;
    }
    if (ms.empty() || ms[0].ttlSimSec > 0.0) {
      std::cerr << "[test_combat] missile should be consumed on detonation. ttl="
                << (ms.empty() ? -1.0 : ms[0].ttlSimSec) << "\n";
      ++fails;
    }
  }



  // --- Missile proximity fuse should detonate near closest approach (no direct collision) ---
  {
    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 2.0;
    m.radiusKm = 0.1;
    m.blastRadiusKm = 2.0;
    m.proximityFuseKm = 0.25; // enough to detonate on a near miss
    m.dmg = 5.0;
    m.fromPlayer = true;
    m.shooterId = 0;
    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 777;
    ms.push_back(m);

    // Offset so we miss the direct collision radius, but pass within fuse range.
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 1;
    tgt.id = 777;
    tgt.centerKm = {0.7, 0, 5};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 0.5;

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (dets.size() != 1) {
      std::cerr << "[test_combat] proximity fuse should detonate once. got=" << dets.size() << "\n";
      ++fails;
    } else {
      // Should detonate around z=5 (closest approach), not at entry.
      if (!approx(dets[0].pointKm.z, 5.0, 1e-6)) {
        std::cerr << "[test_combat] proximity detonation time unexpected. z=" << dets[0].pointKm.z << " expected ~5\n";
        ++fails;
      }
    }

    if (hits.size() != 1 || hits[0].targetId != 777 || !approx(hits[0].dmg, 4.5, 1e-6)) {
      std::cerr << "[test_combat] proximity fuse splash damage unexpected. hits=" << hits.size()
                << " id=" << (hits.empty() ? 0 : hits[0].targetId)
                << " dmg=" << (hits.empty() ? 0.0 : hits[0].dmg) << "\n";
      ++fails;
    }
    if (ms.empty() || ms[0].ttlSimSec > 0.0) {
      std::cerr << "[test_combat] proximity fuse missile should be consumed on detonation. ttl="
                << (ms.empty() ? -1.0 : ms[0].ttlSimSec) << "\n";
      ++fails;
    }
  }


  // --- Missile motor burn should accelerate and increase traveled distance ---
  {
    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 10.0;
    m.radiusKm = 0.1;

    m.thrustAccelKmS2 = 10.0;
    m.maxSpeedKmS = 20.0;
    m.motorBurnRemainingSimSec = 1.0;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, nullptr, 0, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] missile motor test: missile unexpectedly destroyed.\n";
      ++fails;
    } else {
      // With v0=10, a=10 for 1s: distance = 10*1 + 0.5*10*1^2 = 15, speed=20.
      if (!approx(ms[0].posKm.z, 15.0, 1e-9)) {
        std::cerr << "[test_combat] missile motor distance wrong. z=" << ms[0].posKm.z << " expected 15\n";
        ++fails;
      }
      if (!approx(ms[0].velKmS.length(), 20.0, 1e-9)) {
        std::cerr << "[test_combat] missile motor speed wrong. v=" << ms[0].velKmS.length() << " expected 20\n";
        ++fails;
      }
      if (!approx(ms[0].motorBurnRemainingSimSec, 0.0, 1e-9)) {
        std::cerr << "[test_combat] missile motor burn time not consumed. remaining=" << ms[0].motorBurnRemainingSimSec
                  << " expected 0\n";
        ++fails;
      }
    }
  }

  // --- Missile motor should clamp to maxSpeedKmS ---
  {
    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 10.0;
    m.radiusKm = 0.1;

    m.thrustAccelKmS2 = 10.0;
    m.maxSpeedKmS = 15.0;
    m.motorBurnRemainingSimSec = 2.0;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 2.0, nullptr, 0, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] missile motor clamp test: missile unexpectedly destroyed.\n";
      ++fails;
    } else {
      // v0=10, a=10, vmax=15. Reach vmax at t=0.5.
      // distance = 10*0.5 + 0.5*10*0.5^2 + 15*(2-0.5) = 28.75
      if (!approx(ms[0].posKm.z, 28.75, 1e-9)) {
        std::cerr << "[test_combat] missile motor clamp distance wrong. z=" << ms[0].posKm.z << " expected 28.75\n";
        ++fails;
      }
      if (!approx(ms[0].velKmS.length(), 15.0, 1e-9)) {
        std::cerr << "[test_combat] missile motor clamp speed wrong. v=" << ms[0].velKmS.length() << " expected 15\n";
        ++fails;
      }
      if (!approx(ms[0].motorBurnRemainingSimSec, 0.0, 1e-9)) {
        std::cerr << "[test_combat] missile motor clamp burn time not consumed. remaining=" << ms[0].motorBurnRemainingSimSec
                  << " expected 0\n";
        ++fails;
      }
    }
  }




  // --- G-limited turning should respect end-of-step speed when the missile accelerates ---
  {
    // Missile starts pointed +Z, target far along +X.
    // With maxLateralAccelKmS2=2 and the missile accelerating from 10->20 km/s over dt=1,
    // the effective turn rate should be limited by the *end* speed (20), yielding:
    //   turnRate = 2 / 20 = 0.1 rad/s  -> maxAngle = 0.1 rad over dt=1.
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 5;
    tgt.centerKm = {1000, 0, 0};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 0.0;
    m.blastRadiusKm = 0.0;

    m.turnRateRadS = 100.0; // high, so the G-limit dominates
    m.maxLateralAccelKmS2 = 2.0;

    m.thrustAccelKmS2 = 10.0;
    m.maxSpeedKmS = 20.0;
    m.motorBurnRemainingSimSec = 1.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 5;
    m.guidance = sim::MissileGuidance::LeadPursuit;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] G-limit test: missile unexpectedly destroyed.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      if (!approx(d.x, std::sin(0.1), 0.01) || !approx(d.z, std::cos(0.1), 0.01)) {
        std::cerr << "[test_combat] G-limit turn angle unexpected. dir=(" << d.x << "," << d.y << "," << d.z
                  << ") expected approx (" << std::sin(0.1) << ",0," << std::cos(0.1) << ")\n";
        ++fails;
      }
    }
  }

  // --- Countermeasures integrate + expire ---
  {
    std::vector<sim::Countermeasure> cms;
    sim::Countermeasure c{};
    c.id = 1;
    c.posKm = {0, 0, 0};
    c.velKmS = {0, 0, 1};
    c.radiusKm = 0.1;
    c.ttlSimSec = 2.0;
    c.ttlMaxSimSec = 2.0;
    c.heatStrength = 10.0;
    c.radarStrength = 0.0;
    cms.push_back(c);

    sim::stepCountermeasures(cms, 1.0);
    if (cms.size() != 1 || !approx(cms[0].posKm.z, 1.0) || !approx(cms[0].ttlSimSec, 1.0)) {
      std::cerr << "[test_combat] countermeasure integrate failed. size=" << cms.size()
                << " z=" << (cms.empty() ? -1.0 : cms[0].posKm.z)
                << " ttl=" << (cms.empty() ? -1.0 : cms[0].ttlSimSec) << "\n";
      ++fails;
    }

    sim::stepCountermeasures(cms, 1.5);
    if (!cms.empty()) {
      std::cerr << "[test_combat] countermeasure should expire once ttl<=0. remaining=" << cms.size() << "\n";
      ++fails;
    }
  }

  // --- Missile seekers should prefer strong decoy targets (flares) ---
  {
    // Locked target (ship) is far.
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 1;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;

    // A nearby flare-style decoy.
    std::vector<sim::Countermeasure> cms;
    sim::Countermeasure flare{};
    flare.id = 100;
    flare.type = sim::CountermeasureType::Flare;
    flare.posKm = {10, 0, 10};
    flare.velKmS = {0, 0, 0};
    flare.radiusKm = 0.5;
    flare.ttlSimSec = 10.0;
    flare.ttlMaxSimSec = 10.0;
    flare.heatStrength = 50.0; // strong
    flare.radarStrength = 0.0;
    cms.push_back(flare);

    std::vector<sim::SphereTarget> targets;
    targets.push_back(shipT);
    sim::appendCountermeasureTargets(cms, targets);

    // Missile starts pointing at the ship (positive Z), but should steer toward the decoy.
    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 50.0; // allow aggressive steering

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 1;

    m.seeker = sim::MissileSeekerType::Heat;
    m.seekerFovCos = 0.0;
    m.decoyResistance = 1.0;

    const auto toDecoy = (flare.posKm - m.posKm).normalized();

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 0.1, targets.data(), targets.size(), dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] missile unexpectedly destroyed while testing decoy steering.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
      if (dot < 0.95) {
        std::cerr << "[test_combat] missile did not steer toward decoy. dot=" << dot
                  << " vel=(" << d.x << "," << d.y << "," << d.z << ")" << "\n";
        ++fails;
      }
    }
  }

  // --- Heat seekers should weight target thermal signatures when comparing against decoys ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 1;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;

    // A flare-like decoy at a similar range but offset in X so guidance is distinguishable.
    sim::SphereTarget decoyT{};
    decoyT.kind = sim::CombatTargetKind::Decoy;
    decoyT.index = 1;
    decoyT.id = 100;
    decoyT.centerKm = {50, 0, 200};
    decoyT.velKmS = {0, 0, 0};
    decoyT.radiusKm = 0.5;
    decoyT.decoyHeat = 15.0;
    decoyT.decoyRadar = 0.0;

    sim::SphereTarget targets[2]{shipT, decoyT};

    auto makeMissile = [&]() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 50.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 1;

      m.seeker = sim::MissileSeekerType::Heat;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = 0.0;
      m.decoyResistance = 1.0;
      return m;
    };

    // Case A: a hot ship should not be overridden by a moderately hot decoy.
    targets[0].heatSignature = 20.0;
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] heatSignature test (A): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
        const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
        if (dot > 0.985) {
          std::cerr << "[test_combat] heatSignature test (A): missile should favor hot target over decoy. dot="
                    << dot << "\n";
          ++fails;
        }
      }
    }

    // Case B: a cooler ship should be overridden by the same decoy.
    targets[0].heatSignature = 10.0;
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] heatSignature test (B): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
        const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
        if (dot < 0.99) {
          std::cerr << "[test_combat] heatSignature test (B): expected decoy override for cooler target. dot="
                    << dot << "\n";
          ++fails;
        }
      }
    }
  }


  // --- Radar signature weighting (big contacts resist chaff more than small) ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 2;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;
    shipT.radarSignature = 50.0;

    // A chaff-like decoy at a similar range but offset in X so guidance is distinguishable.
    sim::SphereTarget decoyT{};
    decoyT.kind = sim::CombatTargetKind::Decoy;
    decoyT.index = 1;
    decoyT.id = 200;
    decoyT.centerKm = {50, 0, 200};
    decoyT.velKmS = {0, 0, 0};
    decoyT.radiusKm = 1.0;
    decoyT.decoyHeat = 0.0;
    decoyT.decoyRadar = 15.0;

    sim::SphereTarget targets[2]{shipT, decoyT};

    auto makeMissile = [&]() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 100.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 2;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = 0.0;
      m.decoyResistance = 1.0;
      m.decoyCommitSimSec = 0.0;
      return m;
    };

    // Case A: strong radar signature should resist the decoy.
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] radarSignature test (A): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
        const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
        if (dot > 0.985) {
          std::cerr << "[test_combat] radarSignature test (A): expected strong target to beat decoy. dot="
                    << dot << "\n";
          ++fails;
        }
      }
    }

    // Case B: weak radar signature should be overridden by the same decoy.
    targets[0].radarSignature = 10.0;
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] radarSignature test (B): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
        const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
        if (dot < 0.99) {
          std::cerr << "[test_combat] radarSignature test (B): expected decoy override for weak target. dot="
                    << dot << "\n";
          ++fails;
        }
      }
    }
  }

  
  // --- Heat seeker aspect weighting (tail-on hotter than head-on) ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 1;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;
    shipT.heatSignature = 15.0;

    // A flare-like decoy at a similar range but offset in X so guidance is distinguishable.
    sim::SphereTarget decoyT{};
    decoyT.kind = sim::CombatTargetKind::Decoy;
    decoyT.index = 1;
    decoyT.id = 100;
    decoyT.centerKm = {50, 0, 200};
    decoyT.velKmS = {0, 0, 0};
    decoyT.radiusKm = 1.0;
    decoyT.decoyHeat = 15.0;
    decoyT.decoyRadar = 0.0;

    sim::SphereTarget targets[2]{shipT, decoyT};

    auto makeMissile = [&]() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 50.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 1;

      m.seeker = sim::MissileSeekerType::Heat;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = 0.0;
      m.decoyResistance = 1.0;

      // Force strong aspect sensitivity for the test.
      m.heatAspectFrontFactor = 0.55;
      m.heatAspectRearFactor = 1.45;
      m.heatAspectMinSpeedKmS = 0.0;
      m.heatAspectSpeedForFullKmS = 0.01;

      return m;
    };

    // Case A: tail chase (target moving away) should resist the decoy.
    targets[0].velKmS = {0, 0, 0.25};
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] heatAspect test (A): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
        const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
        if (dot > 0.985) {
          std::cerr << "[test_combat] heatAspect test (A): expected tail-on ship to beat decoy. dot="
                    << dot << "\n";
          ++fails;
        }
      }
    }

    // Case B: head-on (target closing) should be easier to decoy.
    targets[0].velKmS = {0, 0, -0.25};
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] heatAspect test (B): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
        const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
        if (dot < 0.99) {
          std::cerr << "[test_combat] heatAspect test (B): expected decoy override head-on. dot="
                    << dot << "\n";
          ++fails;
        }
      }
    }
  }

// --- HOJ should allow radar seekers to ignore doppler notch against jamming targets ---
  {
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 9001;
    // Place the target to the +X side so the missile must turn.
    tgt.centerKm = {10, 0, 0};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;
    // Active jammer emission.
    tgt.jammerPower = 1.0;

    auto makeMissile = [&](bool hoj) {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 20.0; // enough to swing to +X in one step

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 9001;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0; // always inside cone
      m.requireLineOfSight = false;

      // Always notched in this geometry (radial velocity is ~0 along +X).
      m.radarDopplerNotchKmS = 5.0;

      m.guidance = sim::MissileGuidance::LeadPursuit;

      m.homeOnJam = hoj;
      m.homeOnJamMinJammerPower = 0.25;
      m.homeOnJamTrackQualityCap = 0.70;
      return m;
    };

    // Without HOJ: lock should break due to notch and the missile should not turn.
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile(false));

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, &tgt, 1, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] HOJ test (disabled): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        if (std::fabs(d.x) > 0.02 || d.z < 0.98) {
          std::cerr << "[test_combat] HOJ test (disabled): missile turned despite notch. dir=(" << d.x << "," << d.y
                    << "," << d.z << ")\n";
          ++fails;
        }
      }
    }

    // With HOJ: the jammer should allow the missile to keep tracking and turn toward +X.
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile(true));

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, &tgt, 1, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] HOJ test (enabled): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        if (d.x < 0.25) {
          std::cerr << "[test_combat] HOJ test (enabled): missile did not turn toward jammer. dir=(" << d.x << "," << d.y
                    << "," << d.z << ")\n";
          ++fails;
        }
      }
    }
  }


  // --- HOJ received-jamming gating: far target should not bypass notch unless jamming is strong enough ---
  {
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 9100;
    tgt.centerKm = {300, 0, 0};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;
    tgt.jammerPower = 1.0;

    auto makeMissile = [&]() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 20.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 9100;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;

      m.radarDopplerNotchKmS = 5.0;

      m.guidance = sim::MissileGuidance::LeadPursuit;

      m.homeOnJam = true;
      m.homeOnJamMinJammerPower = 0.25;
      // Require a minimum received jamming level, with distance-aware scaling.
      m.radarJammingHalfRangeKm = 100.0;
      m.radarJammingMinDistKm = 1.0;
      m.homeOnJamMinJamming01 = 0.20;

      return m;
    };

    // Case A: far range -> received jamming is weak -> HOJ should NOT bypass the notch.
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, &tgt, 1, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] HOJ jamming01 test (far): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        if (std::fabs(d.x) > 0.02 || d.z < 0.98) {
          std::cerr << "[test_combat] HOJ jamming01 test (far): missile turned despite weak received jamming. dir=(" << d.x
                    << "," << d.y << "," << d.z << ")\n";
          ++fails;
        }
      }
    }

    // Case B: closer range -> received jamming is stronger -> HOJ should bypass the notch and turn.
    tgt.centerKm = {100, 0, 0};
    {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, &tgt, 1, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] HOJ jamming01 test (near): missile unexpectedly destroyed.\n";
        ++fails;
      } else {
        const auto d = ms[0].velKmS.normalized();
        if (d.x < 0.25) {
          std::cerr << "[test_combat] HOJ jamming01 test (near): missile did not turn with strong received jamming. dir=(" << d.x
                    << "," << d.y << "," << d.z << ")\n";
          ++fails;
        }
      }
    }
  }

  // --- ProNav steering should match expected sign/magnitude in a simple case ---
  {
    // Geometry picked so the closed-form PN turn command is easy to reason about.
    // Missile starts pointed +Z; target is offset +X,+Z and stationary.

    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 42;
    tgt.centerKm = {10, 0, 10};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 100.0; // keep unclamped for predictable turn

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 42;

    m.guidance = sim::MissileGuidance::ProNav;
    m.navConstant = 3.0;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 0.1, &tgt, 1, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] ProNav missile unexpectedly destroyed while testing steering.\n";
      ++fails;
    } else {
      // Expected: a modest turn toward +X, ~sin(0.106) ≈ 0.105.
      const auto d = ms[0].velKmS.normalized();
      if (!approx(d.x, 0.105, 0.02) || !approx(d.y, 0.0, 0.02) || !approx(d.z, 0.994, 0.02)) {
        std::cerr << "[test_combat] ProNav steering unexpected. dir=(" << d.x << "," << d.y << "," << d.z
                  << ") expected approx (0.105,0,0.994)\n";
        ++fails;
      }
    }
  }

  // --- APN (augmented ProNav) should respond to target lateral acceleration ---
  {
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 77;
    tgt.centerKm = {0, 0, 20};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;

    auto makeBaseMissile = []() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 100.0; // unclamped

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 77;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0; // always inside cone
      m.decoyResistance = 1.0e9;

      m.guidance = sim::MissileGuidance::ProNav;
      m.navConstant = 3.0;
      return m;
    };

    std::vector<sim::Missile> pn;
    std::vector<sim::Missile> apn;

    pn.push_back(makeBaseMissile());
    {
      auto m = makeBaseMissile();
      m.apnTargetAccelGain = 0.5 * m.navConstant;
      m.apnMaxTargetAccelKmS2 = 1.0e6; // effectively no clamp for this test
      m.apnAccelHalfLifeSimSec = 0.0;  // raw estimate
      apn.push_back(m);
    }

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1: establish baseline velocity sample (no lateral motion).
    sim::stepMissiles(pn, 0.1, &tgt, 1, dets, hits);
    dets.clear(); hits.clear();
    sim::stepMissiles(apn, 0.1, &tgt, 1, dets, hits);
    dets.clear(); hits.clear();

    // Step 2: target picks up a lateral velocity (finite-difference acceleration).
    tgt.velKmS = {1, 0, 0};
    sim::stepMissiles(pn, 0.1, &tgt, 1, dets, hits);
    dets.clear(); hits.clear();
    sim::stepMissiles(apn, 0.1, &tgt, 1, dets, hits);

    if (pn.empty() || apn.empty()) {
      std::cerr << "[test_combat] APN test: missile unexpectedly destroyed.\n";
      ++fails;
    } else {
      const double xPN = pn[0].velKmS.normalized().x;
      const double xAPN = apn[0].velKmS.normalized().x;
      if (!(xAPN > xPN + 0.01)) {
        std::cerr << "[test_combat] APN did not increase lateral turn as expected. xPN=" << xPN
                  << " xAPN=" << xAPN << "\n";
        ++fails;
      }
    }
  }


  // --- Lead pursuit should account for target acceleration when provided ---
  {
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 1234;
    tgt.centerKm = {0, 0, 100};
    tgt.velKmS = {0, 0, 0};
    tgt.accelKmS2 = {1.0, 0, 0}; // constant lateral accel (km/s^2)
    tgt.radiusKm = 1.0;

    auto makeMissile = []() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 60.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;

      // High turn authority so the first step reflects the computed lead direction.
      m.turnRateRadS = 1000.0;
      m.maxTurnAccelRadS2 = 0.0; // disable lag for deterministic snapping in this test

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 1234;

      // Ensure the seeker does not block tracking.
      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;
      m.requireLineOfSight = false;

      m.guidance = sim::MissileGuidance::LeadPursuit;
      return m;
    };

    auto runX = [&](const sim::SphereTarget& t) {
      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, &t, 1, dets, hits);

      if (ms.empty()) return 0.0;
      return ms[0].velKmS.normalized().x;
    };

    // Case A: with zero target acceleration, lead pursuit should aim essentially straight.
    sim::SphereTarget noA = tgt;
    noA.accelKmS2 = {0, 0, 0};
    const double xNoA = runX(noA);

    // Case B: with strong lateral acceleration, the lead point should shift in +X.
    const double xA = runX(tgt);

    if (std::fabs(xNoA) > 0.05) {
      std::cerr << "[test_combat] accel lead test: expected near-zero x without accel. xNoA=" << xNoA << "\n";
      ++fails;
    }
    if (!(xA > 0.55)) {
      std::cerr << "[test_combat] accel lead test: expected strong +X lead with accel. xA=" << xA << "\n";
      ++fails;
    }
  }


  // --- Decoy commitment should hold guidance for a minimum time window ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 1;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;

    sim::SphereTarget decoyT{};
    decoyT.kind = sim::CombatTargetKind::Decoy;
    decoyT.index = 1;
    decoyT.id = 1000;
    decoyT.centerKm = {10, 0, 10};
    decoyT.velKmS = {0, 0, 0};
    decoyT.radiusKm = 0.5;
    decoyT.decoyHeat = 50.0;
    decoyT.decoyRadar = 0.0;

    sim::SphereTarget targets[2]{shipT, decoyT};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 100.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 1;

    m.seeker = sim::MissileSeekerType::Heat;
    m.seekerFovCos = 0.0;
    m.decoyResistance = 1.0;
    m.decoyCommitSimSec = 0.5;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1: strong decoy should override the true target.
    sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] decoy commit test: missile unexpectedly destroyed in step 1.\n";
      ++fails;
    } else {
      const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
      const auto d = ms[0].velKmS.normalized();
      const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
      if (dot < 0.90) {
        std::cerr << "[test_combat] decoy commit test: step 1 did not steer to decoy. dot=" << dot << "\n";
        ++fails;
      }
    }

    // Step 2: decoy becomes weak, so scoring would normally return to the ship.
    // Commitment should keep steering toward the same decoy for a short time.
    targets[1].decoyHeat = 1.0e-6;
    sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] decoy commit test: missile unexpectedly destroyed in step 2.\n";
      ++fails;
    } else {
      const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
      const auto d = ms[0].velKmS.normalized();
      const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
      if (dot < 0.85) {
        std::cerr << "[test_combat] decoy commit test: step 2 should remain committed to decoy. dot=" << dot << "\n";
        ++fails;
      }
    }

    // Step 3: after the commitment window expires, the missile should be free to switch back.
    sim::stepMissiles(ms, 0.5, targets, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] decoy commit test: missile unexpectedly destroyed in step 3.\n";
      ++fails;
    } else {
      const auto toDecoy = (targets[1].centerKm - ms[0].posKm).normalized();
      const auto toShip = (targets[0].centerKm - ms[0].posKm).normalized();
      const auto d = ms[0].velKmS.normalized();

      const double dotDecoy = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
      const double dotShip = d.x * toShip.x + d.y * toShip.y + d.z * toShip.z;

      if (dotShip <= dotDecoy) {
        std::cerr << "[test_combat] decoy commit test: expected switch back toward ship after hold. "
                  << "dotShip=" << dotShip << " dotDecoy=" << dotDecoy << "\n";
        ++fails;
      }
    }
  }

  // --- Track quality: lock breaks should make decoys more effective ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 500;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;

    // A radar decoy aligned with the memory track. We set its strength so that it is
    // *not* strong enough to override a good lock, but *is* strong enough to override
    // once track quality decays.
    sim::SphereTarget decoyT{};
    decoyT.kind = sim::CombatTargetKind::Decoy;
    decoyT.index = 1;
    decoyT.id = 9000;
    decoyT.centerKm = {0, 0, 200};
    decoyT.velKmS = {0, 0, 0};
    decoyT.radiusKm = 0.5;
    decoyT.decoyHeat = 0.0;
    decoyT.decoyRadar = 0.5;

    sim::SphereTarget targets[2]{shipT, decoyT};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 100.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 500;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = std::cos(0.10); // narrow-ish cone
    m.decoyResistance = 1.0;
    m.decoyCommitSimSec = 1.0;
    m.targetMemorySimSec = 5.0;

    // Enable track quality with a fast decay so the test runs quickly.
    m.enableTrackQuality = true;
    m.trackQuality = 1.0;
    m.trackQualityRiseHalfLifeSimSec = 0.10;
    m.trackQualityFallHalfLifeSimSec = 0.20;
    m.trackQualityResistFloor = 0.25;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1: with a good lock, the half-strength decoy should NOT override.
    sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);
    if (ms.empty()) {
      std::cerr << "[test_combat] track quality test: missile destroyed unexpectedly in step 1.\n";
      ++fails;
    } else if (ms[0].committedDecoyId != 0) {
      std::cerr << "[test_combat] track quality test: decoy should not commit with good track. committedId="
                << ms[0].committedDecoyId << "\n";
      ++fails;
    }

    // Break the lock by moving the target outside the seeker's FOV. Target memory should
    // keep a base track, but track quality will decay.
    targets[0].centerKm = {200, 0, 0};

    // Step 2: after a lock break, reduced track quality should allow the same decoy to override.
    sim::stepMissiles(ms, 1.0, targets, 2, dets, hits);
    if (ms.empty()) {
      std::cerr << "[test_combat] track quality test: missile destroyed unexpectedly in step 2.\n";
      ++fails;
    } else if (ms[0].committedDecoyId != 9000) {
      std::cerr << "[test_combat] track quality test: expected decoy commit after lock break. committedId="
                << ms[0].committedDecoyId << "\n";
      ++fails;
    }
  }



  // --- Track quality: noise jamming can degrade radar measurement quality (optional) ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 510;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;
    shipT.jammerPower = 1.0;

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    auto makeMissile = [&]() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 1.0e9; // ensure heading snaps to the guide direction

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 510;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;
      m.decoyResistance = 1.0e9;

      m.enableTrackQuality = true;
      m.trackQuality = 1.0;
      m.trackQualityRiseHalfLifeSimSec = 0.10;
      m.trackQualityFallHalfLifeSimSec = 0.25;
      m.trackQualityResistFloor = 1.0;

      return m;
    };

    // Control: no suppression gain -> track quality remains ~1 while directly tracking.
    {
      std::vector<sim::Missile> ms;
      auto m = makeMissile();
      m.radarJammingTrackSuppressionGain = 0.0;
      ms.push_back(m);

      sim::stepMissiles(ms, 0.5, &shipT, 1, dets, hits);
      if (ms.empty()) {
        std::cerr << "[test_combat] jamming track quality test (control): missile destroyed unexpectedly.\n";
        ++fails;
      } else if (ms[0].trackQuality < 0.95) {
        std::cerr << "[test_combat] jamming track quality test (control): expected ~no degradation. q=" << ms[0].trackQuality
                  << "\n";
        ++fails;
      }
    }

    // With suppression gain + distance scaling: track quality should degrade below 1.
    {
      std::vector<sim::Missile> ms;
      auto m = makeMissile();
      m.radarJammingHalfRangeKm = 200.0; // at 200 km -> jamming01 ~0.5
      m.radarJammingMinDistKm = 1.0;
      m.radarJammingTrackSuppressionGain = 2.0;
      ms.push_back(m);

      sim::stepMissiles(ms, 0.5, &shipT, 1, dets, hits);
      if (ms.empty()) {
        std::cerr << "[test_combat] jamming track quality test (suppressed): missile destroyed unexpectedly.\n";
        ++fails;
      } else if (ms[0].trackQuality > 0.90) {
        std::cerr << "[test_combat] jamming track quality test (suppressed): expected degradation under jamming. q="
                  << ms[0].trackQuality << "\n";
        ++fails;
      }
    }
  }


  // --- Seeker slew: limited gimbal rate should lag behind fast LOS changes ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 700;
    shipT.centerKm = {0, 0, 100};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 100.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 700;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = 0.0; // 90-degree cone
    m.seekerSlewRateRadS = 0.5; // deliberately slow (rad/s)
    m.decoyResistance = 1e9;

    m.enableTrackQuality = true;
    m.trackQuality = 1.0;
    m.trackQualityRiseHalfLifeSimSec = 0.20;
    m.trackQualityFallHalfLifeSimSec = 2.0;
    m.trackQualityResistFloor = 1.0;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 0: seed the seeker direction while the target is straight ahead.
    sim::stepMissiles(ms, 0.1, &shipT, 1, dets, hits);
    if (ms.empty()) {
      std::cerr << "[test_combat] seeker slew test: missile destroyed unexpectedly in step 0.\n";
      ++fails;
    }

    // Abrupt LOS change: target jumps to +X, at the edge of the 90-degree FOV.
    shipT.centerKm = {100, 0, 0};

    sim::stepMissiles(ms, 1.0, &shipT, 1, dets, hits);
    if (ms.empty()) {
      std::cerr << "[test_combat] seeker slew test: missile destroyed unexpectedly in step 1.\n";
      ++fails;
    } else {
      // With a 0.5 rad/s gimbal and dt=1, the seeker should slew by ~0.5 rad.
      const sim::Missile& mm = ms[0];
      const math::Vec3d s = mm.seekerDirWorld.lengthSq() < 1e-12 ? math::Vec3d{0, 0, 1} : mm.seekerDirWorld.normalized();
      const double ex = std::sin(0.5);
      const double ez = std::cos(0.5);
      if (!approx(s.x, ex, 0.06) || !approx(s.z, ez, 0.06)) {
        std::cerr << "[test_combat] seeker slew test: unexpected seeker dir after slew. got=(" << s.x << "," << s.y << "," << s.z
                  << ") expected~(" << ex << ",0," << ez << ")\n";
        ++fails;
      }

      // Track quality should decay because the seeker was initially pointed away from the new LOS.
      const double expectedQ = std::exp(-math::kLn2 * 1.0 / 2.0); // half-life=2s over 1s
      if (!approx(mm.trackQuality, expectedQ, 0.06)) {
        std::cerr << "[test_combat] seeker slew test: unexpected track quality after LOS jump. got=" << mm.trackQuality
                  << " expected~" << expectedQ << "\n";
        ++fails;
      }
    }
  }


  // --- Turn acceleration limiting: missile should not snap instantly to large heading changes ---
  {
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 701;
    tgt.centerKm = {100, 0, 0};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.turnRateRadS = 10.0; // high enough that turn rate is not the limiting factor

    // Enable the lag model: limit how fast turn rate can build.
    m.maxTurnAccelRadS2 = 1.0; // rad/s^2

    // Always track in this test.
    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = -1.0;
    m.decoyResistance = 1e9;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 701;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1: desired heading change is ~90deg (1.57 rad) but the missile can only ramp
    // turn rate to 1 rad/s over this 1s frame -> it rotates by ~1 rad.
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] turn accel test: missile destroyed unexpectedly in step 1.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      const double ex = std::sin(1.0);
      const double ez = std::cos(1.0);
      if (!approx(d.x, ex, 0.06) || !approx(d.z, ez, 0.06)) {
        std::cerr << "[test_combat] turn accel test: unexpected heading after accel-limited turn. dir=(" << d.x
                  << "," << d.y << "," << d.z << ") expected~(" << ex << ",0," << ez << ")\n";
        ++fails;
      }
    }

    // Step 2: should be able to finish the remaining turn and nearly point at the target.
    dets.clear();
    hits.clear();
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] turn accel test: missile destroyed unexpectedly in step 2.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      if (d.x < 0.98 || std::fabs(d.z) > 0.10) {
        std::cerr << "[test_combat] turn accel test: expected near-full turn by step 2. dir=(" << d.x << "," << d.y
                  << "," << d.z << ")\n";
        ++fails;
      }
    }
  }



  // --- Seeker FOV: missile should not steer toward targets outside its seeker cone ---
  {
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 7;
    tgt.centerKm = {0, 0, -100}; // behind the missile
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 20.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 100.0; // if it could see the target, it would turn quickly

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 7;

    // 90 deg FOV; target at 180 deg should be untrackable.
    m.seekerFovCos = 0.0;
    m.targetMemorySimSec = 0.0;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 0.1, &tgt, 1, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] seeker FOV test: missile unexpectedly destroyed.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      if (std::fabs(d.x) > 1e-3 || std::fabs(d.y) > 1e-3 || d.z < 0.999) {
        std::cerr << "[test_combat] seeker FOV test: expected no turn toward rear target. dir=("
                  << d.x << "," << d.y << "," << d.z << ")\n";
        ++fails;
      }
    }
  }

  // --- Target memory: missile should keep turning for a short window after losing target updates ---
  {
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 9;
    tgt.centerKm = {10, 0, 10};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 1.0;

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 0.5; // 0.05 rad per 0.1s step

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 9;

    m.seekerFovCos = 0.0;
    m.targetMemorySimSec = 0.15;

    m.guidance = sim::MissileGuidance::LeadPursuit;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1: target present, establishes last-known state.
    sim::stepMissiles(ms, 0.1, &tgt, 1, dets, hits);
    if (ms.empty()) {
      std::cerr << "[test_combat] target memory test: missile destroyed in step 1.\n";
      ++fails;
    } else {
      const auto d1 = ms[0].velKmS.normalized();

      // Step 2: target missing, memory should keep turning.
      sim::stepMissiles(ms, 0.1, nullptr, 0, dets, hits);
      if (ms.empty()) {
        std::cerr << "[test_combat] target memory test: missile destroyed in step 2.\n";
        ++fails;
      } else {
        const auto d2 = ms[0].velKmS.normalized();
        if (d2.x <= d1.x + 0.02) {
          std::cerr << "[test_combat] target memory test: expected continued turn under memory. "
                    << "d1.x=" << d1.x << " d2.x=" << d2.x << "\n";
          ++fails;
        }

        // Step 3: memory expires, so direction should stop changing (ballistic).
        sim::stepMissiles(ms, 0.1, nullptr, 0, dets, hits);
        if (ms.empty()) {
          std::cerr << "[test_combat] target memory test: missile destroyed in step 3.\n";
          ++fails;
        } else {
          const auto d3 = ms[0].velKmS.normalized();
          if (!approx(d3.x, d2.x, 1e-6) || !approx(d3.y, d2.y, 1e-6) || !approx(d3.z, d2.z, 1e-6)) {
            std::cerr << "[test_combat] target memory test: expected no further turn after memory expiry. "
                      << "d2=(" << d2.x << "," << d2.y << "," << d2.z << ") "
                      << "d3=(" << d3.x << "," << d3.y << "," << d3.z << ")\n";
            ++fails;
          }
        }
      }
    }
  }


  // --- Seeker activation delay: decoys should not affect guidance until the seeker is active ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 10;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;

    sim::SphereTarget decoyT{};
    decoyT.kind = sim::CombatTargetKind::Decoy;
    decoyT.index = 1;
    decoyT.id = 2000;
    decoyT.centerKm = {10, 0, 10};
    decoyT.velKmS = {0, 0, 0};
    decoyT.radiusKm = 0.5;
    decoyT.decoyHeat = 0.0;
    decoyT.decoyRadar = 50.0;

    sim::SphereTarget targets[2]{shipT, decoyT};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 1.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 100.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 10;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerFovCos = 0.0;
    m.decoyResistance = 1.0;

    // Seeker is inactive for the first step, becomes active on the second step.
    m.seekerActivationSimSec = 0.1;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1 (inactive seeker): should not steer toward the decoy.
    sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);
    if (ms.empty()) {
      std::cerr << "[test_combat] seeker activation test: missile destroyed in step 1.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      if (std::fabs(d.x) > 0.02) {
        std::cerr << "[test_combat] seeker activation test: expected no early decoy turn. dir.x=" << d.x << "\n";
        ++fails;
      }
    }

    // Step 2 (active seeker): decoy should override and attract the missile.
    sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);
    if (ms.empty()) {
      std::cerr << "[test_combat] seeker activation test: missile destroyed in step 2.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      const auto toDecoy = (decoyT.centerKm - ms[0].posKm).normalized();
      const double dot = d.x * toDecoy.x + d.y * toDecoy.y + d.z * toDecoy.z;
      if (dot < 0.85) {
        std::cerr << "[test_combat] seeker activation test: expected decoy attraction after activation. dot=" << dot
                  << "\n";
        ++fails;
      }
    }
  }




  // --- Seeker update period: decoy/reacquire decisions should only update on scan ticks ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 10;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;
    shipT.heatSignature = 1.0;

    sim::SphereTarget decoyT{};
    decoyT.kind = sim::CombatTargetKind::Decoy;
    decoyT.index = 1;
    decoyT.id = 30000;
    decoyT.centerKm = {50, 0, 200};
    decoyT.velKmS = {0, 0, 0};
    decoyT.radiusKm = 0.5;
    decoyT.decoyHeat = 250.0; // very strong flare-like signal

    sim::SphereTarget targetsNoDecoy[1]{shipT};
    sim::SphereTarget targetsWithDecoy[2]{shipT, decoyT};

    std::vector<sim::Missile> ms;

    auto makeMissile = [&](double updatePeriodSec) {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.proximityFuseKm = 0.0;
      m.turnRateRadS = 100.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 10;

      m.seeker = sim::MissileSeekerType::Heat;
      m.seekerFovCos = -1.0; // wide cone
      m.seekerActivationSimSec = 0.0;

      m.decoyResistance = 1.0;
      m.decoyCommitSimSec = 5.0;
      m.targetMemorySimSec = 5.0;

      // Feature under test:
      //  - 0.0 => continuous (legacy)
      //  - >0  => scan ticks
      m.seekerUpdatePeriodSimSec = updatePeriodSec;
      return m;
    };

    // Missile A: continuous updates (legacy).
    ms.push_back(makeMissile(0.0));
    // Missile B: 1 Hz seeker update (rate-limited).
    ms.push_back(makeMissile(1.0));

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1: establish baseline memory on the ship target (no decoy present yet).
    sim::stepMissiles(ms, 0.1, targetsNoDecoy, 1, dets, hits);
    if (ms.size() < 2) {
      std::cerr << "[test_combat] seeker update period test: missiles destroyed unexpectedly in step 1.\n";
      ++fails;
    }

    // Step 2: the decoy appears.
    sim::stepMissiles(ms, 0.1, targetsWithDecoy, 2, dets, hits);
    if (ms.size() < 2) {
      std::cerr << "[test_combat] seeker update period test: missiles destroyed unexpectedly in step 2.\n";
      ++fails;
    } else {
      // Continuous-update missile should commit immediately.
      if (ms[0].committedDecoyId != decoyT.id) {
        std::cerr << "[test_combat] seeker update period test: expected immediate decoy commit for continuous seeker. committedId="
                  << ms[0].committedDecoyId << "\n";
        ++fails;
      }
      // Rate-limited seeker should NOT commit until its next scan tick.
      if (ms[1].committedDecoyId != 0) {
        std::cerr << "[test_combat] seeker update period test: expected no early commit for rate-limited seeker. committedId="
                  << ms[1].committedDecoyId << "\n";
        ++fails;
      }
    }

    // Advance until the rate-limited missile gets a scan tick (period=1.0s).
    for (int i = 0; i < 30 && ms.size() >= 2; ++i) {
      sim::stepMissiles(ms, 0.1, targetsWithDecoy, 2, dets, hits);
      if (ms[1].committedDecoyId == decoyT.id) break;
    }

    if (ms.size() < 2) {
      std::cerr << "[test_combat] seeker update period test: missiles destroyed unexpectedly while advancing time.\n";
      ++fails;
    } else if (ms[1].committedDecoyId != decoyT.id) {
      std::cerr << "[test_combat] seeker update period test: expected decoy commit on later scan tick. committedId="
                << ms[1].committedDecoyId << "\n";
      ++fails;
    }
  }

  // --- Midcourse datalink gating: optional LOS/range should control pre-activation target updates ---
  {
    // Shooter platform (launching ship).
    sim::SphereTarget shooter{};
    shooter.kind = sim::CombatTargetKind::Ship;
    shooter.index = 0;
    shooter.id = 99;
    shooter.centerKm = {0, 0, 0};
    shooter.velKmS = {0, 0, 0};
    shooter.radiusKm = 1.0;

    // Locked target.
    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 1;
    tgt.id = 1;
    tgt.centerKm = {1000, 0, 0};
    tgt.velKmS = {0, 1, 0};
    tgt.radiusKm = 1.0;

    sim::SphereTarget targets1[2]{shooter, tgt};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.dmg = 0.0;
    m.blastRadiusKm = 0.0;
    m.turnRateRadS = 10.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 1;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerFovCos = 0.0;

    // Keep the seeker inactive so only the midcourse update path is exercised.
    m.seekerActivationSimSec = 10.0;

    // Enable datalink gating.
    m.datalinkRangeKm = 2000.0;
    m.datalinkRequireLineOfSight = true;
    m.datalinkOcclusionPadKm = 0.0;

    // Keep inertial memory alive long enough to observe a missed update.
    m.targetMemorySimSec = 5.0;

    m.fromPlayer = true;
    m.shooterId = 99;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    // Step 1: datalink has LOS, so the missile should ingest the target velocity sample.
    sim::stepMissiles(ms, 0.1, targets1, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] datalink gating test: missile destroyed in step 1.\n";
      ++fails;
    } else if (!ms[0].hasLastKnownTarget) {
      std::cerr << "[test_combat] datalink gating test: expected lastKnownTarget to be seeded.\n";
      ++fails;
    } else if (!approx(ms[0].lastKnownTargetVelKmS.y, 1.0, 1e-9)) {
      std::cerr << "[test_combat] datalink gating test: expected lastKnownTargetVel.y=1 after update. got "
                << ms[0].lastKnownTargetVelKmS.y << "\n";
      ++fails;
    }

    // Step 2: introduce an occluding asteroid and change the target velocity.
    sim::SphereTarget asteroid{};
    asteroid.kind = sim::CombatTargetKind::Asteroid;
    asteroid.index = 2;
    asteroid.id = 500;
    asteroid.centerKm = {500, 0, 0};
    asteroid.velKmS = {0, 0, 0};
    asteroid.radiusKm = 200.0;

    sim::SphereTarget tgt2 = tgt;
    tgt2.velKmS = {0, -10, 0};

    sim::SphereTarget targets2[3]{shooter, tgt2, asteroid};

    sim::stepMissiles(ms, 0.1, targets2, 3, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] datalink gating test: missile destroyed in step 2.\n";
      ++fails;
    } else {
      // Datalink should be blocked by LOS, so the last-known velocity should remain the old sample.
      if (!approx(ms[0].lastKnownTargetVelKmS.y, 1.0, 1e-9)) {
        std::cerr << "[test_combat] datalink gating test: expected lastKnownTargetVel unchanged under LOS break. got "
                  << ms[0].lastKnownTargetVelKmS.y << "\n";
        ++fails;
      }
      if (ms[0].lastKnownTargetAgeSimSec <= 0.0) {
        std::cerr << "[test_combat] datalink gating test: expected lastKnownTargetAge to advance under LOS break. age="
                  << ms[0].lastKnownTargetAgeSimSec << "\n";
        ++fails;
      }
    }
  }



  // --- Radar decoy discrimination gates: reject implausible decoys by range/angle/doppler ---
  {
    sim::SphereTarget shipT{};
    shipT.kind = sim::CombatTargetKind::Ship;
    shipT.index = 0;
    shipT.id = 55;
    shipT.centerKm = {0, 0, 200};
    shipT.velKmS = {0, 0, 0};
    shipT.radiusKm = 1.0;

    // Decoy A: very strong but far in range/angle relative to the ship track.
    sim::SphereTarget decoyA{};
    decoyA.kind = sim::CombatTargetKind::Decoy;
    decoyA.index = 1;
    decoyA.id = 3000;
    decoyA.centerKm = {10, 0, 10};
    decoyA.velKmS = {0, 0, 0};
    decoyA.radiusKm = 0.5;
    decoyA.decoyHeat = 0.0;
    decoyA.decoyRadar = 50.0;

    // Decoy B: strong and close to the ship track (should pass gates and win once A is rejected).
    sim::SphereTarget decoyB{};
    decoyB.kind = sim::CombatTargetKind::Decoy;
    decoyB.index = 2;
    decoyB.id = 3001;
    decoyB.centerKm = {50, 0, 190};
    decoyB.velKmS = {0, 0, 0};
    decoyB.radiusKm = 0.5;
    decoyB.decoyHeat = 0.0;
    decoyB.decoyRadar = 50.0;

    sim::SphereTarget targets[3]{shipT, decoyA, decoyB};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.turnRateRadS = 100.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 55;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = 0.0;
    m.decoyResistance = 1.0;

    // Tight-ish gates so A is rejected but B passes.
    m.decoyAngleGateCos = std::cos(0.35); // ~20 deg
    m.decoyRangeGateKm = 20.0;
    m.decoyDopplerGateKmS = 1.0;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 0.1, targets, 3, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] decoy gate test: missile destroyed unexpectedly.\n";
      ++fails;
    } else {
      const auto d = ms[0].velKmS.normalized();
      // If decoyA won, d.x would be ~0.707. If no decoy won, d.x ~0.0.
      // With gates, we expect a moderate +X turn toward decoyB (d.x ~0.255).
      if (d.x < 0.20 || d.x > 0.35) {
        std::cerr << "[test_combat] decoy gate test: expected steering toward decoyB (moderate +X). dir=("
                  << d.x << "," << d.y << "," << d.z << ")\n";
        ++fails;
      }
    }
  }


  // --- Radar burn-through: close-range decoy attraction should be reduced when enabled ---
  {
    auto makeMissile = []() {
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.dmg = 0.0;
      m.blastRadiusKm = 0.0;
      m.turnRateRadS = 100.0; // allow full snap for deterministic direction

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 1;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0; // always inside cone
      m.decoyResistance = 1.0;

      // Enable burn-through so decoys become less effective at close range.
      m.decoyBurnThroughRangeKm = 100.0;
      m.decoyBurnThroughMinFactor = 0.05;
      return m;
    };

    auto runCase = [&](double rangeKm) {
      sim::SphereTarget shipT{};
      shipT.kind = sim::CombatTargetKind::Ship;
      shipT.index = 0;
      shipT.id = 1;
      shipT.centerKm = {0, 0, rangeKm};
      shipT.velKmS = {0, 0, 0};
      shipT.radiusKm = 1.0;

      // Place a decoy at the SAME range but offset in angle.
      const double x = 0.2 * rangeKm;
      const double z = std::sqrt(std::max(0.0, rangeKm * rangeKm - x * x));

      sim::SphereTarget decoyT{};
      decoyT.kind = sim::CombatTargetKind::Decoy;
      decoyT.index = 1;
      decoyT.id = 2;
      decoyT.centerKm = {x, 0, z};
      decoyT.velKmS = {0, 0, 0};
      decoyT.radiusKm = 0.5;
      decoyT.decoyHeat = 0.0;
      decoyT.decoyRadar = 1.1; // just strong enough to override without burn-through

      sim::SphereTarget targets[2]{shipT, decoyT};

      std::vector<sim::Missile> ms;
      ms.push_back(makeMissile());

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        return 0.0;
      }
      return ms[0].velKmS.normalized().x;
    };

    const double xFar = runCase(200.0);
    const double xClose = runCase(20.0);

    if (!(xFar > 0.10)) {
      std::cerr << "[test_combat] burn-through test: expected far-range decoy attraction. xFar=" << xFar << "\n";
      ++fails;
    }
    if (std::fabs(xClose) > 0.05) {
      std::cerr << "[test_combat] burn-through test: expected close-range decoy rejection. xClose=" << xClose
                << "\n";
      ++fails;
    }
  }


  // --- Seeker slew: large-angle decoy pulls should be suppressed when the seeker can't re-point quickly ---
  {
    auto runCase = [&](double slewRateRadS) {
      sim::SphereTarget shipT{};
      shipT.kind = sim::CombatTargetKind::Ship;
      shipT.index = 0;
      shipT.id = 1;
      shipT.centerKm = {0, 0, 100};
      shipT.velKmS = {0, 0, 0};
      shipT.radiusKm = 1.0;

      // Decoy placed at the same range but 30 degrees off-boresight.
      const double rangeKm = 100.0;
      const double angRad = 0.5235987755982988; // 30 deg
      const double x = rangeKm * std::sin(angRad);
      const double z = rangeKm * std::cos(angRad);

      sim::SphereTarget decoyT{};
      decoyT.kind = sim::CombatTargetKind::Decoy;
      decoyT.index = 1;
      decoyT.id = 2;
      decoyT.centerKm = {x, 0, z};
      decoyT.velKmS = {0, 0, 0};
      decoyT.radiusKm = 0.5;
      decoyT.decoyHeat = 0.0;
      decoyT.decoyRadar = 1.2; // slightly stronger than the true target (when slew is ignored)

      sim::SphereTarget targets[2]{shipT, decoyT};

      std::vector<sim::Missile> ms;
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.turnRateRadS = 100.0; // snap for deterministic steering

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 1;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;
      m.decoyResistance = 1.0;

      // This is the key: limit how fast the seeker can swing to a new direction.
      m.seekerSlewRateRadS = slewRateRadS;
      m.enableTrackQuality = false;

      ms.push_back(m);

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) return 0.0;
      return ms[0].velKmS.normalized().x;
    };

    const double xFast = runCase(0.0);
    const double xSlow = runCase(0.5);

    if (!(xFast > 0.20)) {
      std::cerr << "[test_combat] seeker slew decoy test: expected decoy pull when slew is unlimited. xFast=" << xFast
                << "\n";
      ++fails;
    }
    if (std::fabs(xSlow) > 0.05) {
      std::cerr << "[test_combat] seeker slew decoy test: expected lock to resist large-angle decoy when slew is slow. xSlow="
                << xSlow << "\n";
      ++fails;
    }
  }


  // --- Auto-acquire: missile without initial lock should acquire a target within range/FOV when enabled ---
  {
    sim::SphereTarget shipA{};
    shipA.kind = sim::CombatTargetKind::Ship;
    shipA.index = 0;
    shipA.id = 101;
    shipA.centerKm = {0, 0, 50};
    shipA.velKmS = {0, 0, 0};
    shipA.radiusKm = 1.0;

    sim::SphereTarget shipB{};
    shipB.kind = sim::CombatTargetKind::Ship;
    shipB.index = 1;
    shipB.id = 102;
    shipB.centerKm = {50, 0, 50};
    shipB.velKmS = {0, 0, 0};
    shipB.radiusKm = 1.0;

    sim::SphereTarget targets[2]{shipA, shipB};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.turnRateRadS = 100.0;

    // No initial target lock.
    m.hasTarget = false;
    m.targetId = 0;

    // Enable reacquisition.
    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.autoAcquireRangeKm = 100.0;
    m.seekerFovCos = 0.0;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] auto-acquire test: missile destroyed unexpectedly.\n";
      ++fails;
    } else {
      if (!ms[0].hasTarget || ms[0].targetId != 101) {
        std::cerr << "[test_combat] auto-acquire test: expected lock on shipA (id 101). got hasTarget="
                  << ms[0].hasTarget << " id=" << ms[0].targetId << "\n";
        ++fails;
      }
    }
  }


  // --- Line-of-sight occlusion: asteroids should prevent seeker tracking when enabled ---
  {
    sim::SphereTarget ship{};
    ship.kind = sim::CombatTargetKind::Ship;
    ship.index = 0;
    ship.id = 500;
    ship.centerKm = {0, 0, 100};
    ship.velKmS = {0, 0, 0};
    ship.radiusKm = 1.0;

    sim::SphereTarget rock{};
    rock.kind = sim::CombatTargetKind::Asteroid;
    rock.index = 1;
    rock.id = 501;
    rock.centerKm = {0, 0, 50};
    rock.velKmS = {0, 0, 0};
    rock.radiusKm = 10.0;

    sim::SphereTarget targets[2]{ship, rock};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {1, 0, 10}; // slightly off-axis
    m.ttlSimSec = 30.0;
    m.radiusKm = 0.1;
    m.turnRateRadS = 100.0;

    // Disable memory + auto-acquire so LOS breaks guidance immediately.
    m.seeker = sim::MissileSeekerType::Heat;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = -1.0;
    m.requireLineOfSight = true;
    m.lineOfSightOcclusionPadKm = 0.0;
    m.targetMemorySimSec = 0.0;
    m.autoAcquireRangeKm = 0.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 500;

    ms.push_back(m);

    const double x0 = ms[0].velKmS.normalized().x;

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, targets, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_combat] LOS occlusion test: missile destroyed unexpectedly.\n";
      ++fails;
    } else {
      const double x1 = ms[0].velKmS.normalized().x;
      if (x1 < 0.05) {
        std::cerr << "[test_combat] LOS occlusion test: expected no guidance turn while occluded. x0="
                  << x0 << " x1=" << x1 << "\n";
        ++fails;
      }

      // Remove the asteroid: guidance should resume and turn toward the ship (negative X).
      sim::SphereTarget clearTargets[1]{ship};
      dets.clear();
      hits.clear();
      sim::stepMissiles(ms, 1.0, clearTargets, 1, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] LOS occlusion test: missile destroyed unexpectedly after clearing LOS.\n";
        ++fails;
      } else {
        const double x2 = ms[0].velKmS.normalized().x;
        if (x2 > 0.0) {
          std::cerr << "[test_combat] LOS occlusion test: expected guidance turn toward ship after clearing LOS. x2="
                    << x2 << "\n";
          ++fails;
        }
      }
    }
  }


  // --- Asteroid avoidance: enabled missiles should bias away from imminent asteroid collisions ---
  {
    sim::SphereTarget ship{};
    ship.kind = sim::CombatTargetKind::Ship;
    ship.index = 0;
    ship.id = 700;
    ship.centerKm = {0, 0, 200};
    ship.velKmS = {0, 0, 0};
    ship.radiusKm = 1.0;

    sim::SphereTarget rock{};
    rock.kind = sim::CombatTargetKind::Asteroid;
    rock.index = 1;
    rock.id = 701;
    rock.centerKm = {0, 0, 50};
    rock.velKmS = {0, 0, 0};
    rock.radiusKm = 10.0;

    sim::SphereTarget targets[2]{ship, rock};

    // Baseline: without avoidance the missile should keep flying straight.
    {
      std::vector<sim::Missile> ms;
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 20};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.turnRateRadS = 2.0;

      m.seeker = sim::MissileSeekerType::Heat;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;
      m.requireLineOfSight = false;
      m.targetMemorySimSec = 0.0;
      m.autoAcquireRangeKm = 0.0;

      m.asteroidAvoidanceStrength = 0.0;
      m.asteroidAvoidanceLookaheadSimSec = 2.0;
      m.asteroidAvoidancePadKm = 0.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 700;

      ms.push_back(m);

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] asteroid avoidance test: baseline missile destroyed unexpectedly.\n";
        ++fails;
      } else {
        const math::Vec3d d = ms[0].velKmS.normalized();
        if (std::abs(d.x) > 1e-6 || std::abs(d.y) > 1e-6) {
          std::cerr << "[test_combat] asteroid avoidance test: expected straight flight when disabled. dir="
                    << d.x << "," << d.y << "," << d.z << "\n";
          ++fails;
        }
      }
    }

    // With avoidance enabled, the missile should introduce a lateral bias away from the rock.
    {
      std::vector<sim::Missile> ms;
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 20};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.turnRateRadS = 2.0;

      m.seeker = sim::MissileSeekerType::Heat;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;
      m.requireLineOfSight = false;
      m.targetMemorySimSec = 0.0;
      m.autoAcquireRangeKm = 0.0;

      m.asteroidAvoidanceStrength = 0.85;
      m.asteroidAvoidanceLookaheadSimSec = 2.0;
      m.asteroidAvoidancePadKm = 0.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 700;

      ms.push_back(m);

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 0.1, targets, 2, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] asteroid avoidance test: missile destroyed unexpectedly.\n";
        ++fails;
      } else {
        const math::Vec3d d = ms[0].velKmS.normalized();
        if (std::abs(d.x) < 0.05 && std::abs(d.y) < 0.05) {
          std::cerr << "[test_combat] asteroid avoidance test: expected lateral bias when enabled. dir="
                    << d.x << "," << d.y << "," << d.z << "\n";
          ++fails;
        }
      }
    }
  }


  // --- Radar doppler notch: near-zero radial velocity should break radar tracking when enabled ---
  {
    sim::SphereTarget ship{};
    ship.kind = sim::CombatTargetKind::Ship;
    ship.index = 0;
    ship.id = 600;
    ship.centerKm = {10, 0, 100};
    ship.velKmS = {0, 0, 10}; // match missile velocity -> vr ~ 0
    ship.radiusKm = 1.0;

    sim::SphereTarget targets[1]{ship};

    // Notch disabled: missile should steer toward the offset target.
    {
      std::vector<sim::Missile> ms;
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.turnRateRadS = 100.0;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;
      m.radarDopplerNotchKmS = 0.0;
      m.targetMemorySimSec = 0.0;
      m.autoAcquireRangeKm = 0.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 600;

      ms.push_back(m);

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 1.0, targets, 1, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] radar notch test: missile destroyed unexpectedly (notch disabled).\n";
        ++fails;
      } else {
        const double x = ms[0].velKmS.normalized().x;
        if (x < 0.05) {
          std::cerr << "[test_combat] radar notch test: expected steering when notch disabled. x=" << x
                    << "\n";
          ++fails;
        }
      }
    }

    // Notch enabled: vr ~ 0 should prevent tracking, so missile should not steer.
    {
      std::vector<sim::Missile> ms;
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 30.0;
      m.radiusKm = 0.1;
      m.turnRateRadS = 100.0;

      m.seeker = sim::MissileSeekerType::Radar;
      m.seekerActivationSimSec = 0.0;
      m.seekerFovCos = -1.0;
      m.radarDopplerNotchKmS = 0.5;
      m.targetMemorySimSec = 0.0;
      m.autoAcquireRangeKm = 0.0;

      m.hasTarget = true;
      m.targetKind = sim::CombatTargetKind::Ship;
      m.targetId = 600;

      ms.push_back(m);

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 1.0, targets, 1, dets, hits);

      if (ms.empty()) {
        std::cerr << "[test_combat] radar notch test: missile destroyed unexpectedly (notch enabled).\n";
        ++fails;
      } else {
        const double x = ms[0].velKmS.normalized().x;
        if (std::fabs(x) > 0.02) {
          std::cerr << "[test_combat] radar notch test: expected no steering when notched. x=" << x
                    << "\n";
          ++fails;
        }
      }
    }
  }




  // --- Directional blast: rear targets should receive reduced splash damage when enabled ---
  {
    sim::SphereTarget front{};
    front.kind = sim::CombatTargetKind::Ship;
    front.index = 0;
    front.id = 800;
    front.centerKm = {0, 0, 5};
    front.velKmS = {0, 0, 0};
    front.radiusKm = 0.0;

    sim::SphereTarget rear{};
    rear.kind = sim::CombatTargetKind::Ship;
    rear.index = 1;
    rear.id = 801;
    rear.centerKm = {0, 0, -5};
    rear.velKmS = {0, 0, 0};
    rear.radiusKm = 0.0;

    sim::SphereTarget detonator{};
    detonator.kind = sim::CombatTargetKind::Asteroid;
    detonator.index = 2;
    detonator.id = 802;
    detonator.centerKm = {0, 0, 0};
    detonator.velKmS = {0, 0, 0};
    detonator.radiusKm = 0.0;

    sim::SphereTarget targets[3]{front, rear, detonator};

    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.prevKm = {0, 0, 0};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 5.0;
    m.radiusKm = 0.1;
    m.dmg = 100.0;
    m.blastRadiusKm = 10.0;
    m.proximityFuseKm = 0.0;
    m.turnRateRadS = 0.0;

    m.blastDirectionalStrength = 1.0;
    m.blastDirectionalMinFactor = 0.2;

    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, targets, 3, dets, hits);

    double dmgFront = -1.0;
    double dmgRear = -1.0;
    for (const auto& h : hits) {
      if (h.targetId == 800) dmgFront = h.dmg;
      if (h.targetId == 801) dmgRear = h.dmg;
    }

    if (dmgFront <= 0.0 || dmgRear <= 0.0) {
      std::cerr << "[test_combat] directional blast test: expected splash hits on both targets. front="
                << dmgFront << " rear=" << dmgRear << "\n";
      ++fails;
    } else {
      const double ratio = dmgRear / dmgFront;
      if (ratio > 0.30) {
        std::cerr << "[test_combat] directional blast test: expected rear damage to be significantly lower. ratio="
                  << ratio << " front=" << dmgFront << " rear=" << dmgRear << "\n";
        ++fails;
      }
    }
  }


  // --- Blast LOS occlusion: asteroids should block (or attenuate) splash damage when enabled ---
  {
    sim::SphereTarget ship{};
    ship.kind = sim::CombatTargetKind::Ship;
    ship.index = 0;
    ship.id = 810;
    ship.centerKm = {0, 0, 10};
    ship.velKmS = {0, 0, 0};
    ship.radiusKm = 0.0;

    sim::SphereTarget occluder{};
    occluder.kind = sim::CombatTargetKind::Asteroid;
    occluder.index = 1;
    occluder.id = 811;
    occluder.centerKm = {0, 0, 5};
    occluder.velKmS = {0, 0, 0};
    occluder.radiusKm = 2.0;

    sim::SphereTarget detonator{};
    detonator.kind = sim::CombatTargetKind::Asteroid;
    detonator.index = 2;
    detonator.id = 812;
    detonator.centerKm = {0, 0, 0};
    detonator.velKmS = {0, 0, 0};
    detonator.radiusKm = 0.0;

    sim::SphereTarget targets[3]{ship, occluder, detonator};

    // Occlusion enabled: expect no hit (occluded factor = 0).
    {
      std::vector<sim::Missile> ms;
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 5.0;
      m.radiusKm = 0.1;
      m.dmg = 100.0;
      m.blastRadiusKm = 20.0;
      m.turnRateRadS = 0.0;

      m.blastRequireLineOfSight = true;
      m.blastOccludedFactor = 0.0;

      ms.push_back(m);

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 1.0, targets, 3, dets, hits);

      bool hitShip = false;
      for (const auto& h : hits) {
        if (h.targetId == 810) hitShip = true;
      }
      if (hitShip) {
        std::cerr << "[test_combat] blast LOS occlusion test: expected ship to be occluded (no hit).\n";
        ++fails;
      }
    }

    // Occlusion disabled: expect a hit.
    {
      std::vector<sim::Missile> ms;
      sim::Missile m{};
      m.prevKm = {0, 0, 0};
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 10};
      m.ttlSimSec = 5.0;
      m.radiusKm = 0.1;
      m.dmg = 100.0;
      m.blastRadiusKm = 20.0;
      m.turnRateRadS = 0.0;

      m.blastRequireLineOfSight = false;

      ms.push_back(m);

      std::vector<sim::MissileDetonation> dets;
      std::vector<sim::MissileHit> hits;
      sim::stepMissiles(ms, 1.0, targets, 3, dets, hits);

      bool hitShip = false;
      for (const auto& h : hits) {
        if (h.targetId == 810) hitShip = true;
      }
      if (!hitShip) {
        std::cerr << "[test_combat] blast LOS occlusion test: expected ship hit when occlusion disabled.\n";
        ++fails;
      }
    }
  }

  // --- Missile swarm coordination: separation bias should deterministically split paths ---
  {
    std::vector<sim::Missile> ms;

    sim::Missile m0{};
    m0.prevKm = {-0.02, 0, 0};
    m0.posKm = m0.prevKm;
    m0.velKmS = {0, 0, 10};
    m0.ttlSimSec = 10.0;
    m0.radiusKm = 0.05;
    m0.dmg = 0.0;
    m0.blastRadiusKm = 0.0;
    m0.turnRateRadS = 1.2;
    m0.maxLateralAccelKmS2 = 0.0;
    m0.maxTurnAccelRadS2 = 0.0;
    m0.fromPlayer = true;
    m0.shooterId = 42;
    m0.hasTarget = true;
    m0.targetKind = sim::CombatTargetKind::Ship;
    m0.targetId = 7001;
    m0.seekerFovCos = -1.0;
    m0.seekerActivationSimSec = 0.0;

    m0.swarmSeparationStrength = 1.0;
    m0.swarmCohesionStrength = 0.0;
    m0.swarmAlignmentStrength = 0.0;
    m0.swarmSeparationKm = 0.20;
    m0.swarmNeighborRangeKm = 1.0;
    m0.swarmMaxSteerRad = 0.0; // use default fraction of turn capability

    sim::Missile m1 = m0;
    m1.prevKm = {0.02, 0, 0};
    m1.posKm = m1.prevKm;

    ms.push_back(m0);
    ms.push_back(m1);

    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 7001;
    tgt.centerKm = {0, 0, 100};
    tgt.radiusKm = 0.5;

    const double sep0 = (ms[0].posKm - ms[1].posKm).length();

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 0.5, &tgt, 1, dets, hits);

    const double sep1 = (ms[0].posKm - ms[1].posKm).length();
    if (!(sep1 > sep0 + 0.25)) {
      std::cerr << "[test_combat] swarm separation should increase spacing. sep0=" << sep0 << " sep1=" << sep1
                << "\n";
      ++fails;
    }

    if (!(ms[0].posKm.x < -0.02 && ms[1].posKm.x > 0.02)) {
      std::cerr << "[test_combat] swarm separation should split laterally. x0=" << ms[0].posKm.x
                << " x1=" << ms[1].posKm.x << "\n";
      ++fails;
    }
  }

  // --- Ballistic projectile aim assist biases shots toward a lead solution ---
  {
    // Arrange: shooter faces +Z, target is slightly offset in X so a straight shot misses
    // (given combined radii), but a small aim correction should score a hit.
    sim::Ship shooter0{};
    shooter0.setPositionKm({0, 0, 0});
    shooter0.setVelocityKmS({0, 0, 0});

    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 0;
    tgt.id = 77;
    tgt.centerKm = {900.0, 0.0, 100000.0};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 50.0;
    tgt.minAimCos = std::cos(math::degToRad(5.0));

    // Fire without providing targets -> no aim assist -> should miss.
    {
      sim::Ship shooter = shooter0;
      const sim::FireResult fr = sim::tryFireWeapon(shooter, sim::WeaponType::Cannon, /*cooldownSimSec*/0.0, /*distributorMk*/1, /*shooterId*/123, /*fromPlayer*/true, nullptr, 0);
      if (!fr.hasProjectile) {
        std::cerr << "[test_combat] expected projectile fire result\n";
        ++fails;
      } else {
        std::vector<sim::Projectile> ps;
        ps.push_back(fr.projectile);
        std::vector<sim::ProjectileHit> hits;
        sim::stepProjectiles(ps, /*dtSimSec*/900.0, &tgt, 1, hits);
        if (!hits.empty()) {
          std::cerr << "[test_combat] non-assisted projectile should miss (hits=" << hits.size() << ")\n";
          ++fails;
        }
      }
    }

    // Fire with provided targets -> aim assist can engage -> should hit.
    {
      sim::Ship shooter = shooter0;
      const sim::FireResult fr = sim::tryFireWeapon(shooter, sim::WeaponType::Cannon, /*cooldownSimSec*/0.0, /*distributorMk*/1, /*shooterId*/123, /*fromPlayer*/true, &tgt, 1);
      if (!fr.hasProjectile) {
        std::cerr << "[test_combat] expected projectile fire result (assisted)\n";
        ++fails;
      } else {
        std::vector<sim::Projectile> ps;
        ps.push_back(fr.projectile);
        std::vector<sim::ProjectileHit> hits;
        sim::stepProjectiles(ps, /*dtSimSec*/900.0, &tgt, 1, hits);
        if (hits.size() != 1 || hits[0].targetId != tgt.id) {
          std::cerr << "[test_combat] assisted projectile should hit (hits=" << hits.size() << ")\n";
          ++fails;
        }
      }
    }

    // Fire with provided targets but aim-assist weight=0 -> should behave like unassisted and miss.
    {
      sim::Ship shooter = shooter0;
      sim::SphereTarget weak = tgt;
      weak.aimAssistWeight01 = 0.0;
      const sim::FireResult fr = sim::tryFireWeapon(shooter, sim::WeaponType::Cannon, /*cooldownSimSec*/0.0, /*distributorMk*/1, /*shooterId*/123, /*fromPlayer*/true, &weak, 1);
      if (!fr.hasProjectile) {
        std::cerr << "[test_combat] expected projectile fire result (weighted)\n";
        ++fails;
      } else {
        std::vector<sim::Projectile> ps;
        ps.push_back(fr.projectile);
        std::vector<sim::ProjectileHit> hits;
        sim::stepProjectiles(ps, /*dtSimSec*/900.0, &weak, 1, hits);
        if (!hits.empty()) {
          std::cerr << "[test_combat] weight=0 projectile should miss (hits=" << hits.size() << ")\n";
          ++fails;
        }
      }
    }
  }
  if (fails == 0) {
    std::cout << "[test_combat] PASS\n";
  }
  return fails;
}
