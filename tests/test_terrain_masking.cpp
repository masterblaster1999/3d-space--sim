#include "stellar/sim/MissileDefense.h"

#include "stellar/math/Geometry.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

int test_terrain_masking() {
  int fails = 0;

  // Basic: pick a nearby asteroid and recommend a direction toward a cover point.
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 0};

    const math::Vec3d tgtPos{10, 0, 0};
    const math::Vec3d tgtVel{0, 0, 0};

    sim::SphereTarget a{};
    a.kind = sim::CombatTargetKind::Asteroid;
    a.id = 111;
    a.centerKm = {5, 3, 0};
    a.radiusKm = 2.0;

    const sim::SphereTarget targets[] = {a};

    sim::MissileTerrainMaskParams p{};
    p.maxCoverTravelKm = 100.0;
    p.lookaheadSec = 0.0;
    p.coverPadKm = 0.01;
    p.ignoreIfAlreadyOccluded = true;

    const auto plan = sim::planMissileTerrainMasking(m, tgtPos, tgtVel, targets, 1, /*seed=*/7, p);

    if (!plan.valid) {
      std::cerr << "FAIL: expected valid terrain mask plan\n";
      fails++;
    } else {
      if (plan.asteroidId != 111) {
        std::cerr << "FAIL: expected asteroidId 111, got " << plan.asteroidId << "\n";
        fails++;
      }

      const double expectedTravel = (plan.coverPointKm - tgtPos).length();
      if (!approx(plan.travelKm, expectedTravel, 1e-9)) {
        std::cerr << "FAIL: travelKm mismatch. got=" << plan.travelKm
                  << " expected=" << expectedTravel << "\n";
        fails++;
      }

      const double dirLen = plan.dirWorld.length();
      if (!approx(dirLen, 1.0, 1e-9)) {
        std::cerr << "FAIL: expected dirWorld normalized, len=" << dirLen << "\n";
        fails++;
      }

      // Ensure the chosen cover point is actually in the asteroid shadow from the missile.
      if (!math::segmentHitsSphere(m.posKm, plan.coverPointKm, a.centerKm, a.radiusKm)) {
        std::cerr << "FAIL: expected cover point to be occluded by asteroid sphere\n";
        fails++;
      }
    }
  }

  // If the target is already occluded (and ignoreIfAlreadyOccluded is true), no plan is needed.
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 0};

    const math::Vec3d tgtPos{10, 0, 0};
    const math::Vec3d tgtVel{0, 0, 0};

    sim::SphereTarget a{};
    a.kind = sim::CombatTargetKind::Asteroid;
    a.id = 222;
    a.centerKm = {5, 0, 0};
    a.radiusKm = 2.0;

    const sim::SphereTarget targets[] = {a};

    const auto plan = sim::planMissileTerrainMasking(m, tgtPos, tgtVel, targets, 1, /*seed=*/9);
    if (plan.valid) {
      std::cerr << "FAIL: expected invalid plan when already occluded\n";
      fails++;
    }

    sim::MissileTerrainMaskParams p{};
    p.ignoreIfAlreadyOccluded = false;
    p.maxCoverTravelKm = 100.0;
    p.lookaheadSec = 0.0;

    const auto plan2 = sim::planMissileTerrainMasking(m, tgtPos, tgtVel, targets, 1, /*seed=*/9, p);
    if (!plan2.valid) {
      std::cerr << "FAIL: expected valid plan when ignoreIfAlreadyOccluded is false\n";
      fails++;
    }
  }

  return fails;
}
