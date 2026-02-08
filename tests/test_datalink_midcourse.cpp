#include "stellar/sim/Combat.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_datalink_midcourse() {
  int fails = 0;

  // Locked target.
  sim::SphereTarget tgt{};
  tgt.kind = sim::CombatTargetKind::Ship;
  tgt.index = 0;
  tgt.id = 1;
  tgt.centerKm = {1000, 0, 0};
  tgt.velKmS = {0, 1, 0};
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
  m.turnRateRadS = 0.0;

  m.hasTarget = true;
  m.targetKind = sim::CombatTargetKind::Ship;
  m.targetId = 1;

  // Keep the seeker inactive so only midcourse datalink logic is exercised.
  m.seeker = sim::MissileSeekerType::Radar;
  m.seekerFovCos = 0.0;
  m.seekerActivationSimSec = 10.0;

  // Enable discrete datalink updates with latency.
  //
  // Convention:
  //  - Initial guidance fix at t=0 (launch)
  //  - Subsequent update arrivals at t = latency + n*period
  //
  // With period=1.0 and latency=0.2, the first post-launch update arrives at t=0.2,
  // and the next at t=1.2.
  m.datalinkRangeKm = 0.0; // perfect range (no shooter needed)
  m.datalinkUpdatePeriodSimSec = 1.0;
  m.datalinkLatencySimSec = 0.2;

  // Intentionally leave targetMemorySimSec at the default (0) to verify that discrete
  // midcourse datalink updates still bridge between ticks.

  ms.push_back(m);

  std::vector<sim::MissileDetonation> dets;
  std::vector<sim::MissileHit> hits;

  // Step 0 -> 0.1: should ingest the initial launch fix immediately.
  {
    sim::SphereTarget targets[1]{tgt};
    sim::stepMissiles(ms, 0.1, targets, 1, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_datalink_midcourse] missile destroyed unexpectedly in first step.\n";
      return 1;
    }

    if (!ms[0].hasLastKnownTarget) {
      std::cerr << "[test_datalink_midcourse] expected lastKnownTarget to be seeded at launch.\n";
      ++fails;
    } else if (!approx(ms[0].lastKnownTargetVelKmS.y, 1.0, 1e-9)) {
      std::cerr << "[test_datalink_midcourse] unexpected initial lastKnownTargetVel.y. got="
                << ms[0].lastKnownTargetVelKmS.y << " expected 1.0\n";
      ++fails;
    }
  }

  // Advance time to t=0.6.
  // Change target velocity at t>=0.5. The missile should not see the new velocity
  // until the next datalink tick at t=1.2.
  double t = 0.1;
  sim::SphereTarget cur = tgt;

  for (int i = 0; i < 5; ++i) {
    if (t >= 0.5 - 1e-12) {
      cur.velKmS = {0, -10, 0};
    }
    sim::SphereTarget targets[1]{cur};
    sim::stepMissiles(ms, 0.1, targets, 1, dets, hits);
    t += 0.1;

    if (ms.empty()) {
      std::cerr << "[test_datalink_midcourse] missile destroyed unexpectedly while advancing to t=0.6.\n";
      return 1;
    }
  }

  if (!ms[0].hasLastKnownTarget) {
    std::cerr << "[test_datalink_midcourse] expected lastKnownTarget to persist between datalink ticks.\n";
    ++fails;
  } else {
    // The last datalink tick should have occurred at t=0.2, so by t=0.6 the age should be > 0.
    if (ms[0].lastKnownTargetAgeSimSec <= 0.0) {
      std::cerr << "[test_datalink_midcourse] expected lastKnownTargetAge to advance between ticks. age="
                << ms[0].lastKnownTargetAgeSimSec << "\n";
      ++fails;
    }

    // Velocity should remain the old sample until the 1.2s tick.
    if (!approx(ms[0].lastKnownTargetVelKmS.y, 1.0, 1e-9)) {
      std::cerr << "[test_datalink_midcourse] expected lastKnownTargetVel.y to remain 1.0 before the 1.2s tick. got="
                << ms[0].lastKnownTargetVelKmS.y << "\n";
      ++fails;
    }
  }

  // Advance to just after t=1.2 so the delayed update arrives.
  while (t < 1.25) {
    cur.velKmS = {0, -10, 0};
    sim::SphereTarget targets[1]{cur};
    sim::stepMissiles(ms, 0.1, targets, 1, dets, hits);
    t += 0.1;

    if (ms.empty()) {
      std::cerr << "[test_datalink_midcourse] missile destroyed unexpectedly while advancing to t=1.25.\n";
      return 1;
    }
  }

  if (!ms[0].hasLastKnownTarget) {
    std::cerr << "[test_datalink_midcourse] expected lastKnownTarget to remain after the delayed tick.\n";
    ++fails;
  } else if (!approx(ms[0].lastKnownTargetVelKmS.y, -10.0, 1e-9)) {
    std::cerr << "[test_datalink_midcourse] expected lastKnownTargetVel.y to update on the delayed tick. got="
              << ms[0].lastKnownTargetVelKmS.y << " expected -10\n";
    ++fails;
  }

  return fails;
}
