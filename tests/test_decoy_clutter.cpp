#include "stellar/sim/Combat.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_decoy_clutter() {
  int fails = 0;

  // --- Decoy clutter suppression should reduce trackQuality on a direct track ---
  //
  // We set trackQuality half-lives to 0 so the update is instantaneous and the
  // expected value is easy to reason about.
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 10.0;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = std::cos(0.6);
    m.seekerSlewRateRadS = 0.0; // simplify measurement quality

    m.enableTrackQuality = true;
    m.trackQuality = 1.0;
    m.trackQualityRiseHalfLifeSimSec = 0.0;
    m.trackQualityFallHalfLifeSimSec = 0.0;
    m.trackQualityResistFloor = 0.0;

    m.decoyClutterTrackSuppressionGain = 1.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 42;

    sim::SphereTarget ship{};
    ship.kind = sim::CombatTargetKind::Ship;
    ship.id = 42;
    ship.centerKm = {0, 0, 1000};
    ship.velKmS = {0, 0, 0};
    ship.radarSignature = 1.0;
    ship.radiusKm = 1.0;

    sim::SphereTarget chaff{};
    chaff.kind = sim::CombatTargetKind::Decoy;
    chaff.id = 9001;
    chaff.centerKm = {0, 0, 1000};
    chaff.velKmS = {0, 0, 0};
    chaff.decoyRadar = 1.0;
    chaff.radiusKm = 0.1;

    sim::SphereTarget targets[2]{ship, chaff};

    std::vector<sim::Missile> ms;
    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    sim::stepMissiles(ms, /*dtSimSec*/0.25, targets, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_decoy_clutter] missile unexpectedly expired\n";
      ++fails;
    } else {
      // With identical target/decoy score, clutter ratio ~= 1 and measQ becomes:
      //   1 / (1 + gain * 1) = 0.5
      if (!approx(ms[0].trackQuality, 0.5, 1e-6)) {
        std::cerr << "[test_decoy_clutter] expected trackQuality ~= 0.5 with clutter. got "
                  << ms[0].trackQuality << "\n";
        ++fails;
      }
    }
  }

  // --- Gain=0 should preserve legacy behavior (trackQuality stays 1 for perfect lock) ---
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 10.0;

    m.seeker = sim::MissileSeekerType::Radar;
    m.seekerActivationSimSec = 0.0;
    m.seekerFovCos = std::cos(0.6);
    m.seekerSlewRateRadS = 0.0;

    m.enableTrackQuality = true;
    m.trackQuality = 1.0;
    m.trackQualityRiseHalfLifeSimSec = 0.0;
    m.trackQualityFallHalfLifeSimSec = 0.0;

    m.decoyClutterTrackSuppressionGain = 0.0;

    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 42;

    sim::SphereTarget ship{};
    ship.kind = sim::CombatTargetKind::Ship;
    ship.id = 42;
    ship.centerKm = {0, 0, 1000};
    ship.velKmS = {0, 0, 0};
    ship.radarSignature = 1.0;
    ship.radiusKm = 1.0;

    sim::SphereTarget chaff{};
    chaff.kind = sim::CombatTargetKind::Decoy;
    chaff.id = 9001;
    chaff.centerKm = {0, 0, 1000};
    chaff.velKmS = {0, 0, 0};
    chaff.decoyRadar = 1.0;
    chaff.radiusKm = 0.1;

    sim::SphereTarget targets[2]{ship, chaff};

    std::vector<sim::Missile> ms;
    ms.push_back(m);

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;

    sim::stepMissiles(ms, /*dtSimSec*/0.25, targets, 2, dets, hits);

    if (ms.empty()) {
      std::cerr << "[test_decoy_clutter] missile unexpectedly expired (gain=0)\n";
      ++fails;
    } else if (!approx(ms[0].trackQuality, 1.0, 1e-9)) {
      std::cerr << "[test_decoy_clutter] expected trackQuality == 1 with gain=0. got "
                << ms[0].trackQuality << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_decoy_clutter] PASS\n";
  }
  return fails;
}
