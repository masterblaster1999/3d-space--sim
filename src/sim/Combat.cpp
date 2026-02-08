#include "stellar/sim/Combat.h"

#include "stellar/core/Random.h"

#include "stellar/math/Geometry.h"
#include "stellar/math/Math.h"
#include "stellar/sim/Ballistics.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

bool raySphereIntersectKm(const math::Vec3d& originKm,
                          const math::Vec3d& dirNormalized,
                          const math::Vec3d& centerKm,
                          double radiusKm,
                          double& outTEnterKm) {
  // Delegate the actual math to the shared geometry helper so ray/sphere behavior
  // stays consistent across sim modules.
  return math::raySphereIntersectEnter(originKm,
                                       dirNormalized,
                                       centerKm,
                                       std::max(0.0, radiusKm),
                                       outTEnterKm);
}

RaycastHit raycastNearestSphereKm(const math::Vec3d& originKm,
                                  const math::Vec3d& dirNormalized,
                                  double maxRangeKm,
                                  const SphereTarget* targets,
                                  std::size_t targetCount) {
  RaycastHit best{};
  best.hit = false;
  best.tKm = std::max(0.0, maxRangeKm);
  best.pointKm = originKm + dirNormalized * best.tKm;

  if (!targets || targetCount == 0) return best;

  const double maxR = std::max(0.0, maxRangeKm);
  double bestT = maxR;

  for (std::size_t i = 0; i < targetCount; ++i) {
    const SphereTarget& t = targets[i];
    const math::Vec3d toC = t.centerKm - originKm;
    const double dist2 = toC.lengthSq();
    if (dist2 < 1e-12) continue;

    // Aim cone filter (optional).
    if (t.minAimCos > -0.5) {
      const double dist = std::sqrt(dist2);
      const double aim = math::dot(dirNormalized, toC / dist);
      if (aim < t.minAimCos) continue;
    }

    // Cheap early reject: center projection outside [0, maxRange + radius].
    const double tProj = math::dot(toC, dirNormalized);
    if (tProj < -t.radiusKm) continue;
    if (tProj > maxR + t.radiusKm) continue;

    double tEnter = 0.0;
    if (!raySphereIntersectKm(originKm, dirNormalized, t.centerKm, t.radiusKm, tEnter)) continue;
    if (tEnter > maxR) continue;

    if (tEnter < bestT) {
      bestT = tEnter;
      best.hit = true;
      best.tKm = tEnter;
      best.pointKm = originKm + dirNormalized * tEnter;
      best.kind = t.kind;
      best.index = t.index;
      best.id = t.id;
    }
  }

  if (!best.hit) {
    best.tKm = maxR;
    best.pointKm = originKm + dirNormalized * maxR;
  }
  return best;
}

bool segmentHitsSphereKm(const math::Vec3d& aKm,
                         const math::Vec3d& bKm,
                         const math::Vec3d& centerKm,
                         double radiusKm) {
  return math::segmentHitsSphere(aKm, bKm, centerKm, std::max(0.0, radiusKm));
}

static bool occludedByAsteroidSpheres(const math::Vec3d& fromKm,
                                     const math::Vec3d& toKm,
                                     const SphereTarget* targets,
                                     std::size_t targetCount,
                                     double radiusPadKm) {
  if (!targets || targetCount == 0) return false;

  const double pad = std::max(0.0, radiusPadKm);

  for (std::size_t i = 0; i < targetCount; ++i) {
    const SphereTarget& t = targets[i];
    if (t.kind != CombatTargetKind::Asteroid) continue;

    const double r = std::max(0.0, t.radiusKm) + pad;
    if (r <= 0.0) continue;

    if (math::segmentHitsSphere(fromKm, toKm, t.centerKm, r)) {
      return true;
    }
  }

  return false;
}

static double heatSeekerAspectFactor(const Missile& m,
                                     const math::Vec3d& targetVelKmS,
                                     const math::Vec3d& toTargetDirNorm) {
  const double front = std::max(0.0, m.heatAspectFrontFactor);
  const double rear  = std::max(0.0, m.heatAspectRearFactor);

  // Fast path: defaults preserve legacy behavior (no aspect weighting).
  if (std::fabs(front - 1.0) < 1e-12 && std::fabs(rear - 1.0) < 1e-12) {
    return 1.0;
  }

  const double minSpd = std::max(0.0, m.heatAspectMinSpeedKmS);
  const double spd = targetVelKmS.length();
  if (!(spd > minSpd)) return 1.0;

  const double fullSpd = std::max(minSpd + 1e-6, m.heatAspectSpeedForFullKmS);
  double spd01 = (spd - minSpd) / (fullSpd - minSpd);
  spd01 = std::clamp(spd01, 0.0, 1.0);

  math::Vec3d velDir = targetVelKmS / spd;
  double aspectCos = math::dot(velDir, toTargetDirNorm);
  aspectCos = std::clamp(aspectCos, -1.0, 1.0);

  // aspectCos = -1 => head-on, +1 => tail chase.
  const double t = 0.5 * (aspectCos + 1.0);
  const double full = front * (1.0 - t) + rear * t;

  // Blend in by speed so near-stationary targets behave like the old model.
  return 1.0 + (full - 1.0) * spd01;
}

static double effectiveHeatSignature(const Missile& m,
                                     double baseSignature,
                                     const math::Vec3d& targetVelKmS,
                                     const math::Vec3d& toTargetDirNorm) {
  const double sig = std::max(0.0, baseSignature);
  if (!(sig > 0.0)) return 0.0;
  if (toTargetDirNorm.lengthSq() < 1e-12) return sig;
  return sig * heatSeekerAspectFactor(m, targetVelKmS, toTargetDirNorm);
}

static double effectiveRadarSignature(double baseSignature) {
  return std::max(0.0, baseSignature);
}

static double receivedJamming01(const Missile& m, double distKm, double jammerPower) {
  jammerPower = std::max(0.0, jammerPower);
  if (!(jammerPower > 0.0)) return 0.0;

  const double halfRangeKm = std::max(0.0, m.radarJammingHalfRangeKm);
  if (halfRangeKm <= 0.0) {
    // Back-compat / simple mode: interpret jammerPower directly as a [0,1] level.
    return std::clamp(jammerPower, 0.0, 1.0);
  }

  const double minDistKm = std::max(1.0e-6, m.radarJammingMinDistKm);
  distKm = std::max(distKm, minDistKm);

  const double r2 = halfRangeKm * halfRangeKm;
  const double d2 = distKm * distKm;
  const double snr = jammerPower * (r2 / d2);
  if (!std::isfinite(snr) || snr <= 0.0) return 0.0;

  const double j01 = snr / (1.0 + snr);
  return std::clamp(j01, 0.0, 1.0);
}



// Optional: received seeker signal level in [0,1] (distance falloff).
//
// This is a deliberately lightweight SNR model intended to modulate seeker
// track quality without adding random noise. It uses the same stable
// `snr/(1+snr)` mapping as our radar jamming reception helper.
static double receivedSeekerSignal01(const Missile& miss, double distKm, double signalStrength) {
  const double halfRangeKm = std::max(0.0, miss.seekerSignalHalfRangeKm);
  if (halfRangeKm <= 0.0) {
    // Disabled: preserve legacy behavior (perfect signal).
    return 1.0;
  }

  signalStrength = std::max(0.0, signalStrength);
  if (!(signalStrength > 0.0)) return 0.0;

  const double minDistKm = std::max(1.0e-6, miss.seekerSignalMinDistKm);
  distKm = std::max(distKm, minDistKm);

  const double r2 = halfRangeKm * halfRangeKm;
  const double d2 = distKm * distKm;
  const double snr = signalStrength * (r2 / d2);
  if (!std::isfinite(snr) || snr <= 0.0) return 0.0;

  const double s01 = snr / (1.0 + snr);
  return std::clamp(s01, 0.0, 1.0);
}

static bool blastOccludedByAsteroidSpheres(const math::Vec3d& fromKm,
                                          const math::Vec3d& toKm,
                                          const SphereTarget* targets,
                                          std::size_t targetCount,
                                          double radiusPadKm,
                                          core::u64 ignoreAsteroidId) {
  if (!targets || targetCount == 0) return false;

  const double pad = std::max(0.0, radiusPadKm);
  constexpr double kEps = 1.0e-6;

  for (std::size_t i = 0; i < targetCount; ++i) {
    const SphereTarget& t = targets[i];
    if (t.kind != CombatTargetKind::Asteroid) continue;
    if (ignoreAsteroidId != 0 && t.id == ignoreAsteroidId) continue;

    const double r = std::max(0.0, t.radiusKm) + pad;
    if (r <= 0.0) continue;

    double tEnter01 = 0.0;
    double tExit01 = 0.0;
    if (!math::segmentSphereIntersectionT(fromKm, toKm, t.centerKm, r, tEnter01, tExit01)) {
      continue;
    }

    // Treat as an occluder only when the segment actually passes through a
    // non-trivial portion of the asteroid volume (not just touching at an
    // endpoint).
    const double e = std::clamp(tEnter01, 0.0, 1.0);
    const double x = std::clamp(tExit01, 0.0, 1.0);

    if (x <= kEps) continue;            // touch at start
    if (e >= 1.0 - kEps) continue;      // touch at end
    if (x - e <= kEps) continue;        // tangent / negligible thickness

    return true;
  }

  return false;
}



static math::Vec3d rotateTowards(const math::Vec3d& fromDir,
                                 const math::Vec3d& toDir,
                                 double maxAngleRad) {
  const math::Vec3d f = fromDir.normalized();
  const math::Vec3d t = toDir.normalized();
  if (f.lengthSq() < 1e-12) return t;
  if (t.lengthSq() < 1e-12) return f;

  const double c = std::clamp(math::dot(f, t), -1.0, 1.0);
  const double ang = std::acos(c);
  if (ang <= 1e-9) return t;
  if (maxAngleRad <= 0.0) return f;
  if (ang <= maxAngleRad) return t;

  math::Vec3d axis = math::cross(f, t);
  const double al2 = axis.lengthSq();
  if (al2 < 1e-18) {
    // Nearly opposite; just keep current direction.
    return f;
  }
  axis = axis / std::sqrt(al2);

  const double ca = std::cos(maxAngleRad);
  const double sa = std::sin(maxAngleRad);

  // Rodrigues rotation formula.
  const math::Vec3d v = f;
  const math::Vec3d v2 = math::cross(axis, v);
  return (v * ca + v2 * sa).normalized();
}

static math::Vec3d rotateAroundAxis(const math::Vec3d& v,
                                   const math::Vec3d& axisUnit,
                                   double angleRad) {
  // Rodrigues rotation formula; assumes axisUnit is normalized.
  const double ca = std::cos(angleRad);
  const double sa = std::sin(angleRad);
  const math::Vec3d a = axisUnit;
  const math::Vec3d v2 = math::cross(a, v);
  return (v * ca + v2 * sa + a * (math::dot(a, v) * (1.0 - ca)));
}

static void makePerpBasisFromForward(const math::Vec3d& forwardWorld,
                                    math::Vec3d& outRight,
                                    math::Vec3d& outUp) {
  math::Vec3d f = forwardWorld;
  if (f.lengthSq() <= 1e-12) f = {0, 0, 1};
  f = f.normalized();

  math::Vec3d up = {0, 1, 0};
  if (std::abs(math::dot(up, f)) > 0.95) up = {1, 0, 0};

  math::Vec3d r = math::cross(up, f);
  if (r.lengthSq() <= 1e-12) {
    up = {0, 0, 1};
    if (std::abs(math::dot(up, f)) > 0.95) up = {1, 0, 0};
    r = math::cross(up, f);
  }

  r = math::safeNormalized(r, {1, 0, 0}, 1e-18);
  math::Vec3d u = math::cross(f, r);
  u = math::safeNormalized(u, {0, 1, 0}, 1e-18);

  outRight = r;
  outUp = u;
}

static math::Vec3d applyTerminalWeaveDir(Missile& m,
                                        const math::Vec3d& baseDir,
                                        const SphereTarget* guide,
                                        const math::Vec3d& missilePosKm,
                                        double speedKmS,
                                        double ageSimSec,
                                        bool seekerActive) {
  const double amp = std::max(0.0, m.terminalWeaveAmplitudeRad);
  const double hz  = std::max(0.0, m.terminalWeaveFrequencyHz);

  if (amp <= 0.0 || hz <= 0.0) return baseDir;
  if (!seekerActive) return baseDir;
  if (!guide) return baseDir;

  speedKmS = std::max(0.0, speedKmS);
  if (!(speedKmS > 1e-9)) return baseDir;

  const math::Vec3d to = guide->centerKm - missilePosKm;
  const double dist = std::sqrt(std::max(0.0, to.lengthSq()));
  const double tti  = dist / speedKmS;

  const double startTti = std::max(0.0, m.terminalWeaveStartTtiSec);
  if (startTti > 1e-9 && tti > startTti) return baseDir;

  double ramp = 1.0;
  if (startTti > 1e-9) {
    ramp = std::clamp(1.0 - tti / startTti, 0.0, 1.0);
  }

  double ampEff = amp * ramp;
  const double clampRad = std::max(0.0, m.terminalWeaveMaxRad);
  if (clampRad > 0.0) ampEff = std::min(ampEff, clampRad);
  if (!(ampEff > 1e-9)) return baseDir;

  // Seed a deterministic phase offset once.
  if (!m.hasTerminalWeavePhase) {
    core::u64 seed = core::hashCombine(m.shooterId, m.targetId);
    seed = core::hashCombine(seed,
                             static_cast<core::u64>(m.fromPlayer ? 0xA5A5A5A5ull : 0x5A5A5A5Aull));
    core::SplitMix64 rng(seed);
    m.terminalWeavePhaseRad = rng.nextDouble() * (2.0 * math::kPi);
    m.hasTerminalWeavePhase = true;
  }

  const double phase = m.terminalWeavePhaseRad + (2.0 * math::kPi) * hz * std::max(0.0, ageSimSec);

  math::Vec3d f = baseDir;
  if (f.lengthSq() <= 1e-18) f = {0, 0, 1};
  else f = f.normalized();

  math::Vec3d right{1, 0, 0}, up{0, 1, 0};
  makePerpBasisFromForward(f, right, up);

  const math::Vec3d lateral = right * std::cos(phase) + up * std::sin(phase);

  const double ca = std::cos(ampEff);
  const double sa = std::sin(ampEff);
  math::Vec3d out = f * ca + lateral * sa;
  if (out.lengthSq() <= 1e-18) return f;
  return out.normalized();
}


static math::Vec3d applyTrackErrorDir(Missile& m,
                                      const math::Vec3d& baseDir,
                                      double ageSimSec,
                                      bool seekerActive,
                                      const SphereTarget* guide) {
  const double maxErr = std::max(0.0, m.trackErrorMaxRad);
  const double hz = std::max(0.0, m.trackErrorFrequencyHz);

  if (maxErr <= 0.0 || hz <= 0.0) return baseDir;
  if (!seekerActive) return baseDir;
  if (!guide) return baseDir;

  // Decoy tracks are already "bright"; do not apply additional sensor error.
  if (guide->kind == CombatTargetKind::Decoy) return baseDir;

  // Track error is driven by track quality. If the feature is disabled, preserve legacy behavior.
  if (!m.enableTrackQuality) return baseDir;

  const double q = std::clamp(m.trackQuality, 0.0, 1.0);
  const double amp = maxErr * (1.0 - q);
  if (!(amp > 1e-9)) return baseDir;

  // Seed a deterministic phase offset once.
  if (!m.hasTrackErrorPhase) {
    core::u64 seed = core::hashCombine(m.shooterId, m.targetId);
    seed = core::hashCombine(seed,
                             static_cast<core::u64>(m.fromPlayer ? 0x13579BDFull : 0x2468ACE0ull));
    seed = core::hashCombine(seed, 0xC001D00Dull);
    core::SplitMix64 rng(seed);
    m.trackErrorPhaseRad = rng.nextDouble() * (2.0 * math::kPi);
    m.hasTrackErrorPhase = true;
  }

  const double phase = m.trackErrorPhaseRad + (2.0 * math::kPi) * hz * std::max(0.0, ageSimSec);

  math::Vec3d f = baseDir;
  if (f.lengthSq() <= 1e-18) f = {0, 0, 1};
  else f = f.normalized();

  math::Vec3d right{1, 0, 0}, up{0, 1, 0};
  makePerpBasisFromForward(f, right, up);

  math::Vec3d lateral = right * std::cos(phase) + up * std::sin(phase);
  if (!(lateral.lengthSq() > 1e-18)) return f;
  lateral = lateral.normalized();

  const double ca = std::cos(amp);
  const double sa = std::sin(amp);
  math::Vec3d out = f * ca + lateral * sa;
  if (!(out.lengthSq() > 1e-18)) return f;
  return out.normalized();
}

static math::Vec3d applySeekerLagGuidanceDir(const Missile& m,
                                             const math::Vec3d& vHat,
                                             double maxTurnRad,
                                             const math::Vec3d& cmdDir,
                                             bool seekerActive,
                                             bool seekerUpdateThisStep,
                                             const SphereTarget* lockedTarget,
                                             const SphereTarget* guide,
                                             double seekerFovCos) {
  const double gain = std::max(0.0, m.seekerLagGuidanceGain);
  if (gain <= 0.0) return cmdDir;
  if (!seekerActive) return cmdDir;
  if (!seekerUpdateThisStep) return cmdDir;
  if (!(m.seekerSlewRateRadS > 0.0)) return cmdDir;

  // Only apply when we have a direct seeker measurement of our locked target.
  if (!lockedTarget || !guide || guide != lockedTarget) return cmdDir;
  if (lockedTarget->kind == CombatTargetKind::Decoy) return cmdDir;
  if (!m.hasSeekerDir) return cmdDir;

  math::Vec3d sDir = m.seekerDirWorld;
  if (sDir.lengthSq() < 1e-12) return cmdDir;
  sDir = sDir.normalized();

  const math::Vec3d to = lockedTarget->centerKm - m.posKm;
  if (to.lengthSq() < 1e-12) return cmdDir;
  const math::Vec3d toDir = to.normalized();

  // Misalignment measure in [0,1] using the same "distance from FOV edge" mapping
  // used by the track-quality gimbal metric.
  const double align = std::clamp(math::dot(sDir, toDir), -1.0, 1.0);

  seekerFovCos = std::clamp(seekerFovCos, -1.0, 1.0);
  const double denom = 1.0 - seekerFovCos;

  double gimbalQ = 1.0;
  if (denom > 1e-6) {
    gimbalQ = std::clamp((align - seekerFovCos) / denom, 0.0, 1.0);
  } else {
    // Extremely narrow cone: treat perfect alignment as "good".
    gimbalQ = (align >= seekerFovCos - 1e-9) ? 1.0 : 0.0;
  }

  const double mis01 = 1.0 - gimbalQ;
  const double w = std::clamp(gain * mis01, 0.0, 1.0);
  if (!(w > 1e-9)) return cmdDir;

  math::Vec3d base = cmdDir;
  if (base.lengthSq() < 1e-12) base = vHat;
  base = base.normalized();

  math::Vec3d blended = base * (1.0 - w) + sDir * w;
  blended = math::safeNormalized(blended, base, 1e-12);

  // Ensure we don't violate the missile's per-step turn limit.
  if (maxTurnRad > 0.0) {
    blended = rotateTowards(vHat, blended, maxTurnRad);
  }
  return blended;
}

static double seekerSlewPenaltyFactor(const Missile& m,
                                     double dtSim,
                                     const math::Vec3d& seekerDirNorm,
                                     const math::Vec3d& toDirNorm) {
  const double slew = std::max(0.0, m.seekerSlewRateRadS);
  if (slew <= 0.0 || dtSim <= 0.0) return 1.0;

  const double c = std::clamp(math::dot(seekerDirNorm, toDirNorm), -1.0, 1.0);
  const double ang = std::acos(c);
  const double maxAng = slew * dtSim;

  if (ang <= maxAng + 1.0e-9) return 1.0;

  const double ratio = maxAng / std::max(ang, 1.0e-6);

  // Track-quality coupling: a firm lock makes it harder to yank the seeker
  // across large angular jumps; poor track quality allows wider excursions.
  double q = 1.0;
  if (m.enableTrackQuality) {
    q = std::clamp(m.trackQuality, 0.0, 1.0);
  }
  const double p = 0.75 + 1.25 * q; // 0.75..2.0

  return std::pow(std::clamp(ratio, 0.0, 1.0), p);
}

static double effectiveTurnRateRadS(const Missile& m, double speedKmS) {
  double tr = std::max(0.0, m.turnRateRadS);
  if (tr <= 0.0) return 0.0;

  const double maxLat = std::max(0.0, m.maxLateralAccelKmS2);
  if (maxLat > 0.0 && speedKmS > 1e-9) {
    tr = std::min(tr, maxLat / speedKmS);
  }
  return tr;
}

static double guidanceSpeedKmS(const Missile& m, double currentSpeedKmS) {
  double s = std::max(0.0, currentSpeedKmS);
  const double accel = std::max(0.0, m.thrustAccelKmS2);
  const double maxSpeed = std::max(0.0, m.maxSpeedKmS);
  const double burn = std::max(0.0, m.motorBurnRemainingSimSec);

  if (accel > 0.0 && burn > 0.0 && maxSpeed > s + 1e-12) {
    // Small look-ahead so lead pursuit doesn't consistently over-lead missiles
    // that will accelerate immediately after launch.
    const double horizon = std::min(burn, 1.0);
    s = std::min(maxSpeed, s + accel * horizon);
  }
  return s;
}

static math::Vec3d applyAsteroidAvoidanceDir(const Missile& m,
                                              const math::Vec3d& missilePosKm,
                                              const math::Vec3d& vHat,
                                              double speedKmS,
                                              double turnRateRadS,
                                              double dtSim,
                                              const math::Vec3d& cmdDir,
                                              const SphereTarget* targets,
                                              std::size_t targetCount) {
  const double strength = std::max(0.0, m.asteroidAvoidanceStrength);
  if (strength <= 0.0) return cmdDir;
  if (!targets || targetCount == 0) return cmdDir;

  // Don't avoid the target if the missile is explicitly aimed at an asteroid.
  const bool targetIsAsteroid = (m.hasTarget && m.targetKind == CombatTargetKind::Asteroid);

  // Prediction horizon: default to ~1s when not specified.
  double horizonSim = m.asteroidAvoidanceLookaheadSimSec;
  if (!(horizonSim > 0.0)) horizonSim = 1.0;
  horizonSim = std::max(horizonSim, std::max(0.0, dtSim));

  speedKmS = std::max(0.0, speedKmS);
  const double lookDistKm = speedKmS * horizonSim;
  if (!(lookDistKm > 1e-9)) return cmdDir;

  const double pad = std::max(0.0, m.asteroidAvoidancePadKm);
  const double missileR = std::max(0.0, m.radiusKm);

  math::Vec3d accum{0, 0, 0};

  // A small, deterministic repulsion field: for asteroids that lie close to the
  // current velocity ray within the lookahead distance, accumulate a lateral
  // "push" direction away from the sphere.
  for (std::size_t i = 0; i < targetCount; ++i) {
    const SphereTarget& t = targets[i];
    if (t.kind != CombatTargetKind::Asteroid) continue;
    if (targetIsAsteroid && t.id == m.targetId) continue;

    const double r = std::max(0.0, t.radiusKm) + pad + missileR;
    if (!(r > 0.0)) continue;

    const math::Vec3d rel = t.centerKm - missilePosKm;

    // Project onto the missile's current travel ray.
    const double along = math::dot(rel, vHat);
    if (!(along > 0.0)) continue;
    if (along > lookDistKm + r) continue;

    const math::Vec3d lateral = rel - vHat * along;
    const double d2 = lateral.lengthSq();
    const double r2 = r * r;
    if (d2 > r2) continue;

    const double d = std::sqrt(std::max(0.0, d2));
    const double invR = 1.0 / std::max(1e-6, r);

    // Weight grows rapidly as we approach the collision tube.
    const double kDist = std::clamp(1.0 - d * invR, 0.0, 1.0);
    const double kAlong = std::clamp(1.0 - along / std::max(1e-6, lookDistKm + r), 0.0, 1.0);
    const double w = (kDist * kDist) * (0.25 + 0.75 * kAlong);
    if (!(w > 0.0)) continue;

    // Away direction is perpendicular to the travel ray when possible.
    math::Vec3d away = {-lateral.x, -lateral.y, -lateral.z};
    if (!(away.lengthSq() > 1e-18)) {
      // Degenerate: if we're aimed dead-center, choose a stable perpendicular.
      away = math::cross(vHat, cmdDir);
      if (!(away.lengthSq() > 1e-18)) away = math::cross(vHat, math::Vec3d{0, 1, 0});
      if (!(away.lengthSq() > 1e-18)) away = math::cross(vHat, math::Vec3d{1, 0, 0});
    }
    away = math::safeNormalized(away, math::Vec3d{1, 0, 0});

    accum = accum + away * w;
  }

  if (!(accum.lengthSq() > 1e-18)) return cmdDir;

  const math::Vec3d push = accum.normalized();
  math::Vec3d base = cmdDir;
  if (!(base.lengthSq() > 1e-18)) base = vHat;
  base = base.normalized();

  // Blend toward a direction that includes a lateral push.
  math::Vec3d desired = base + push * strength;
  if (!(desired.lengthSq() > 1e-18)) return base;
  desired = desired.normalized();

  // Re-limit by max turn for stability.
  const double maxTurn = std::max(0.0, turnRateRadS) * std::max(0.0, dtSim);
  if (maxTurn <= 0.0) return vHat;

  return rotateTowards(vHat, desired, maxTurn);
}

struct MissileSwarmSnapshot {
  bool alive{false};
  bool fromPlayer{false};
  core::u64 shooterId{0};

  bool hasTarget{false};
  CombatTargetKind targetKind{CombatTargetKind::Ship};
  core::u64 targetId{0};

  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velDir{0, 0, 1};
};

static math::Vec3d applySwarmBiasDir(const Missile& m,
                                    std::size_t missileIndex,
                                    const std::vector<MissileSwarmSnapshot>& snap,
                                    const math::Vec3d& vHat,
                                    const math::Vec3d& cmdDir,
                                    double turnRateRadS,
                                    double dtSim) {
  const double sepStrength = std::max(0.0, m.swarmSeparationStrength);
  const double cohStrength = std::max(0.0, m.swarmCohesionStrength);
  const double aliStrength = std::max(0.0, m.swarmAlignmentStrength);
  if (sepStrength <= 0.0 && cohStrength <= 0.0 && aliStrength <= 0.0) return cmdDir;

  const double sepKm = std::max(0.0, m.swarmSeparationKm);
  if (sepStrength > 0.0 && sepKm <= 1e-12) {
    // No meaningful separation distance configured.
    return cmdDir;
  }

  double neighborRangeKm = std::max(0.0, m.swarmNeighborRangeKm);
  if (!(neighborRangeKm > 0.0)) {
    neighborRangeKm = (sepKm > 0.0) ? (4.0 * sepKm) : 0.0;
  }
  if (!(neighborRangeKm > 0.0)) return cmdDir;

  if (!m.hasTarget) return cmdDir;
  if (missileIndex >= snap.size()) return cmdDir;

  const MissileSwarmSnapshot& self = snap[missileIndex];
  if (!self.alive) return cmdDir;

  math::Vec3d sep{0, 0, 0};
  math::Vec3d sumPos{0, 0, 0};
  math::Vec3d sumVel{0, 0, 0};
  int count = 0;

  const double n2 = neighborRangeKm * neighborRangeKm;
  const double s2 = sepKm * sepKm;

  for (std::size_t j = 0; j < snap.size(); ++j) {
    if (j == missileIndex) continue;
    const MissileSwarmSnapshot& other = snap[j];
    if (!other.alive) continue;

    // Only coordinate with missiles from the same shooter when possible.
    if (m.shooterId != 0 && other.shooterId != 0) {
      if (other.shooterId != m.shooterId) continue;
    } else {
      if (other.fromPlayer != m.fromPlayer) continue;
    }

    // Only coordinate when pursuing the same target (prevents cross-target coupling).
    if (!other.hasTarget) continue;
    if (other.targetKind != m.targetKind) continue;
    if (other.targetId != m.targetId) continue;

    const math::Vec3d d = self.posKm - other.posKm;
    const double distSq = d.lengthSq();
    if (distSq < 1e-18) continue;
    if (distSq > n2) continue;

    ++count;
    sumPos = sumPos + other.posKm;
    sumVel = sumVel + other.velDir;

    if (sepStrength > 0.0 && distSq < s2) {
      const double dist = std::sqrt(distSq);
      const double k = std::clamp(1.0 - dist / std::max(1e-9, sepKm), 0.0, 1.0);
      sep = sep + (d / dist) * k;
    }
  }

  if (count <= 0) return cmdDir;

  math::Vec3d bias{0, 0, 0};

  if (sepStrength > 0.0 && sep.lengthSq() > 1e-18) {
    bias = bias + sep.normalized() * sepStrength;
  }

  if (cohStrength > 0.0) {
    const math::Vec3d center = sumPos * (1.0 / (double)count);
    math::Vec3d toCenter = center - self.posKm;
    toCenter = toCenter - vHat * math::dot(toCenter, vHat);
    if (toCenter.lengthSq() > 1e-18) {
      bias = bias + toCenter.normalized() * cohStrength;
    }
  }

  if (aliStrength > 0.0 && sumVel.lengthSq() > 1e-18) {
    math::Vec3d align = sumVel.normalized();
    align = align - vHat * math::dot(align, vHat);
    if (align.lengthSq() > 1e-18) {
      bias = bias + align.normalized() * aliStrength;
    }
  }

  // Keep bias lateral so it doesn't act like a speed change.
  bias = bias - vHat * math::dot(bias, vHat);
  if (!(bias.lengthSq() > 1e-18)) return cmdDir;

  math::Vec3d base = cmdDir;
  if (!(base.lengthSq() > 1e-18)) base = vHat;
  base = base.normalized();

  math::Vec3d desired = base + bias;
  if (!(desired.lengthSq() > 1e-18)) return base;
  desired = desired.normalized();

  double maxSteer = m.swarmMaxSteerRad;
  if (!(maxSteer > 0.0)) {
    // Default: use a fraction of per-step turn capability so it doesn't dominate guidance.
    maxSteer = 0.45 * std::max(0.0, turnRateRadS) * std::max(0.0, dtSim);
  }
  if (!(maxSteer > 0.0)) return base;

  return rotateTowards(base, desired, maxSteer);
}

static double estimateBoostCoastTtlSimSec(double rangeKm,
                                         double speed0KmS,
                                         double accelKmS2,
                                         double maxSpeedKmS,
                                         double burnTimeSimSec) {
  rangeKm = std::max(0.0, rangeKm);
  speed0KmS = std::max(1e-6, speed0KmS);
  accelKmS2 = std::max(0.0, accelKmS2);
  maxSpeedKmS = std::max(0.0, maxSpeedKmS);
  burnTimeSimSec = std::max(0.0, burnTimeSimSec);

  if (rangeKm <= 0.0) return 0.0;

  // Default: constant speed.
  if (accelKmS2 <= 0.0 || burnTimeSimSec <= 0.0 || maxSpeedKmS <= speed0KmS + 1e-12) {
    return rangeKm / speed0KmS;
  }

  const double tb = burnTimeSimSec;
  const double tToMax = std::max(0.0, (maxSpeedKmS - speed0KmS) / accelKmS2);

  // Accelerate up to either max speed or the end of burn.
  const double tAccel = std::min(tb, tToMax);
  const double dAccel = speed0KmS * tAccel + 0.5 * accelKmS2 * tAccel * tAccel;

  double dBurn = dAccel;
  double speedEnd = speed0KmS + accelKmS2 * tAccel;

  // If we hit max speed before burn ends, cruise at max speed for the remainder of the burn.
  if (tToMax < tb) {
    speedEnd = maxSpeedKmS;
    dBurn += maxSpeedKmS * (tb - tToMax);
  }

  // Range reached during burn.
  if (rangeKm <= dBurn + 1e-9) {
    // Within the acceleration portion.
    if (rangeKm <= dAccel + 1e-9) {
      // Solve: 0.5*a*t^2 + v0*t - range = 0
      const double a = 0.5 * accelKmS2;
      const double b = speed0KmS;
      const double c = -rangeKm;
      const double disc = b * b - 4.0 * a * c;
      if (disc <= 0.0) return 0.0;
      const double t = (-b + std::sqrt(disc)) / (2.0 * a);
      return std::max(0.0, t);
    }

    // Within the max-speed cruise portion of the burn.
    // This only happens when we reach max speed.
    const double remainKm = std::max(0.0, rangeKm - dAccel);
    return tToMax + remainKm / maxSpeedKmS;
  }

  // Coast after burn.
  const double remainKm = rangeKm - dBurn;
  return tb + remainKm / std::max(1e-6, speedEnd);
}

// Simple 3D proportional navigation steering step.
//
// This is an "energy conserving" PN variant that produces an acceleration command
// normal to the missile velocity (good fit for our constant-speed steering model).
//
// Reference (notation):
//  - R = Rt - Rm
//  - Vr = Vt - Vm
//  - Omega = (R x Vr) / |R|^2   (line-of-sight rotation vector)
//  - Vc = -dot(R, Vr) / |R|     (closing speed along LOS; positive when closing)
//  - a = N * Vc * (Omega x Vm_hat)
//
// When the PN geometry becomes degenerate (opening targets / tiny range / etc),
// we fall back to a simple turn toward the provided fallback aim direction.
static math::Vec3d proNavSteerDir(const math::Vec3d& vHat,
                                 const math::Vec3d& missileVelKmS,
                                 double missileSpeedKmS,
                                 double turnRateRadS,
                                 double dtSim,
                                 const math::Vec3d& missilePosKm,
                                 const SphereTarget& target,
                                 double navConstant,
                                 const math::Vec3d& targetAccelKmS2,
                                 double apnTargetAccelGain,
                                 double apnMaxTargetAccelKmS2,
                                 const math::Vec3d& fallbackAimDir) {
  const double sp = std::max(1e-9, missileSpeedKmS);
  const double N = std::clamp(navConstant, 0.0, 12.0);

  // Normalize the fallback direction defensively.
  math::Vec3d fb = fallbackAimDir;
  if (fb.lengthSq() < 1e-12) fb = vHat;
  fb = fb.normalized();

  // No steering authority -> keep current heading.
  const double maxTurnRate = std::max(0.0, turnRateRadS);
  if (!(dtSim > 0.0) || maxTurnRate <= 0.0) {
    return vHat;
  }

  const double maxTurn = maxTurnRate * dtSim;

  // PN disabled -> just turn toward the fallback aim direction.
  if (N <= 0.0) {
    return rotateTowards(vHat, fb, maxTurn);
  }

  const math::Vec3d R = target.centerKm - missilePosKm;
  const double r2 = R.lengthSq();
  if (!(r2 > 1e-9)) {
    return rotateTowards(vHat, fb, maxTurn);
  }
  const double r = std::sqrt(r2);

  const math::Vec3d Vr = target.velKmS - missileVelKmS;
  const double closing = -math::dot(R, Vr) / r;  // positive when range is decreasing

  // Opening targets (or near-zero closing) are numerically ugly for PN; fall back to
  // a simple lead/pursuit direction.
  if (!(closing > 1e-9)) {
    return rotateTowards(vHat, fb, maxTurn);
  }

  // Line-of-sight rotation vector.
  const math::Vec3d omega = math::cross(R, Vr) / r2;

  // Energy-conserving PN acceleration (normal to missile velocity).
  // a = N * Vc * (Omega x Vm_hat)
  math::Vec3d aCmd = math::cross(omega, vHat) * (N * closing);

  // Optional APN (augmented proportional navigation) feed-forward:
  // add a component proportional to the estimated target acceleration that is
  // perpendicular to the current line-of-sight, projected to remain normal
  // to the missile velocity (energy-conserving steering).
  const double apnGain = std::clamp(apnTargetAccelGain, 0.0, 24.0);
  if (apnGain > 0.0 && math::isFinite(targetAccelKmS2)) {
    const math::Vec3d losHat = R / r;
    // Remove LOS-parallel component (does not affect cross-range miss).
    math::Vec3d aPerp = targetAccelKmS2 - losHat * math::dot(targetAccelKmS2, losHat);
    // Keep steering energy-conserving: remove any component along the missile velocity.
    aPerp -= vHat * math::dot(aPerp, vHat);
    const double maxA = std::max(0.0, apnMaxTargetAccelKmS2);
    if (maxA > 0.0) aPerp = math::clampMagnitude(aPerp, maxA);
    aCmd += aPerp * apnGain;
  }

  // Convert to heading rate (rad/s): dvHat/dt = aCmd / |Vm|.
  math::Vec3d dvHatDt = aCmd / sp;

  // Clamp to the missile's steering limit.
  const double turn = dvHatDt.length();
  if (turn > maxTurnRate) {
    dvHatDt *= (maxTurnRate / turn);
  }

  math::Vec3d vNew = vHat + dvHatDt * dtSim;
  if (vNew.lengthSq() < 1e-12) return vHat;
  return vNew.normalized();
}

static bool movingSphereSegmentIntersectionEnter01(const math::Vec3d& segA,
                                                   const math::Vec3d& segB,
                                                   const SphereTarget& target,
                                                   double hitRadiusKm,
                                                   double dtSim,
                                                   double& outTEnter01) {
  outTEnter01 = 0.0;

  // Linear relative motion: subtract target's swept center from the projectile segment.
  const math::Vec3d c0 = target.centerKm;
  const math::Vec3d c1 = target.centerKm + target.velKmS * dtSim;

  const math::Vec3d relA = segA - c0;
  const math::Vec3d relB = segB - c1;

  return math::segmentSphereIntersectionEnterT(relA,
                                               relB,
                                               /*center=*/math::Vec3d{0, 0, 0},
                                               std::max(0.0, hitRadiusKm),
                                               outTEnter01);
}

void stepProjectiles(std::vector<Projectile>& projectiles,
                     double dtSim,
                     const SphereTarget* targets,
                     std::size_t targetCount,
                     std::vector<ProjectileHit>& outHits) {
  if (dtSim <= 0.0) return;
  if (projectiles.empty()) return;

  outHits.clear();

  for (auto& p : projectiles) {
    if (p.ttlSimSec <= 0.0) continue;

    const math::Vec3d a = p.posKm;
    const math::Vec3d b = p.posKm + p.velKmS * dtSim;

    p.prevKm = a;
    p.posKm = b;
    p.ttlSimSec -= dtSim;

    if (!targets || targetCount == 0) continue;

    // Decide which target kinds this projectile can hit.
    const bool allowHitPlayer = !p.fromPlayer;
    const bool allowHitShips = true;
    const bool allowHitAsteroids = p.fromPlayer; // fun: player slugs can smack rocks

    // Pick the earliest intersection along the segment (prevents order-dependent hits).
    double bestT = 2.0; // [0,1]
    const SphereTarget* bestTarget = nullptr;

    for (std::size_t ti = 0; ti < targetCount; ++ti) {
      const SphereTarget& t = targets[ti];

      if (!allowHitPlayer && t.kind == CombatTargetKind::Player) continue;
      if (!allowHitShips && t.kind == CombatTargetKind::Ship) continue;
      if (!allowHitAsteroids && t.kind == CombatTargetKind::Asteroid) continue;

      // Avoid immediate self-hits when the shooter has a real id.
      if (p.shooterId != 0 && t.id != 0 && t.id == p.shooterId) continue;

      const double hitRadiusKm = std::max(0.0, t.radiusKm) + std::max(0.0, p.radiusKm);

      double tEnter01 = 0.0;
      if (!movingSphereSegmentIntersectionEnter01(a, b, t, hitRadiusKm, dtSim, tEnter01)) continue;

      if (tEnter01 < bestT) {
        bestT = tEnter01;
        bestTarget = &t;

        // Can't do better than an immediate hit.
        if (bestT <= 0.0) break;
      }
    }

    if (!bestTarget) continue;

    const math::Vec3d hitPoint = a + (b - a) * bestT;

    ProjectileHit h{};
    h.kind = bestTarget->kind;
    h.targetIndex = bestTarget->index;
    h.targetId = bestTarget->id;
    h.pointKm = hitPoint;
    h.dmg = p.dmg;
    h.fromPlayer = p.fromPlayer;
    h.shooterId = p.shooterId;
    outHits.push_back(h);

    // Consume the projectile and keep its end position consistent with the hit.
    p.posKm = hitPoint;
    p.ttlSimSec = 0.0;
  }
}

static void stepMissilesImpl(std::vector<Missile>& missiles,
                             std::vector<Projectile>* interceptProjectiles,
                             double dtSim,
                             const SphereTarget* targets,
                             std::size_t targetCount,
                             std::vector<MissileDetonation>& outDetonations,
                             std::vector<MissileHit>& outHits) {
  if (dtSim <= 0.0) return;
  if (missiles.empty()) return;

  outDetonations.clear();
  outHits.clear();

  // Snapshot missile state at the start of the step for deterministic
  // swarm/formation steering. This avoids order dependence inside the
  // per-missile integration loop.
  std::vector<MissileSwarmSnapshot> swarmSnap;
  swarmSnap.reserve(missiles.size());
  for (const auto& m : missiles) {
    MissileSwarmSnapshot s{};
    s.alive = (m.ttlSimSec > 0.0);
    s.fromPlayer = m.fromPlayer;
    s.shooterId = m.shooterId;
    s.hasTarget = m.hasTarget;
    s.targetKind = m.targetKind;
    s.targetId = m.targetId;
    s.posKm = m.posKm;
    math::Vec3d v = m.velKmS;
    if (v.lengthSq() < 1e-12) v = {0, 0, 1};
    s.velDir = v.normalized();
    swarmSnap.push_back(s);
  }

  for (std::size_t mi = 0; mi < missiles.size(); ++mi) {
    auto& m = missiles[mi];
    if (m.ttlSimSec <= 0.0) continue;

    const double ageStart = std::max(0.0, m.ageSimSec);
    const double activateAfter = std::max(0.0, m.seekerActivationSimSec);
    const bool seekerActive = (activateAfter <= 0.0) || (ageStart >= activateAfter);

    // Optional: discrete seeker measurement updates.
    //
    // When enabled (period > 0), we only treat the active seeker as having a fresh measurement
    // on specific "update ticks". Between ticks, guidance falls back to inertial target memory
    // and decoy/reacquire decisions are deferred.
    bool seekerUpdateThisStep = true;
    if (seekerActive) {
      const double period = std::max(0.0, m.seekerUpdatePeriodSimSec);
      if (period > 1e-9) {
        // Time since seeker activation (clamped to >= 0).
        const double t0 = std::max(0.0, ageStart - activateAfter);
        const double t1 = std::max(0.0, ageStart + dtSim - activateAfter);

        // Measurement ticks occur at t = n*period (including t=0). Add a tiny epsilon so
        // boundaries like t==1.0 land on the expected tick under floating-point accumulation.
        const double eps = 1e-9;
        const long long k0 = (long long)std::floor((t0 + eps) / period);
        const long long k1 = (long long)std::floor((t1 + eps) / period);

        seekerUpdateThisStep = (t0 <= 1e-12) || (k1 != k0);
      }
    }

    // Optional: discrete midcourse datalink target updates.
    //
    // When enabled (period > 0) and the seeker is not yet active, we only ingest
    // fresh target updates on specific "datalink ticks". Between ticks, guidance
    // falls back to inertial target memory.
    bool datalinkUpdateThisStep = true;
    if (!seekerActive) {
      const double period = std::max(0.0, m.datalinkUpdatePeriodSimSec);
      if (period > 1e-9) {
        const double latency = std::max(0.0, m.datalinkLatencySimSec);
        const double t0 = std::max(0.0, ageStart);
        const double t1 = std::max(0.0, ageStart + dtSim);

        // Datalink update arrival times are modeled as:
        //   t = latency + n*period
        // We detect if the current integration step crosses one or more ticks.
        const double eps = 1e-9;
        const long long k0 = (long long)std::floor(((t0 - latency) + eps) / period);
        const long long k1 = (long long)std::floor(((t1 - latency) + eps) / period);

        // Always allow an initial "launch fix" at t=0 even if latency > 0.
        datalinkUpdateThisStep = (t0 <= 1e-12) || (k1 != k0);
      }
    }

    const math::Vec3d a = m.posKm;

    // --- Guidance ---
    math::Vec3d v = m.velKmS;
    double speed = v.length();
    if (speed < 1e-6) {
      speed = 1.0;
      v = {0, 0, 1};
    }

    math::Vec3d desiredDir = v.normalized();
    // When missiles have a motor model, their speed can increase during dtSim. When we also
    // have a lateral-acceleration ("G") cap, we must compute the effective turn rate using a
    // conservative (high) speed estimate so we don't violate the cap after accelerating.
    double speedForTurn = speed;
    {
      const double accel = std::max(0.0, m.thrustAccelKmS2);
      const double maxSpeed = std::max(0.0, m.maxSpeedKmS);
      const double burnRem = std::max(0.0, m.motorBurnRemainingSimSec);

      if (accel > 0.0 && burnRem > 0.0 && maxSpeed > speedForTurn + 1e-12) {
        const double burnDt = std::min(dtSim, burnRem);
        speedForTurn = std::min(maxSpeed, speedForTurn + accel * burnDt);
      }
    }

    const double turnRateRadS = effectiveTurnRateRadS(m, speedForTurn);
    const SphereTarget* guide = nullptr;

    // When true, the active seeker is maintaining track via Home-On-Jam (HOJ)
    // rather than a clean doppler return (i.e. inside the notch).
    bool lockViaJam = false;

    // Optional: decoy commitment ("minimum dwell") timer.
    if (m.decoyCommitRemainingSimSec > 0.0) {
      m.decoyCommitRemainingSimSec = std::max(0.0, m.decoyCommitRemainingSimSec - dtSim);
    }

    // --- Target lock, seeker FOV, and target memory ---
    const math::Vec3d fwd = v.normalized();
    const double seekerFovCos = std::clamp(m.seekerFovCos, -1.0, 1.0);

    // Age / expire the last-known target memory.
    {
      double memMax = std::max(0.0, m.targetMemorySimSec);

      // If the seeker is not yet active and midcourse datalink updates are discrete,
      // ensure target memory persists long enough to bridge the next expected update.
      //
      // This avoids the surprising failure mode where datalinkUpdatePeriodSimSec>0 but
      // targetMemorySimSec==0, causing the missile to "forget" its last update before
      // the next tick.
      if (!seekerActive) {
        const double dlPeriod = std::max(0.0, m.datalinkUpdatePeriodSimSec);
        if (dlPeriod > 1e-9) {
          const double dlLatency = std::max(0.0, m.datalinkLatencySimSec);
          const double implied = dlPeriod + dlLatency + 2.0 * dtSim;
          memMax = std::max(memMax, implied);
        }
      }

      if (m.hasLastKnownTarget) {
        m.lastKnownTargetAgeSimSec = std::max(0.0, m.lastKnownTargetAgeSimSec + dtSim);
        if (memMax <= 0.0 || m.lastKnownTargetAgeSimSec > memMax) {
          m.hasLastKnownTarget = false;
        }
      }
    }

    // Optional APN (augmented proportional navigation) state updates.
    //
    // We keep a finite-difference estimate of the target acceleration from
    // direct velocity samples, and decay it toward zero when we don't have
    // fresh measurements.
    if (m.apnTargetAccelGain > 0.0) {
      if (m.hasTargetVelSample) {
        m.targetVelSampleAgeSimSec = std::min(5.0, std::max(0.0, m.targetVelSampleAgeSimSec + dtSim));
      }
      const double hl = std::max(0.0, m.apnAccelHalfLifeSimSec);
      if (hl > 0.0) {
        const double decay = math::halfLifeDecayFactor(dtSim, hl);
        m.targetAccelEstKmS2 *= decay;
      } else {
        // Without smoothing, treat the estimate as per-step only.
        m.targetAccelEstKmS2 = {0, 0, 0};
      }
    }

    const SphereTarget* lockedTarget = nullptr;

    // Base target guide: either the currently trackable locked target, or a synthetic
    // memory target if tracking is temporarily lost.
    SphereTarget memoryTarget{};
    const SphereTarget* baseTarget = nullptr;

    // "Lock" score used to compare against decoys (inverse-square falloff).
    double targetScore = 0.0;

    if (m.hasTarget && targets && targetCount > 0) {
      for (std::size_t i = 0; i < targetCount; ++i) {
        const SphereTarget& t = targets[i];
        if (t.kind != m.targetKind) continue;
        if (t.id != m.targetId) continue;
        lockedTarget = &t;
        break;
      }
    }

    // Track the locked target.
    //
    // When the seeker is inactive (activation delay), the missile can optionally receive
    // midcourse target updates via a simple "datalink" model (range + optional LOS).
    // If datalinkRangeKm is 0, behavior matches the legacy perfect-update midcourse.
    // Once active, we enforce the seeker's FOV and allow decoys to override.
    if (lockedTarget) {
      const math::Vec3d toTgt = lockedTarget->centerKm - m.posKm;
      const double distSq = toTgt.lengthSq();
      if (distSq > 1e-12) {
        const double dist = std::sqrt(distSq);
        const math::Vec3d toDir = toTgt / dist;
        const double cosAng = math::dot(fwd, toDir);

        bool trackable = false;

        if (seekerActive) {
          trackable = (cosAng >= seekerFovCos);

          // If the seeker only updates intermittently, treat the target as untrackable between
          // update ticks so we rely on inertial target memory for guidance.
          if (trackable && !seekerUpdateThisStep) {
            trackable = false;
          }
        } else {
          // Midcourse (pre-activation): either legacy "perfect" updates (datalinkRangeKm==0),
          // or a range/LOS gated datalink from the launching platform.
          const double dlRangeKm = std::max(0.0, m.datalinkRangeKm);
          if (dlRangeKm <= 0.0) {
            trackable = true;
          } else if (targets && targetCount > 0) {
            const SphereTarget* shooter = nullptr;
            for (std::size_t si = 0; si < targetCount; ++si) {
              const SphereTarget& s = targets[si];
              if (s.id != m.shooterId) continue;
              if (s.kind != CombatTargetKind::Ship && s.kind != CombatTargetKind::Player) continue;
              shooter = &s;
              break;
            }

            if (shooter) {
              const double dlRangeSq = dlRangeKm * dlRangeKm;
              const double distLinkSq = (lockedTarget->centerKm - shooter->centerKm).lengthSq();
              if (distLinkSq <= dlRangeSq + 1.0e-9) {
                trackable = true;

                if (m.datalinkRequireLineOfSight) {
                  const double pad = std::max(0.0, m.datalinkOcclusionPadKm);
                  if (occludedByAsteroidSpheres(shooter->centerKm,
                                               lockedTarget->centerKm,
                                               targets,
                                               targetCount,
                                               pad)) {
                    trackable = false;
                  }
                }
              }
            }
          }

          // If midcourse datalink updates are discrete, treat the target as untrackable between
          // update ticks so guidance falls back to inertial memory.
          if (trackable) {
            const double dlPeriod = std::max(0.0, m.datalinkUpdatePeriodSimSec);
            if (dlPeriod > 1e-9 && !datalinkUpdateThisStep) {
              trackable = false;
            }
          }
        }

        // Optional LOS + (radar) doppler notch filters apply only during the active seeker phase.
        if (trackable && seekerActive) {
          if (m.requireLineOfSight) {
            const double pad = std::max(0.0, m.lineOfSightOcclusionPadKm);
            if (occludedByAsteroidSpheres(m.posKm, lockedTarget->centerKm, targets, targetCount, pad)) {
              trackable = false;
            }
          }

          if (trackable && m.seeker == MissileSeekerType::Radar) {
            const double notch = std::max(0.0, m.radarDopplerNotchKmS);
            if (notch > 0.0) {
              // Radial velocity along the line-of-sight. Near-zero magnitude implies a "beam"/notch.
              const double vrKmS = math::dot(lockedTarget->velKmS - v, toDir);
              if (std::fabs(vrKmS) < notch) {
                // Inside the doppler notch.
                //
                // If HOJ is enabled and the target is actively jamming, allow the
                // missile to maintain a coarse track by homing on the jammer emission.
                bool hoj = false;
                if (m.homeOnJam) {
                  const double minJam = std::max(0.0, m.homeOnJamMinJammerPower);
                  const double jp = std::max(0.0, lockedTarget->jammerPower);
                  if (jp > 0.0 && jp >= minJam) {
                    hoj = true;

                    // Optional: gate HOJ on *received* jamming level when configured.
                    const double minJ01 = std::clamp(m.homeOnJamMinJamming01, 0.0, 1.0);
                    if (minJ01 > 0.0) {
                      const double j01 = receivedJamming01(m, dist, jp);
                      if (j01 < minJ01) hoj = false;
                    }
                  }
                }

                if (!hoj) {
                  trackable = false;
                } else {
                  lockViaJam = true;
                }
              }
            }
          }
        }

        if (trackable) {
          // Seed / maintain an internal seeker pointing direction.
          //
          // During the pre-activation phase (midcourse datalink), keep this aligned to the locked
          // target so the active seeker can "come up" already pointed near the correct track.
          if (!m.hasSeekerDir) {
            m.hasSeekerDir = true;
            m.seekerDirWorld = toDir;
          } else if (!seekerActive) {
            m.seekerDirWorld = toDir;
          }

          // Update target memory even if a decoy later overrides guidance.
          m.hasLastKnownTarget = true;
          m.lastKnownTargetPosKm = lockedTarget->centerKm;
          m.lastKnownTargetVelKmS = lockedTarget->velKmS;
          m.lastKnownTargetAccelKmS2 = lockedTarget->accelKmS2;
          m.lastKnownTargetHeatSignature = lockedTarget->heatSignature;
          m.lastKnownTargetRadarSignature = lockedTarget->radarSignature;
          m.lastKnownTargetAgeSimSec = 0.0;

          // APN target acceleration estimation (optional):
          // finite-difference acceleration from direct velocity samples.
          if (m.apnTargetAccelGain > 0.0) {
            if (!m.hasTargetVelSample) {
              m.hasTargetVelSample = true;
              m.targetVelSampleKmS = lockedTarget->velKmS;
              m.targetVelSampleAgeSimSec = 0.0;
              m.targetAccelEstKmS2 = {0, 0, 0};
            } else {
              const double dtSample = std::max(1e-6, m.targetVelSampleAgeSimSec);
              math::Vec3d aRaw = (lockedTarget->velKmS - m.targetVelSampleKmS) / dtSample;
              const double maxA = std::max(0.0, m.apnMaxTargetAccelKmS2);
              if (maxA > 0.0) aRaw = math::clampMagnitude(aRaw, maxA);
              const double hl = std::max(0.0, m.apnAccelHalfLifeSimSec);
              if (hl > 0.0) {
                const double alpha = math::halfLifeAlpha(dtSim, hl);
                m.targetAccelEstKmS2 = math::lerp(m.targetAccelEstKmS2, aRaw, alpha);
              } else {
                m.targetAccelEstKmS2 = aRaw;
              }
              m.targetVelSampleKmS = lockedTarget->velKmS;
              m.targetVelSampleAgeSimSec = 0.0;
            }
          }

          baseTarget = lockedTarget;
          double sig = 1.0;
          if (m.seeker == MissileSeekerType::Heat) {
            sig = effectiveHeatSignature(m, lockedTarget->heatSignature, lockedTarget->velKmS, toDir);
          } else if (m.seeker == MissileSeekerType::Radar) {
            sig = effectiveRadarSignature(lockedTarget->radarSignature);
          }
          targetScore = sig / (distSq + 1.0e-9);
        }
      }
    }

    // If the target is temporarily untrackable (e.g. leaves the seeker cone), fall back to
    // a simple inertial-memory prediction for a short window.
    if (!baseTarget && m.hasLastKnownTarget) {
      memoryTarget.kind = m.targetKind;
      memoryTarget.index = 0;
      memoryTarget.id = m.targetId;
      const double tAge = std::max(0.0, m.lastKnownTargetAgeSimSec);
      const math::Vec3d acc = m.lastKnownTargetAccelKmS2;
      memoryTarget.centerKm = m.lastKnownTargetPosKm +
                            m.lastKnownTargetVelKmS * tAge +
                            acc * (0.5 * tAge * tAge);
      memoryTarget.velKmS = m.lastKnownTargetVelKmS + acc * tAge;
      memoryTarget.accelKmS2 = acc;
      memoryTarget.heatSignature = m.lastKnownTargetHeatSignature;
      memoryTarget.radarSignature = m.lastKnownTargetRadarSignature;
      memoryTarget.radiusKm = 0.0;

      baseTarget = &memoryTarget;

      const math::Vec3d toMem = memoryTarget.centerKm - m.posKm;
      double sig = 1.0;
      if (m.seeker == MissileSeekerType::Heat) {
        const math::Vec3d toMemDir = math::safeNormalized(toMem, {0, 0, 1}, 1e-12);
        sig = effectiveHeatSignature(m, m.lastKnownTargetHeatSignature, memoryTarget.velKmS, toMemDir);
      } else if (m.seeker == MissileSeekerType::Radar) {
        sig = effectiveRadarSignature(m.lastKnownTargetRadarSignature);
      }
      targetScore = sig / (toMem.lengthSq() + 1.0e-9);
    }

    // --- Seeker track quality (optional) ---
    //
    // Track quality rises while we have a direct measurement of the locked target and
    // decays while guiding on memory or without a track. It is used to modulate decoy
    // resistance so countermeasures become more effective after a lock break.
    if (m.enableTrackQuality && seekerActive) {
      m.trackQuality = std::clamp(m.trackQuality, 0.0, 1.0);

      const double riseHl = std::max(0.0, m.trackQualityRiseHalfLifeSimSec);
      const double fallHl = std::max(0.0, m.trackQualityFallHalfLifeSimSec);

      // Measurement quality in [0,1].
      //
      // With seeker slew enabled, lock quality can drop even while the target remains inside the
      // FOV: high line-of-sight rates (close merges / sharp maneuvers) can outpace the seeker.
      // This reduces decoy resistance and makes evasive flying + countermeasures matter.
      double measQ = 0.0;
      const bool directTrack = (lockedTarget && baseTarget == lockedTarget);
      if (directTrack) {
        // Default (no slew): behave like a perfect lock.
        measQ = 1.0;

        if (m.seekerSlewRateRadS > 0.0) {
          const math::Vec3d to = lockedTarget->centerKm - m.posKm;
          if (to.lengthSq() > 1e-12) {
            const math::Vec3d toDir = to.normalized();

            // If the state wasn't seeded (e.g. missiles spawned mid-flight), assume it starts aligned.
            if (!m.hasSeekerDir) {
              m.hasSeekerDir = true;
              m.seekerDirWorld = toDir;
            }

            math::Vec3d sDir = m.seekerDirWorld;
            if (sDir.lengthSq() < 1e-12) sDir = toDir;
            else sDir = sDir.normalized();

            const double align = std::clamp(math::dot(sDir, toDir), -1.0, 1.0);
            const double denom = std::max(1e-6, 1.0 - seekerFovCos);
            measQ = std::clamp((align - seekerFovCos) / denom, 0.0, 1.0);
          }
        }
        // Optional: range-based received-signal falloff.
        //
        // When enabled (seekerSignalHalfRangeKm > 0), this reduces track quality at long range
        // (or against low-signature targets) without adding randomness. It also naturally
        // increases susceptibility to countermeasures when the seeker is "struggling" to see the
        // target.
        if (m.seekerSignalHalfRangeKm > 0.0) {
          const math::Vec3d to = lockedTarget->centerKm - m.posKm;
          const double dist = std::sqrt(std::max(0.0, to.lengthSq()));
          const math::Vec3d toDir = math::safeNormalized(to, math::Vec3d{0, 0, 1}, 1e-12);

          double sig = 1.0;
          if (m.seeker == MissileSeekerType::Heat) {
            sig = effectiveHeatSignature(m, lockedTarget->heatSignature, lockedTarget->velKmS, toDir);
          } else if (m.seeker == MissileSeekerType::Radar) {
            sig = effectiveRadarSignature(lockedTarget->radarSignature);
          }

          const double s01 = receivedSeekerSignal01(m, dist, sig);
          measQ *= s01;
        }

      }

      // HOJ lock has weaker angular accuracy than a clean radar track.
      if (directTrack && lockViaJam && m.seeker == MissileSeekerType::Radar && m.homeOnJam) {
        const double cap = std::clamp(m.homeOnJamTrackQualityCap, 0.0, 1.0);
        measQ = std::min(measQ, cap);
      }

      // Optional: noise jamming can degrade radar track quality even outside the notch.
      //
      // This is a simple deterministic suppression model controlled by the missile:
      //   measQ *= 1 / (1 + gain * jamming01)
      // where jamming01 is an (optional) distance-scaled received jamming level.
      if (directTrack && m.seeker == MissileSeekerType::Radar) {
        const double gain = std::max(0.0, m.radarJammingTrackSuppressionGain);
        if (gain > 0.0 && lockedTarget) {
          const double jp = std::max(0.0, lockedTarget->jammerPower);
          if (jp > 0.0) {
            const math::Vec3d to = lockedTarget->centerKm - m.posKm;
            const double dist = std::sqrt(std::max(0.0, to.lengthSq()));
            const double j01 = receivedJamming01(m, dist, jp);
            const double denom = 1.0 + gain * j01;
            if (std::isfinite(denom) && denom > 0.0) {
              measQ = std::clamp(measQ / denom, 0.0, 1.0);
            }
          }
        }
      }



      // Optional: countermeasure clutter can suppress direct track quality even if the
      // missile does not fully switch lock to a decoy.
      //
      // This models seeker confusion/sidelobe clutter in a deterministic way by comparing the
      // strongest in-FOV decoy score to the current target score.
      if (directTrack && seekerUpdateThisStep) {
        const double gain = std::max(0.0, m.decoyClutterTrackSuppressionGain);
        if (gain > 0.0 && targetScore > 0.0 && targets && targetCount > 0 && lockedTarget) {
          // Gate/reference values based on the currently tracked target.
          const math::Vec3d toTgt = lockedTarget->centerKm - m.posKm;
          const double tgtR2 = toTgt.lengthSq();
          if (tgtR2 > 1e-12) {
            const double tgtR = std::sqrt(tgtR2);
            const math::Vec3d toTargetDir = toTgt / tgtR;
            const double targetVrKmS = math::dot(lockedTarget->velKmS - m.velKmS, toTargetDir);

            const double angleGateCos = std::clamp(m.decoyAngleGateCos, -1.0, 1.0);
            const double rangeGateKm = std::max(0.0, m.decoyRangeGateKm);
            const double dopplerGateKmS = std::max(0.0, m.decoyDopplerGateKmS);

            // Optional close-range burn-through (radar seekers).
            double burnThroughFactor = 1.0;
            if (m.seeker == MissileSeekerType::Radar) {
              const double btRange = std::max(0.0, m.decoyBurnThroughRangeKm);
              if (btRange > 1e-9) {
                const double minF = std::clamp(m.decoyBurnThroughMinFactor, 0.0, 1.0);
                const double t = std::clamp(tgtR / btRange, 0.0, 1.0);
                burnThroughFactor = minF * (1.0 - t) + t;
              }
            }

            // Current seeker pointing direction used for slew-aware clutter scoring.
            math::Vec3d seekerDir = fwd;
            if (m.hasSeekerDir) seekerDir = m.seekerDirWorld;
            seekerDir = math::safeNormalized(seekerDir, fwd, 1e-12);

            double bestDecoyScore = 0.0;

            for (std::size_t i = 0; i < targetCount; ++i) {
              const SphereTarget& t = targets[i];
              if (t.kind != CombatTargetKind::Decoy) continue;

              const double strength = (m.seeker == MissileSeekerType::Radar) ? t.decoyRadar : t.decoyHeat;
              if (strength <= 0.0) continue;

              const math::Vec3d to = t.centerKm - m.posKm;
              const double r2 = to.lengthSq();
              if (r2 < 1e-12) continue;

              const double r = std::sqrt(r2);
              const math::Vec3d toDir = to / r;

              const double cosAng = math::dot(fwd, toDir);
              if (cosAng < seekerFovCos) continue;

              if (m.requireLineOfSight) {
                const double pad = std::max(0.0, m.lineOfSightOcclusionPadKm);
                if (occludedByAsteroidSpheres(m.posKm, t.centerKm, targets, targetCount, pad)) {
                  continue;
                }
              }

              double score = (strength * std::max(0.0, cosAng)) / (r2 + 1.0e-9);
              score *= burnThroughFactor;
              score *= seekerSlewPenaltyFactor(m, dtSim, seekerDir, toDir);

              // Apply optional gates relative to the locked track.
              if (angleGateCos > -0.5) {
                const double sepCos = math::dot(toTargetDir, toDir);
                if (sepCos < angleGateCos) continue;
              }

              if (rangeGateKm > 0.0) {
                const double sep = std::fabs(r - tgtR);
                const double k = 1.0 - sep / rangeGateKm;
                if (k <= 0.0) continue;
                score *= std::clamp(k, 0.0, 1.0);
              }

              if (dopplerGateKmS > 0.0) {
                const double vr = math::dot(t.velKmS - m.velKmS, toDir);
                const double dv = std::fabs(vr - targetVrKmS);
                const double k = 1.0 - dv / dopplerGateKmS;
                if (k <= 0.0) continue;
                score *= std::clamp(k, 0.0, 1.0);
              }

              if (score > bestDecoyScore) bestDecoyScore = score;
            }

            if (bestDecoyScore > 0.0) {
              const double clutter = bestDecoyScore / std::max(1.0e-12, targetScore);
              const double denom = 1.0 + gain * std::max(0.0, clutter);
              if (std::isfinite(denom) && denom > 0.0) {
                measQ = std::clamp(measQ / denom, 0.0, 1.0);
              }
            }
          }
        }
      }
      const bool rising = (measQ >= m.trackQuality);
      const double hl = rising ? riseHl : fallHl;
      const double alpha = math::halfLifeAlpha(dtSim, hl);
      m.trackQuality = std::clamp(math::lerp(m.trackQuality, measQ, alpha), 0.0, 1.0);
    } else {
      // Clamp defensively so callers can safely read the value.
      m.trackQuality = std::clamp(m.trackQuality, 0.0, 1.0);
    }

    // Effective decoy resistance can be reduced when track quality is poor.
    double effectiveDecoyResistance = std::max(0.0, m.decoyResistance);
    if (m.enableTrackQuality && seekerActive) {
      const double floor = std::clamp(m.trackQualityResistFloor, 0.0, 1.0);
      effectiveDecoyResistance *= math::lerp(floor, 1.0, m.trackQuality);
    }

    // --- Countermeasure / decoy attraction ---
    const SphereTarget* candidateGuide = baseTarget;

    const SphereTarget* bestDecoy = nullptr;
    double bestDecoyScore = 0.0;

    // Only consider decoys when we have a meaningful locked target score (real or memory)
    // *and* the seeker's active phase has begun.
    if (seekerActive && seekerUpdateThisStep && targetScore > 0.0 && targets && targetCount > 0) {

      // Optional range/angle/doppler discrimination gates around the currently tracked target.
      // These help (radar) seekers ignore decoys that are obviously not part of the true track.
      const double angleGateCos = std::clamp(m.decoyAngleGateCos, -1.0, 1.0);
      const double rangeGateKm = std::max(0.0, m.decoyRangeGateKm);
      const double dopplerGateKmS = std::max(0.0, m.decoyDopplerGateKmS);

      bool haveGateRef = false;
      math::Vec3d toTargetDir{0, 0, 1};
      double targetRangeKm = 0.0;
      double targetVrKmS = 0.0;

      if (baseTarget) {
        const math::Vec3d toTgt = baseTarget->centerKm - m.posKm;
        const double r2 = toTgt.lengthSq();
        if (r2 > 1e-12) {
          targetRangeKm = std::sqrt(r2);
          toTargetDir = toTgt / targetRangeKm;

          // Positive vr means the track is opening (range increasing).
          targetVrKmS = math::dot(baseTarget->velKmS - m.velKmS, toTargetDir);
          haveGateRef = true;
        }
      }

      // Optional close-range burn-through: reduce decoy attraction when the
      // missile is near its current target track.
      double burnThroughFactor = 1.0;
      if (m.seeker == MissileSeekerType::Radar) {
        const double btRange = std::max(0.0, m.decoyBurnThroughRangeKm);
        if (btRange > 1e-9 && haveGateRef) {
          const double minF = std::clamp(m.decoyBurnThroughMinFactor, 0.0, 1.0);
          const double t = std::clamp(targetRangeKm / btRange, 0.0, 1.0);
          burnThroughFactor = minF * (1.0 - t) + t;
        }
      }

      // Current seeker pointing direction used for slew-aware decoy scoring.
      math::Vec3d seekerDir = fwd;
      if (m.hasSeekerDir) seekerDir = m.seekerDirWorld;
      seekerDir = math::safeNormalized(seekerDir, fwd, 1e-12);

      for (std::size_t i = 0; i < targetCount; ++i) {
        const SphereTarget& t = targets[i];
        if (t.kind != CombatTargetKind::Decoy) continue;

        const double strength = (m.seeker == MissileSeekerType::Radar) ? t.decoyRadar : t.decoyHeat;
        if (strength <= 0.0) continue;

        const math::Vec3d to = t.centerKm - m.posKm;
        const double distSq = to.lengthSq();
        if (distSq < 1e-12) continue;

        const double dist = std::sqrt(distSq);
        const math::Vec3d toDir = to / dist;

        const double cosAng = math::dot(fwd, toDir);
        if (cosAng < seekerFovCos) continue;

        if (m.requireLineOfSight) {
          const double pad = std::max(0.0, m.lineOfSightOcclusionPadKm);
          if (occludedByAsteroidSpheres(m.posKm, t.centerKm, targets, targetCount, pad)) {
            continue;
          }
        }

        double score = (strength * std::max(0.0, cosAng)) / (distSq + 1.0e-9);

        // Burn-through is applied after basic geometry and FOV checks.
        score *= burnThroughFactor;

        // With seeker gimbal slew enabled, large off-axis jumps are harder: the seeker
        // may not be able to re-point to a decoy instantly.
        score *= seekerSlewPenaltyFactor(m, dtSim, seekerDir, toDir);

        // Apply optional gates relative to the locked track (if available).
        if (haveGateRef) {
          if (angleGateCos > -0.5) {
            const double sepCos = math::dot(toTargetDir, toDir);
            if (sepCos < angleGateCos) continue;
          }

          if (rangeGateKm > 0.0) {
            const double sep = std::fabs(dist - targetRangeKm);
            const double k = 1.0 - sep / rangeGateKm;
            if (k <= 0.0) continue;
            score *= std::clamp(k, 0.0, 1.0);
          }

          if (dopplerGateKmS > 0.0) {
            const double vr = math::dot(t.velKmS - m.velKmS, toDir);
            const double dv = std::fabs(vr - targetVrKmS);
            const double k = 1.0 - dv / dopplerGateKmS;
            if (k <= 0.0) continue;
            score *= std::clamp(k, 0.0, 1.0);
          }
        }

        if (score > bestDecoyScore) {
          bestDecoyScore = score;
          bestDecoy = &t;
        }
      }

      const double resist = effectiveDecoyResistance;
      if (bestDecoy && bestDecoyScore > targetScore * resist) {
        candidateGuide = bestDecoy;
      }
    }

    // --- Autonomous target acquisition (optional) ---
    //
    // If we don't have a trackable target (real or memory) and the seeker is active, allow
    // the missile to acquire a new Ship/Player target within a limited range.
    if (!candidateGuide && seekerActive && seekerUpdateThisStep) {
      const double acquireRangeKm = std::max(0.0, m.autoAcquireRangeKm);
      if (acquireRangeKm > 0.0 && targets && targetCount > 0) {
        const bool allowPlayer = !m.fromPlayer;

        const SphereTarget* best = nullptr;
        double bestScore = 0.0;

        // Current seeker pointing direction used for slew-aware acquisition scoring.
        math::Vec3d seekerDir = fwd;
        if (m.hasSeekerDir) seekerDir = m.seekerDirWorld;
        seekerDir = math::safeNormalized(seekerDir, fwd, 1e-12);

        for (std::size_t i = 0; i < targetCount; ++i) {
          const SphereTarget& t = targets[i];
          if (t.kind != CombatTargetKind::Ship && t.kind != CombatTargetKind::Player) continue;
          if (!allowPlayer && t.kind == CombatTargetKind::Player) continue;
          if (m.shooterId != 0 && t.id != 0 && t.id == m.shooterId) continue;

          const math::Vec3d to = t.centerKm - m.posKm;
          const double distSq = to.lengthSq();
          if (distSq < 1e-12) continue;

          const double dist = std::sqrt(distSq);
          if (dist > acquireRangeKm) continue;

          const math::Vec3d toDir = to / dist;
          const double cosAng = math::dot(fwd, toDir);

          double req = seekerFovCos;
          if (t.minAimCos > -0.5) req = std::max(req, t.minAimCos);
          if (cosAng < req) continue;

          if (m.requireLineOfSight) {
            const double pad = std::max(0.0, m.lineOfSightOcclusionPadKm);
            if (occludedByAsteroidSpheres(m.posKm, t.centerKm, targets, targetCount, pad)) {
              continue;
            }
          }

          if (m.seeker == MissileSeekerType::Radar) {
            const double notch = std::max(0.0, m.radarDopplerNotchKmS);
            if (notch > 0.0) {
              const double vrKmS = math::dot(t.velKmS - v, toDir);
              if (std::fabs(vrKmS) < notch) {
                bool hoj = false;
                if (m.homeOnJam) {
                  const double minJam = std::max(0.0, m.homeOnJamMinJammerPower);
                  const double jp = std::max(0.0, t.jammerPower);
                  if (jp > 0.0 && jp >= minJam) {
                    hoj = true;

                    const double minJ01 = std::clamp(m.homeOnJamMinJamming01, 0.0, 1.0);
                    if (minJ01 > 0.0) {
                      const double j01 = receivedJamming01(m, dist, jp);
                      if (j01 < minJ01) hoj = false;
                    }
                  }
                }
                if (!hoj) {
                  continue;
                }
              }
            }
          }

          // Prefer targets close to boresight and nearby.
          // Heat seekers also bias toward hotter targets.
          double sig = 1.0;
          if (m.seeker == MissileSeekerType::Heat) {
            sig = effectiveHeatSignature(m, t.heatSignature, t.velKmS, toDir);
          } else if (m.seeker == MissileSeekerType::Radar) {
            sig = effectiveRadarSignature(t.radarSignature);
          }
          double score = sig * std::max(0.0, cosAng) / (distSq + 1.0e-9);

          // With seeker gimbal slew enabled, prefer targets that don't require a huge
          // instantaneous re-point.
          score *= seekerSlewPenaltyFactor(m, dtSim, seekerDir, toDir);

          // If HOJ is enabled, bias reacquisition toward active jammers.
          if (m.homeOnJam && m.seeker == MissileSeekerType::Radar) {
            const double minJam = std::max(0.0, m.homeOnJamMinJammerPower);
            if (t.jammerPower >= minJam && t.jammerPower > 0.0) {
              const double bias = std::max(0.0, m.homeOnJamAcquireBias);
              if (bias > 0.0) {
                score *= (1.0 + bias * t.jammerPower);
              }
            }
          }

          if (score > bestScore) {
            bestScore = score;
            best = &t;
          }
        }

        if (best) {
          // Capture the new lock.
          m.hasTarget = true;
          m.targetKind = best->kind;
          m.targetId = best->id;

          // Refresh memory and seed score for potential decoy comparison.
          m.hasLastKnownTarget = true;
          m.lastKnownTargetPosKm = best->centerKm;
          m.lastKnownTargetVelKmS = best->velKmS;
          m.lastKnownTargetAccelKmS2 = best->accelKmS2;
          m.lastKnownTargetHeatSignature = best->heatSignature;
          m.lastKnownTargetRadarSignature = best->radarSignature;
          m.lastKnownTargetAgeSimSec = 0.0;

          // Reset APN estimator state for the new lock (avoids cross-target spikes).
          m.hasTargetVelSample = false;
          m.targetVelSampleAgeSimSec = 0.0;
          m.targetAccelEstKmS2 = {0, 0, 0};

          candidateGuide = best;
          {
            const math::Vec3d to = best->centerKm - m.posKm;
            const double distSq = to.lengthSq();
            double sig = 1.0;
            if (m.seeker == MissileSeekerType::Heat) {
              const math::Vec3d toDir = math::safeNormalized(to, {0, 0, 1}, 1e-12);
              sig = effectiveHeatSignature(m, best->heatSignature, best->velKmS, toDir);
            } else if (m.seeker == MissileSeekerType::Radar) {
              sig = effectiveRadarSignature(best->radarSignature);
            }
            targetScore = sig / (distSq + 1.0e-9);
          }
        }
      }
    }

    // Optional: if we recently switched to a decoy, keep guiding to it for a minimum
    // time window (prevents rapid thrashing and makes countermeasures "commit" the missile).
    const SphereTarget* committedDecoy = nullptr;
    if (m.decoyCommitRemainingSimSec > 0.0 && m.committedDecoyId != 0 && targets && targetCount > 0) {

      for (std::size_t i = 0; i < targetCount; ++i) {
        const SphereTarget& t = targets[i];
        if (t.kind != CombatTargetKind::Decoy) continue;
        if (t.id != m.committedDecoyId) continue;

        const math::Vec3d to = t.centerKm - m.posKm;
        if (to.lengthSq() < 1e-12) {
          committedDecoy = &t;
          break;
        }

        const math::Vec3d toDir = to.normalized();
        const double cosAng = math::dot(fwd, toDir);
        if (cosAng >= seekerFovCos) {
          bool visible = true;
          if (m.requireLineOfSight) {
            const double pad = std::max(0.0, m.lineOfSightOcclusionPadKm);
            if (occludedByAsteroidSpheres(m.posKm, t.centerKm, targets, targetCount, pad)) {
              visible = false;
            }
          }
          if (visible) {
            committedDecoy = &t;
          }
        }
        break;
      }

      if (!committedDecoy) {
        // Lost the committed decoy (expired/out of FOV): clear the latch early.
        m.decoyCommitRemainingSimSec = 0.0;
        m.committedDecoyId = 0;
      }
    }

    guide = committedDecoy ? committedDecoy : candidateGuide;

    // Start a new commitment window when we first switch to a decoy.
    if (seekerActive && !committedDecoy && guide && guide->kind == CombatTargetKind::Decoy) {
      const double hold = std::max(0.0, m.decoyCommitSimSec);
      if (hold > 0.0 && guide->id != m.committedDecoyId) {
        m.committedDecoyId = guide->id;
        m.decoyCommitRemainingSimSec = hold;
      }
    }

    // If we're not guiding to a decoy and no commitment is active, clear the latch so a
    // future decoy can start a fresh commitment window.
    if (!committedDecoy && (!guide || guide->kind != CombatTargetKind::Decoy) && m.decoyCommitRemainingSimSec <= 0.0) {
      m.committedDecoyId = 0;
    }


    // --- Seeker gimbal dynamics (optional) ---
    //
    // Maintain an internal seeker pointing direction that slews toward the currently guided target.
    // This direction feeds into track quality (next frame) and helps model lock breaks during
    // high-LOS-rate situations.
    if (!m.hasSeekerDir) {
      m.hasSeekerDir = true;
      m.seekerDirWorld = fwd;
    }
    if (m.seekerDirWorld.lengthSq() < 1e-12) m.seekerDirWorld = fwd;
    else m.seekerDirWorld = m.seekerDirWorld.normalized();

    if (seekerActive && guide) {
      const math::Vec3d to = guide->centerKm - m.posKm;
      if (to.lengthSq() > 1e-12) {
        const math::Vec3d toDir = to.normalized();

        // Only slew when the guide is within the seeker cone. Memory guidance can point anywhere,
        // but the seeker itself should not "snap" to off-FOV predictions.
        const double cosAng = math::dot(fwd, toDir);
        if (cosAng >= seekerFovCos) {
          const double slew = std::max(0.0, m.seekerSlewRateRadS);
          if (slew <= 0.0) {
            m.seekerDirWorld = toDir;
          } else {
            m.seekerDirWorld = rotateTowards(m.seekerDirWorld, toDir, slew * dtSim);
            if (m.seekerDirWorld.lengthSq() < 1e-12) m.seekerDirWorld = toDir;
            else m.seekerDirWorld = m.seekerDirWorld.normalized();
          }
        }
      }
    }

    // Desired direction used by LeadPursuit mode (and as a fallback).
    if (guide) {
      const double leadSpeed = guidanceSpeedKmS(m, speed);

      const math::Vec3d tAcc = guide->accelKmS2;
      const bool useAccelLead = (math::isFinite(tAcc) && tAcc.lengthSq() > 1e-18);

      const auto lead =
          useAccelLead
              ? solveProjectileLeadAccel(m.posKm,
                                        /*shooterVelKmS=*/{0, 0, 0},
                                        guide->centerKm,
                                        guide->velKmS,
                                        tAcc,
                                        leadSpeed,
                                        /*maxTimeSec=*/1.0e6,
                                        /*minTimeSec=*/1.0e-3)
              : solveProjectileLead(m.posKm,
                                  /*shooterVelKmS=*/{0, 0, 0},
                                  guide->centerKm,
                                  guide->velKmS,
                                  leadSpeed,
                                  /*maxTimeSec=*/1.0e6,
                                  /*minTimeSec=*/1.0e-3);
      if (lead && lead->aimDirWorld.lengthSq() > 1e-12) {
        desiredDir = lead->aimDirWorld;
      } else {
        const math::Vec3d to = guide->centerKm - m.posKm;
        if (to.lengthSq() > 1e-12) desiredDir = to.normalized();
      }
    }


    const math::Vec3d vHat = v.normalized();
    const bool allowProNav =
        (m.guidance == MissileGuidance::ProNav) && guide && (guide->kind != CombatTargetKind::Decoy);

    const double maxTurnThisStep = std::max(0.0, turnRateRadS) * dtSim;

    // Compute the commanded direction as if the missile could instantly achieve the required turn
    // within this time step. If maxTurnAccelRadS2 > 0, an additional "autopilot lag" model will
    // limit how fast the missile can actually change turn rate.
    math::Vec3d newDirCmd = vHat;
    if (allowProNav) {
      math::Vec3d apnAccelKmS2 = m.targetAccelEstKmS2;
      // If the sensor provides a direct target acceleration estimate, prefer it over
      // the finite-difference APN estimate (more stable against noisy samples).
      if (math::isFinite(guide->accelKmS2) && guide->accelKmS2.lengthSq() > 1e-18) {
        apnAccelKmS2 = guide->accelKmS2;
      }
      double apnGain = std::max(0.0, m.apnTargetAccelGain);
      if (m.enableTrackQuality && seekerActive) apnGain *= m.trackQuality;
      const double apnMaxA = std::max(0.0, m.apnMaxTargetAccelKmS2);
      newDirCmd = proNavSteerDir(vHat,
                                /*missileVelKmS=*/m.velKmS,
                                /*missileSpeedKmS=*/speed,
                                /*turnRateRadS=*/turnRateRadS,
                                dtSim,
                                m.posKm,
                                *guide,
                                m.navConstant,
                                /*targetAccelKmS2=*/apnAccelKmS2,
                                /*apnTargetAccelGain=*/apnGain,
                                /*apnMaxTargetAccelKmS2=*/apnMaxA,
                                /*fallbackAimDir=*/desiredDir);
    } else {
      newDirCmd = rotateTowards(vHat, desiredDir, maxTurnThisStep);
    }

    // Optional: seeker lag can directly bias guidance (in addition to its effect on
    // track quality). This helps gimbal slew produce real miss distance in close-range
    // high-LOS-rate fights without introducing randomness.
    newDirCmd = applySeekerLagGuidanceDir(m,
                                         vHat,
                                         maxTurnThisStep,
                                         newDirCmd,
                                         seekerActive,
                                         seekerUpdateThisStep,
                                         lockedTarget,
                                         guide,
                                         seekerFovCos);


    // Optional: deterministic track error (sensor noise).
    // When enabled, poor track quality can translate into real miss distance even without decoys.
    newDirCmd = applyTrackErrorDir(m, newDirCmd, ageStart, seekerActive, guide);

    // Optional: deterministic swarm/formation coordination.
    // Apply a small bounded heading bias derived from nearby friendly missiles
    // pursuing the same target (reduces stacking and creates multi-axis attacks).
    newDirCmd = applySwarmBiasDir(m, mi, swarmSnap, vHat, newDirCmd, turnRateRadS, dtSim);

    // Optional: terminal weave (endgame jink) to complicate point defense.
    newDirCmd = applyTerminalWeaveDir(m, newDirCmd, guide, m.posKm, speedForTurn, ageStart, seekerActive);

    // Optional: deterministic asteroid avoidance. Bias the commanded direction
    // away from impending asteroid collisions inside a short look-ahead horizon.
    newDirCmd = applyAsteroidAvoidanceDir(m,
                                         m.posKm,
                                         vHat,
                                         speedForTurn,
                                         turnRateRadS,
                                         dtSim,
                                         newDirCmd,
                                         targets,
                                         targetCount);

    const double maxTurnAccel = std::max(0.0, m.maxTurnAccelRadS2);
    if (maxTurnAccel > 0.0 && dtSim > 1e-9) {
      if (!m.hasTurnOmega) {
        m.hasTurnOmega = true;
        m.turnOmegaWorld = {0, 0, 0};
      }

      // Convert the commanded direction into an angular-velocity command over this frame.
      const double c = std::clamp(math::dot(vHat, newDirCmd), -1.0, 1.0);
      const double ang = std::acos(c);
      math::Vec3d omegaCmd{0, 0, 0};
      if (ang > 1e-9) {
        math::Vec3d axis = math::cross(vHat, newDirCmd);
        const double al2 = axis.lengthSq();
        if (al2 > 1e-18) {
          axis = axis / std::sqrt(al2);
          omegaCmd = axis * (ang / dtSim);
        }
      }

      // Clamp commanded turn rate to the (possibly g-limited) maximum.
      const double maxTR = std::max(0.0, turnRateRadS);
      const double oc = omegaCmd.length();
      if (oc > maxTR && oc > 1e-12) omegaCmd *= (maxTR / oc);

      // Apply a per-frame angular acceleration limit.
      math::Vec3d omega = m.turnOmegaWorld;
      const math::Vec3d dOmega = omegaCmd - omega;
      const double dMag = dOmega.length();
      const double maxDelta = maxTurnAccel * dtSim;
      if (dMag > maxDelta && dMag > 1e-12) {
        omega += dOmega * (maxDelta / dMag);
      } else {
        omega = omegaCmd;
      }

      // Safety: clamp the resulting turn rate again.
      const double om = omega.length();
      if (om > maxTR && om > 1e-12) omega *= (maxTR / om);

      m.turnOmegaWorld = omega;

      // Integrate direction by applying the angular velocity.
      const double omegaMag = omega.length();
      math::Vec3d newDir = vHat;
      if (omegaMag > 1e-12) {
        const math::Vec3d axis = omega / omegaMag;
        newDir = rotateAroundAxis(vHat, axis, omegaMag * dtSim);
        if (newDir.lengthSq() < 1e-12) newDir = newDirCmd;
        else newDir = newDir.normalized();
      } else {
        newDir = newDirCmd;
      }

      m.velKmS = newDir * speed;
    } else {
      // If the lag model is disabled, keep behavior fully deterministic and clear state.
      m.hasTurnOmega = false;
      m.turnOmegaWorld = {0, 0, 0};
      m.velKmS = newDirCmd * speed;
    }

    // --- Motor model (optional) ---
    // Treat this as a 1D boost/coast along the current velocity direction.
    // Guidance determines direction; motor determines distance and next-frame speed.
    math::Vec3d dir = m.velKmS.normalized();
    if (dir.lengthSq() < 1e-12) dir = {0, 0, 1};

    const double accel = std::max(0.0, m.thrustAccelKmS2);
    const double maxSpeed = std::max(0.0, m.maxSpeedKmS);
    double burnRem = std::max(0.0, m.motorBurnRemainingSimSec);

    double stepDistKm = speed * dtSim;
    double speedEnd = speed;

    if (accel > 0.0 && burnRem > 0.0 && maxSpeed > speed + 1e-12) {
      const double burnDt = std::min(dtSim, burnRem);
      const double coastDt = dtSim - burnDt;

      // Distance traveled during the burn.
      const double tToMax = std::max(0.0, (maxSpeed - speed) / accel);
      if (tToMax >= burnDt) {
        // Never reaches max speed during this burn interval.
        speedEnd = speed + accel * burnDt;
        stepDistKm = speed * burnDt + 0.5 * accel * burnDt * burnDt;
      } else {
        // Accelerate until max speed, then cruise at max speed.
        speedEnd = maxSpeed;
        stepDistKm = speed * tToMax + 0.5 * accel * tToMax * tToMax;
        stepDistKm += maxSpeed * (burnDt - tToMax);
      }

      // Coast after the burn (if the burn ends before dtSim).
      if (coastDt > 0.0) {
        stepDistKm += speedEnd * coastDt;
      }

      burnRem -= burnDt;
      m.motorBurnRemainingSimSec = burnRem;
    }

    const math::Vec3d b = m.posKm + dir * stepDistKm;
    m.velKmS = dir * speedEnd;

    m.prevKm = a;
    m.posKm = b;
    m.ttlSimSec -= dtSim;
    m.ageSimSec = ageStart + dtSim;

    // --- Collision & detonation ---
    bool detonated = false;
    double detT = 1.0; // normalized time along [a,b] where the detonation occurs
    math::Vec3d detPoint = b;

    const SphereTarget* detTarget = nullptr;
    bool detOnContact = false;

    // Optional: point-defense projectile interception.
    bool detIntercepted = false;
    core::u64 detInterceptorShooterId = 0;
    bool detInterceptorFromPlayer = false;
    std::size_t detProjectileIndex = (std::size_t)-1;

    double bestT = 2.0; // [0,1]

    // --- Projectile interception (if enabled) ---
    if (interceptProjectiles && !interceptProjectiles->empty()) {
      for (std::size_t pi = 0; pi < interceptProjectiles->size(); ++pi) {
        Projectile& p = (*interceptProjectiles)[pi];
        if (p.ttlSimSec <= 0.0) continue;

        const double hitRadiusKm = std::max(0.0, m.radiusKm) + std::max(0.0, p.radiusKm);
        if (!(hitRadiusKm > 0.0)) continue;

        // Relative motion between the swept missile segment and the swept projectile segment.
        const math::Vec3d pa = p.prevKm;
        const math::Vec3d pb = p.posKm;

        const math::Vec3d relA = a - pa;
        const math::Vec3d relB = b - pb;

        double tEnter = 0.0;
        double tExit = 0.0;
        if (!math::segmentSphereIntersectionT(relA,
                                              relB,
                                              /*center=*/math::Vec3d{0, 0, 0},
                                              hitRadiusKm,
                                              tEnter,
                                              tExit)) {
          continue;
        }

        if (tEnter < bestT) {
          bestT = tEnter;
          detonated = true;
          detT = tEnter;
          detTarget = nullptr;
          detOnContact = false;

          detIntercepted = true;
          detInterceptorShooterId = p.shooterId;
          detInterceptorFromPlayer = p.fromPlayer;
          detProjectileIndex = pi;

          if (bestT <= 0.0) break;
        }
      }
    }

    // --- Target collision / proximity fuse ---
    if (targets && targetCount > 0) {
      const bool allowHitPlayer = !m.fromPlayer;
      const bool allowHitShips = true;
      const bool allowHitAsteroids = true;

      const double fuseExtraKm = std::max(0.0, m.proximityFuseKm);

      for (std::size_t ti = 0; ti < targetCount; ++ti) {
        const SphereTarget& t = targets[ti];

        if (!allowHitPlayer && t.kind == CombatTargetKind::Player) continue;
        if (!allowHitShips && t.kind == CombatTargetKind::Ship) continue;
        if (!allowHitAsteroids && t.kind == CombatTargetKind::Asteroid) continue;

        // Avoid immediate self-hits when the shooter has a real id.
        if (m.shooterId != 0 && t.id != 0 && t.id == m.shooterId) continue;

        const double baseHitRadiusKm = std::max(0.0, t.radiusKm) + std::max(0.0, m.radiusKm);

        // Linear relative motion: subtract target's swept center from the missile segment.
        const math::Vec3d c0 = t.centerKm;
        const math::Vec3d c1 = t.centerKm + t.velKmS * dtSim;

        const math::Vec3d relA = a - c0;
        const math::Vec3d relB = b - c1;

        // ---- Direct collision ----
        double tEnter = 0.0;
        double tExit = 0.0;
        if (math::segmentSphereIntersectionT(relA,
                                            relB,
                                            /*center=*/math::Vec3d{0, 0, 0},
                                            baseHitRadiusKm,
                                            tEnter,
                                            tExit)) {
          if (tEnter < bestT) {
            bestT = tEnter;
            detonated = true;
            detT = tEnter;
            detTarget = &t;
            detOnContact = true;

            detIntercepted = false;
            detInterceptorShooterId = 0;
            detInterceptorFromPlayer = false;
            detProjectileIndex = (std::size_t)-1;

            if (bestT <= 0.0) break;
          }
          continue;
        }

        // ---- Proximity fuse (standoff detonation) ----
        if (fuseExtraKm > 0.0) {
          const double fuseRadiusKm = baseHitRadiusKm + fuseExtraKm;

          if (math::segmentSphereIntersectionT(relA,
                                              relB,
                                              /*center=*/math::Vec3d{0, 0, 0},
                                              fuseRadiusKm,
                                              tEnter,
                                              tExit)) {
            // Detonate near closest approach (inside the intersection window) for maximum effect.
            double tClosest = math::segmentClosestT(relA, relB, math::Vec3d{0, 0, 0});
            const double tDet = std::clamp(tClosest, tEnter, tExit);

            if (tDet < bestT) {
              bestT = tDet;
              detonated = true;
              detT = tDet;
              detTarget = &t;
              detOnContact = false;

              detIntercepted = false;
              detInterceptorShooterId = 0;
              detInterceptorFromPlayer = false;
              detProjectileIndex = (std::size_t)-1;

              if (bestT <= 0.0) break;
            }
          }
        }
      }
    }

    if (detonated) {
      const math::Vec3d missileAtDet = a + (b - a) * detT;

      // If this was a projectile interception, consume the projectile and place the detonation at
      // the missile's position at the intercept time.
      if (detIntercepted && interceptProjectiles && detProjectileIndex != (std::size_t)-1 &&
          detProjectileIndex < interceptProjectiles->size()) {
        Projectile& p = (*interceptProjectiles)[detProjectileIndex];
        const math::Vec3d pAtDet = p.prevKm + (p.posKm - p.prevKm) * detT;
        p.posKm = pAtDet;
        p.ttlSimSec = 0.0;
        detPoint = missileAtDet;
      } else if (detTarget && detOnContact) {
        // For direct collisions, place the detonation on the target surface (instead of at the
        // missile center) so contact hits deliver full damage to the impacted target.
        const math::Vec3d centerAtDet = detTarget->centerKm + detTarget->velKmS * (dtSim * detT);
        const math::Vec3d toMissile = missileAtDet - centerAtDet;
        const math::Vec3d dir = math::safeNormalized(toMissile, math::Vec3d{0, 0, 1}, 1e-18);
        detPoint = centerAtDet + dir * std::max(0.0, detTarget->radiusKm);
      } else {
        detPoint = missileAtDet;
      }
    }

    // Expired missiles disappear without a detonation (keeps background clutter down).
    if (!detonated) {
      if (m.ttlSimSec <= 0.0) {
        m.ttlSimSec = 0.0;
      }
      continue;
    }

    // Detonate.
    MissileDetonation det{};
    det.pointKm = detPoint;
    det.blastRadiusKm = std::max(0.0, m.blastRadiusKm);
    det.baseDmg = std::max(0.0, m.dmg);
    det.fromPlayer = m.fromPlayer;
    det.shooterId = m.shooterId;
    det.intercepted = detIntercepted;
    det.interceptorShooterId = detInterceptorShooterId;
    det.interceptorFromPlayer = detInterceptorFromPlayer;
    outDetonations.push_back(det);

    // Splash hits.
    if (targets && targetCount > 0 && det.blastRadiusKm > 1e-6 && det.baseDmg > 1e-9) {
      const bool allowHitPlayer = !m.fromPlayer;
      const bool allowHitShips = true;
      const bool allowHitAsteroids = true;

      const math::Vec3d detFwd = math::safeNormalized(dir, math::Vec3d{0, 0, 1}, 1e-12);
      const double dirStrength = std::clamp(m.blastDirectionalStrength, 0.0, 1.0);
      const double dirMin = std::clamp(m.blastDirectionalMinFactor, 0.0, 1.0);

      const bool blastLos = m.blastRequireLineOfSight;
      const double blastPad = std::max(0.0, m.blastLineOfSightOcclusionPadKm);
      const double blastOccFactor = std::clamp(m.blastOccludedFactor, 0.0, 1.0);


      for (std::size_t ti = 0; ti < targetCount; ++ti) {
        const SphereTarget& t = targets[ti];

        if (!allowHitPlayer && t.kind == CombatTargetKind::Player) continue;
        if (!allowHitShips && t.kind == CombatTargetKind::Ship) continue;
        if (!allowHitAsteroids && t.kind == CombatTargetKind::Asteroid) continue;

        // Avoid self splash when shooter has a real id.
        if (m.shooterId != 0 && t.id != 0 && t.id == m.shooterId) continue;

        const math::Vec3d centerAtDet = t.centerKm + t.velKmS * (dtSim * detT);
        const math::Vec3d toVec = centerAtDet - detPoint;
        const double dist = toVec.length();
        const double surface = std::max(0.0, dist - std::max(0.0, t.radiusKm));
        double k = std::clamp(1.0 - surface / det.blastRadiusKm, 0.0, 1.0);
        if (k <= 0.0) continue;

        // Directional blast shaping (optional).
        if (dirStrength > 1e-9 && dist > 1e-9) {
          const math::Vec3d toDir = toVec / dist;
          const double cosAng = std::clamp(math::dot(detFwd, toDir), -1.0, 1.0);
          const double forward01 = std::clamp(0.5 * (1.0 + cosAng), 0.0, 1.0);
          const double shaped = dirMin + (1.0 - dirMin) * forward01;
          const double factor = std::clamp((1.0 - dirStrength) + dirStrength * shaped, 0.0, 1.0);
          k *= factor;
          if (k <= 0.0) continue;
        }

        // Asteroid line-of-sight occlusion for blast damage (optional).
        if (blastLos && dist > 1e-9) {
          const core::u64 ignoreAsteroidId = (t.kind == CombatTargetKind::Asteroid) ? t.id : 0;
          if (blastOccludedByAsteroidSpheres(detPoint, centerAtDet, targets, targetCount, blastPad, ignoreAsteroidId)) {
            k *= blastOccFactor;
            if (k <= 0.0) continue;
          }
        }

        MissileHit h{};
        h.kind = t.kind;
        h.targetIndex = t.index;
        h.targetId = t.id;
        h.pointKm = detPoint;
        h.dmg = det.baseDmg * k;
        h.fromPlayer = m.fromPlayer;
        h.shooterId = m.shooterId;
        outHits.push_back(h);
      }
    }

    m.ttlSimSec = 0.0;
  }
}


void stepMissiles(std::vector<Missile>& missiles,
                  double dtSim,
                  const SphereTarget* targets,
                  std::size_t targetCount,
                  std::vector<MissileDetonation>& outDetonations,
                  std::vector<MissileHit>& outHits) {
  stepMissilesImpl(missiles, /*interceptProjectiles=*/nullptr, dtSim, targets, targetCount, outDetonations, outHits);
}

void stepMissilesWithProjectileInterception(std::vector<Missile>& missiles,
                                            std::vector<Projectile>& projectiles,
                                            double dtSim,
                                            const SphereTarget* targets,
                                            std::size_t targetCount,
                                            std::vector<MissileDetonation>& outDetonations,
                                            std::vector<MissileHit>& outHits) {
  stepMissilesImpl(missiles, &projectiles, dtSim, targets, targetCount, outDetonations, outHits);
}


double weaponHeatDelta(const WeaponDef& w, int distributorMk) {
  // Better distributors run cooler.
  const int dMk = std::clamp(distributorMk, 1, 3);
  const double heatFactor = std::clamp(1.08 - 0.08 * (double)dMk, 0.80, 1.10);
  return std::max(0.0, w.heatPerShot * heatFactor);
}

FireResult tryFireWeapon(Ship& shooter,
                         WeaponType weapon,
                         double cooldownSimSec,
                         int distributorMk,
                         core::u64 shooterId,
                         bool fromPlayer,
                         const SphereTarget* beamTargets,
                         std::size_t beamTargetCount) {
  FireResult out{};
  out.fired = false;
  out.newCooldownSimSec = cooldownSimSec;

  if (cooldownSimSec > 0.0) return out;

  const WeaponDef& w = weaponDef(weapon);
  out.fired = true;
  out.newCooldownSimSec = w.cooldownSimSec;
  out.heatDelta = weaponHeatDelta(w, distributorMk);

  if (w.beam) {
    const math::Vec3d origin = shooter.positionKm();
    math::Vec3d dir = shooter.forward();
    if (dir.lengthSq() < 1e-12) dir = {0, 0, 1};
    dir = dir.normalized();

    const double rangeKm = std::max(0.0, w.rangeKm);
    const RaycastHit hit = raycastNearestSphereKm(origin, dir, rangeKm, beamTargets, beamTargetCount);

    out.hasBeam = true;
    out.beam.aKm = origin;
    out.beam.bKm = hit.pointKm;
    out.beam.r = w.r;
    out.beam.g = w.g;
    out.beam.b = w.b;

    out.hit = hit.hit;
    if (hit.hit) {
      out.hitKind = hit.kind;
      out.hitIndex = hit.index;
      out.hitId = hit.id;
      out.hitPointKm = hit.pointKm;
    }
    return out;
  }

  // Guided weapon (missile).
  if (w.guided) {
    const double muzzleSpeedKmS = std::max(1e-6, w.projSpeedKmS);
    const double rangeKm = std::max(0.0, w.rangeKm);

    math::Vec3d fwd = shooter.forward();
    if (fwd.lengthSq() < 1e-12) fwd = {0, 0, 1};
    fwd = fwd.normalized();

    // Lock onto a target in front of the shooter.
    bool hasLock = false;
    CombatTargetKind lockKind = CombatTargetKind::Ship;
    core::u64 lockId = 0;

    if (beamTargets && beamTargetCount > 0) {
      double bestAim = -1.0;
      double bestDist = 1.0e30;
      const math::Vec3d origin = shooter.positionKm();
      const math::Vec3d dir = fwd;

      for (std::size_t i = 0; i < beamTargetCount; ++i) {
        const SphereTarget& t = beamTargets[i];
        if (t.kind != CombatTargetKind::Ship && t.kind != CombatTargetKind::Player) continue;
        if (fromPlayer && t.kind == CombatTargetKind::Player) continue;
        if (shooterId != 0 && t.id != 0 && t.id == shooterId) continue;

        const math::Vec3d to = t.centerKm - origin;
        const double dist2 = to.lengthSq();
        if (dist2 < 1e-12) continue;
        const double dist = std::sqrt(dist2);
        if (dist > rangeKm) continue;

        const double aim = math::dot(dir, to / dist);
        // Don't lock behind.
        if (aim <= 0.0) continue;

        const double req = (t.minAimCos > -0.5) ? t.minAimCos : -1.0;
        if (aim < req) continue;

        // Prefer highest aim, then nearest.
        if (aim > bestAim + 1e-6 || (std::abs(aim - bestAim) <= 1e-6 && dist < bestDist)) {
          bestAim = aim;
          bestDist = dist;
          hasLock = true;
          lockKind = t.kind;
          lockId = t.id;
        }
      }
    }

    const math::Vec3d spawnKm = shooter.positionKm() + fwd * 520.0;

    Missile m{};
    m.prevKm = spawnKm;
    m.posKm = spawnKm;
    m.velKmS = shooter.velocityKmS() + fwd * muzzleSpeedKmS;
    m.r = w.r;
    m.g = w.g;
    m.b = w.b;
    m.radiusKm = 760.0;
    m.dmg = w.dmg;
    m.blastRadiusKm = std::max(0.0, w.blastRadiusKm);
    // Default proximity fuse: detonate slightly before direct impact when we have a blast radius.
    // This makes near-misses meaningful without changing pure direct-hit weapons.
    m.proximityFuseKm = (m.blastRadiusKm > 0.0) ? (0.45 * m.blastRadiusKm) : 0.0;

    // Warhead shaping: optional directional fragmentation + asteroid shadowing for blast damage.
    if (m.blastRadiusKm > 0.0) {
      // Allow asteroids to fully block splash damage (gameplay: "hide behind rocks").
      m.blastRequireLineOfSight = true;
      m.blastLineOfSightOcclusionPadKm = 0.0;
      m.blastOccludedFactor = 0.0;

      // Mild forward fragmentation bias (directional warhead approximation).
      if (weapon == WeaponType::RadarMissile) {
        m.blastDirectionalStrength = 0.65;
        m.blastDirectionalMinFactor = 0.20;
      } else {
        m.blastDirectionalStrength = 0.45;
        m.blastDirectionalMinFactor = 0.25;
      }
    }
    m.turnRateRadS = std::max(0.0, w.turnRateRadS);

    // Motor model: short boost phase (boost/coast) for a snappier intercept.
    // Defaults are conservative; missiles remain fully deterministic.
    if (weapon == WeaponType::RadarMissile) {
      m.thrustAccelKmS2 = 45.0;
      m.motorBurnRemainingSimSec = 1.40;
      const math::Vec3d maxVel = shooter.velocityKmS() + fwd * (muzzleSpeedKmS * 1.35);
      m.maxSpeedKmS = std::max(m.velKmS.length(), maxVel.length());
    } else {
      m.thrustAccelKmS2 = 35.0;
      m.motorBurnRemainingSimSec = 1.10;
      const math::Vec3d maxVel = shooter.velocityKmS() + fwd * (muzzleSpeedKmS * 1.25);
      m.maxSpeedKmS = std::max(m.velKmS.length(), maxVel.length());
    }

    // Turning "G-limit": as the missile accelerates, cap turn rate so lateral
    // acceleration doesn't implicitly grow with speed.
    m.maxLateralAccelKmS2 = muzzleSpeedKmS * m.turnRateRadS;

    // Autopilot lag: limit turn acceleration so missiles do not snap to max turn rate
    // instantly. This makes close-range evasive maneuvers (jink/beam) more meaningful.
    //
    // The value is expressed as rad/s^2 and scales with the missile's max turn rate.
    m.maxTurnAccelRadS2 = m.turnRateRadS * ((weapon == WeaponType::RadarMissile) ? 18.0 : 12.0);

    // Lifetime is computed from the motor profile so range remains stable even
    // when missiles accelerate.
    m.ttlSimSec = estimateBoostCoastTtlSimSec(rangeKm,
                                             std::max(1e-6, m.velKmS.length()),
                                             m.thrustAccelKmS2,
                                             m.maxSpeedKmS,
                                             m.motorBurnRemainingSimSec);

    // Seeker tuning (decoys + field-of-view).
    // Heat seekers can be lured by flares; radar seekers by chaff.
    // These parameters intentionally remain simple and deterministic.
    if (weapon == WeaponType::RadarMissile) {
      m.seeker = MissileSeekerType::Radar;
      // Simple two-phase behavior: midcourse (datalink) then terminal active seeker.
      // While the seeker is inactive, the missile ignores decoys and does not enforce FOV.
      m.seekerActivationSimSec = 0.35;
      // Midcourse datalink: while the active seeker is inactive, continue to receive
      // target updates from the launching platform (range + optional LOS gate).
      m.datalinkRangeKm = rangeKm;
      m.datalinkRequireLineOfSight = true;
      m.datalinkOcclusionPadKm = 0.0;
      // Once active, allow limited autonomous reacquisition if the original lock is lost.
      m.autoAcquireRangeKm = rangeKm;
      // Narrower seeker cone (radians -> cosine).
      m.seekerFovCos = std::cos(0.80);

      // Seeker gimbal slew: fast but not instantaneous.
      m.seekerSlewRateRadS = 9.0;
      // Slightly resistant to chaff, but can still be fooled.
      m.decoyResistance = 0.92;

      // LOS + radar notch: simple deterministic terrain-masking hooks.
      m.requireLineOfSight = true;
      m.lineOfSightOcclusionPadKm = 0.0;
      // Targets with very low radial velocity can be "notched" / beamed.
      m.radarDopplerNotchKmS = 0.75;

      // Home-On-Jam (HOJ): radar seekers can maintain a coarse track against
      // actively jamming targets even while inside the doppler notch.
      m.homeOnJam = true;
      m.homeOnJamMinJammerPower = 0.25;
      m.homeOnJamTrackQualityCap = 0.70;
      m.homeOnJamAcquireBias = 0.35;

      // Optional decoy discrimination gates: approximate range/angle/doppler track gates.
      // Radar seekers are fairly good at rejecting decoys not close to the true track.
      m.decoyAngleGateCos = std::cos(0.35);
      m.decoyRangeGateKm = 15000.0;
      m.decoyDopplerGateKmS = 3.0;

      // Track quality: makes chaff more effective after a lock break (e.g. notch/terrain).
      m.enableTrackQuality = true;
      m.trackQuality = 1.0;
      m.trackQualityRiseHalfLifeSimSec = 0.18;
      m.trackQualityFallHalfLifeSimSec = 0.35;
      m.trackQualityResistFloor = 0.35;

      // Countermeasure clutter suppression: strong chaff can degrade track quality
      // even before it fully steals lock (makes timed bursts matter).
      m.decoyClutterTrackSuppressionGain = 1.0;

      // Close-range burn-through: at short range, chaff/jamming becomes less
      // effective as the seeker can better discriminate the true target.
      m.decoyBurnThroughRangeKm = 0.12 * rangeKm;
      m.decoyBurnThroughMinFactor = 0.05;

      // Asteroid avoidance: deterministic steering bias away from imminent collisions.
      // Keeps missiles from suiciding into dense asteroid belts without a full path planner.
      m.asteroidAvoidanceStrength = 0.85;
      m.asteroidAvoidanceLookaheadSimSec = 1.0;
      m.asteroidAvoidancePadKm = 0.0;

      // Swarm coordination: encourage multi-axis attacks and avoid missiles
      // stacking perfectly on the same path.
      m.swarmSeparationStrength = 0.85;
      m.swarmCohesionStrength = 0.12;
      m.swarmAlignmentStrength = 0.18;
      m.swarmSeparationKm = 3.0 * std::max(0.0, m.radiusKm);
      m.swarmNeighborRangeKm = 14.0 * std::max(0.0, m.radiusKm);

      // Terminal weave: small endgame jink to complicate point defense.
      m.terminalWeaveAmplitudeRad = math::degToRad(2.25);
      m.terminalWeaveFrequencyHz = 2.2;
      m.terminalWeaveStartTtiSec = 1.6;
      m.terminalWeaveMaxRad = math::degToRad(3.0);
    } else {
      // Default guided weapon: heat seeker.
      m.seeker = MissileSeekerType::Heat;
      m.seekerActivationSimSec = 0.0;
      // Heat seekers require a clean launch lock; no autonomous reacquisition by default.
      m.autoAcquireRangeKm = 0.0;
      // Wider cone than radar.
      m.seekerFovCos = std::cos(0.95);

      // Heat seekers typically have a slower gimbal than active radar.
      m.seekerSlewRateRadS = 6.0;
      // Easier to decoy with flares.
      m.decoyResistance = 0.80;

      // IR aspect: tail-chase is hotter, head-on is cooler (engine plume).
      m.heatAspectFrontFactor = 0.75;
      m.heatAspectRearFactor = 1.25;
      m.heatAspectMinSpeedKmS = 0.02;
      m.heatAspectSpeedForFullKmS = 0.20;

      // Heat seekers still need a clear LOS to maintain track.
      m.requireLineOfSight = true;
      m.lineOfSightOcclusionPadKm = 0.0;

      // Track quality: makes flares more effective after brief FOV/LOS breaks.
      m.enableTrackQuality = true;
      m.trackQuality = 1.0;
      m.trackQualityRiseHalfLifeSimSec = 0.12;
      m.trackQualityFallHalfLifeSimSec = 0.22;
      m.trackQualityResistFloor = 0.20;

      // Countermeasure clutter suppression: flares can reduce effective track quality
      // even when they do not immediately win the decoy comparison.
      m.decoyClutterTrackSuppressionGain = 0.65;

      // Asteroid avoidance: heat seekers are a bit less aggressive.
      m.asteroidAvoidanceStrength = 0.65;
      m.asteroidAvoidanceLookaheadSimSec = 0.8;
      m.asteroidAvoidancePadKm = 0.0;

      // Swarm coordination: lighter than radar missiles (heat seekers are more
      // easily decoyed and shouldn't over-coordinate).
      m.swarmSeparationStrength = 0.55;
      m.swarmCohesionStrength = 0.08;
      m.swarmAlignmentStrength = 0.10;
      m.swarmSeparationKm = 2.6 * std::max(0.0, m.radiusKm);
      m.swarmNeighborRangeKm = 11.0 * std::max(0.0, m.radiusKm);
    }

    // Decoy commitment: once fooled by a countermeasure, keep steering toward it
    // for a short, weapon-tuned minimum time window (prevents jittery lock switching).
    m.decoyCommitSimSec = (weapon == WeaponType::RadarMissile) ? 0.15 : 0.45;

    // Target memory window (inertial midcourse): if the locked target briefly leaves
    // the seeker cone, keep guiding toward the last known kinematics for a short time.
    m.targetMemorySimSec = (weapon == WeaponType::RadarMissile) ? 2.25 : 0.85;

    // Guidance: radar missiles use proportional navigation for cleaner intercepts
    // against maneuvering targets, while heat seekers keep the older lead-pursuit
    // behavior (and remain more "flareable").
    if (weapon == WeaponType::RadarMissile) {
      m.guidance = MissileGuidance::ProNav;
      m.navConstant = 4.0;

      // APN (augmented proportional navigation): target-acceleration feed-forward
      // improves intercept quality against hard-turning targets. Classical scaling
      // is ~N/2; we also filter/clamp the finite-difference estimate for stability.
      m.apnTargetAccelGain = 0.5 * m.navConstant;
      m.apnMaxTargetAccelKmS2 = 2.5;
      m.apnAccelHalfLifeSimSec = 0.12;
    } else {
      m.guidance = MissileGuidance::LeadPursuit;
      m.navConstant = 3.5;
    }
    m.fromPlayer = fromPlayer;
    m.shooterId = shooterId;

    m.hasTarget = hasLock;
    m.targetKind = lockKind;
    m.targetId = lockId;

    out.hasMissile = true;
    out.missile = m;

    // Minimal recoil.
    shooter.setVelocityKmS(shooter.velocityKmS() - fwd * 0.0015);

    return out;
  }

  // Projectile weapon.
  const double muzzleSpeedKmS = std::max(1e-6, w.projSpeedKmS);
  const double rangeKm = std::max(0.0, w.rangeKm);
  const double ttlSim = rangeKm / muzzleSpeedKmS;

  math::Vec3d fwd = shooter.forward();
  if (fwd.lengthSq() < 1e-12) fwd = {0, 0, 1};
  fwd = fwd.normalized();

  math::Vec3d shotDir = fwd;

  // Soft aim assist for ballistic projectiles (Cannon / Railgun):
  // If the caller supplies candidate targets with a minAimCos cone filter,
  // gently bias the shot toward a lead solution for the most-aligned target.
  if (beamTargets && beamTargetCount > 0) {
    const math::Vec3d originKm = shooter.positionKm();
    const math::Vec3d shooterVelKmS = shooter.velocityKmS();

    const SphereTarget* best = nullptr;
    math::Vec3d bestDirTo{0, 0, 0};
    double bestAim = -1.0;
    double bestDist = 0.0;
    double bestCone = -1.0;

    for (std::size_t i = 0; i < beamTargetCount; ++i) {
      const SphereTarget& t = beamTargets[i];
      if (t.kind == CombatTargetKind::Decoy) continue;
      if (fromPlayer && t.kind == CombatTargetKind::Player) continue;

      // Avoid self-lock when ids are available.
      if (shooterId != 0 && t.id != 0 && t.id == shooterId) continue;

      const math::Vec3d to = t.centerKm - originKm;
      const double dist = to.length();
      if (!(dist > 1e-9)) continue;
      if (rangeKm > 0.0 && dist > rangeKm) continue;

      const math::Vec3d dirTo = to / dist;
      const double aim = math::dot(fwd, dirTo);
      if (!(aim > 0.0)) continue;

      const double cone = t.minAimCos;
      // Require an explicit cone filter to engage aim assist.
      if (cone < -0.5) continue;
      if (aim < cone) continue;

      if (!best || aim > bestAim + 1e-9 || (std::abs(aim - bestAim) <= 1e-9 && dist < bestDist)) {
        best = &t;
        bestDirTo = dirTo;
        bestAim = aim;
        bestDist = dist;
        bestCone = cone;
      }
    }

    if (best) {
      math::Vec3d aimDir = bestDirTo;
      if (const auto lead = solveProjectileLeadAccel(originKm,
                                               shooterVelKmS,
                                               best->centerKm,
                                               best->velKmS,
                                               best->accelKmS2,
                                               muzzleSpeedKmS,
                                               ttlSim)) {
        aimDir = lead->aimDirWorld;
      }

      // Smoothstep blend weight based on how centered the target already is.
      // (Prevents snaps at the edge of the assist cone.)
      const double denom = std::max(1e-9, 1.0 - bestCone);
      double t0 = (bestAim - bestCone) / denom;
      t0 = math::clamp(t0, 0.0, 1.0);
      double w0 = t0 * t0 * (3.0 - 2.0 * t0);

      // Also respect the cone for the lead direction (fast lateral targets can
      // push the lead point outside the reticle even if the center is inside).
      const double aimLead = math::dot(fwd, aimDir);
      double t1 = (aimLead - bestCone) / denom;
      t1 = math::clamp(t1, 0.0, 1.0);
      double w1 = t1 * t1 * (3.0 - 2.0 * t1);

      double w = std::min(w0, w1);

      // Allow the caller to scale aim assist strength per-target (e.g. sensor quality).
      w *= math::clamp(best->aimAssistWeight01, 0.0, 1.0);
      if (w > 1e-6 && math::isFinite(aimDir)) {
        const math::Vec3d blended = fwd * (1.0 - w) + aimDir * w;
        if (blended.lengthSq() > 1e-12) shotDir = blended.normalized();
      }
    }
  }

  const math::Vec3d spawnKm = shooter.positionKm() + shotDir * 400.0;

  Projectile p{};
  p.prevKm = spawnKm;
  p.posKm = spawnKm;
  p.velKmS = shooter.velocityKmS() + shotDir * muzzleSpeedKmS;
  p.r = w.r; p.g = w.g; p.b = w.b;
  p.ttlSimSec = ttlSim;
  p.radiusKm = (weapon == WeaponType::Railgun) ? 520.0 : 700.0;
  p.dmg = w.dmg;
  p.fromPlayer = fromPlayer;
  p.shooterId = shooterId;

  out.hasProjectile = true;
  out.projectile = p;

  // Small recoil impulse (matches historical tuning in stellar_game).
  const double recoil = (weapon == WeaponType::Railgun) ? 0.003 : 0.002;
  shooter.setVelocityKmS(shooter.velocityKmS() - shotDir * recoil);

  return out;
}

} // namespace stellar::sim
