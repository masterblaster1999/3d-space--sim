#include "stellar/sim/MissileDefense.h"

#include "stellar/math/Geometry.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace stellar::sim {

static double timeToSeekerActiveSec(const Missile& missile) {
  const double activateAfter = std::max(0.0, missile.seekerActivationSimSec);
  const double age = std::max(0.0, missile.ageSimSec);
  if (activateAfter <= 0.0) return 0.0;
  if (age >= activateAfter) return 0.0;
  return activateAfter - age;
}

static double timeToNextSeekerTickSec(const Missile& missile) {
  const double period = std::max(0.0, missile.seekerUpdatePeriodSimSec);
  if (!(period > 1e-9)) return 0.0; // continuous/legacy

  const double activateAfter = std::max(0.0, missile.seekerActivationSimSec);
  const double age = std::max(0.0, missile.ageSimSec);

  if (age < activateAfter) {
    // First tick coincides with seeker activation.
    return std::max(0.0, activateAfter - age);
  }

  // Seeker is already active; compute phase in the scan cycle.
  const double tSince = std::max(0.0, age - activateAfter);
  const double phase = std::fmod(tSince, period);
  const double eps = 1e-9;
  if (!std::isfinite(phase)) return 0.0;
  if (phase <= eps || period - phase <= eps) return 0.0;
  return std::max(0.0, period - phase);
}

CountermeasureResponsePlan planCountermeasureResponse(const Missile& missile,
                                                      const math::Vec3d& targetPosKm,
                                                      const math::Vec3d& targetVelKmS,
                                                      const CountermeasureResponseParams& params) {
  CountermeasureResponsePlan out{};

  // Pick countermeasure type based on seeker.
  out.type = (missile.seeker == MissileSeekerType::Radar) ? CountermeasureType::Chaff
                                                          : CountermeasureType::Flare;

  // Estimate time-to-impact from relative kinematics.
  const math::Vec3d toTarget = targetPosKm - missile.posKm;
  const double distKm = toTarget.length();
  const math::Vec3d toDir = math::safeNormalized(toTarget, math::Vec3d{0, 0, 1}, 1e-12);
  const math::Vec3d relVel = missile.velKmS - targetVelKmS;
  const double closingKmS = math::dot(relVel, toDir);

  if (!(closingKmS > 1e-9) || !(distKm >= 0.0) || !std::isfinite(distKm)) {
    return out;
  }

  double tti = distKm / closingKmS;
  if (!std::isfinite(tti) || tti < 0.0) tti = 0.0;
  out.timeToImpactSec = tti;

  // If the missile expires sooner than our simple range-rate estimate, clamp.
  if (missile.ttlSimSec > 0.0) {
    out.timeToImpactSec = std::min(out.timeToImpactSec, std::max(0.0, missile.ttlSimSec));
    tti = out.timeToImpactSec;
  }

  const double tActive = timeToSeekerActiveSec(missile);
  out.timeToSeekerActiveSec = tActive;

  // If the seeker won't be active before impact, countermeasures can't influence guidance.
  if (tti <= tActive + 1e-9) {
    return out;
  }

  // Panic: deploy immediately when impact is imminent and the seeker is active.
  const double panic = std::max(0.0, params.panicTtiSec);
  if (tti <= panic && tActive <= 1e-9) {
    out.valid = true;
    out.firstReleaseDelaySec = 0.0;
    out.repeatEverySec = 0.0;
    out.burstCount = 1;
    out.timeToFirstSeekerTickSec = 0.0;
    return out;
  }

  // Don't start deploying too early: use a simple "window before impact" heuristic.
  const double window = std::max(0.0, params.decoyWindowSec);
  double earliestInfluence = 0.0;
  if (window > 1e-9) {
    earliestInfluence = std::max(0.0, tti - window);
  }

  // Can't influence before seeker activation.
  earliestInfluence = std::max(earliestInfluence, tActive);

  const double period = std::max(0.0, missile.seekerUpdatePeriodSimSec);
  const double lead = std::max(0.0, params.releaseLeadSec);

  double firstTick = 0.0;
  if (!(period > 1e-9)) {
    // Continuous seeker: release at the desired influence start time.
    firstTick = earliestInfluence;
  } else {
    // Discrete seeker: align to the first scan tick at/after earliestInfluence.
    const double t0 = timeToNextSeekerTickSec(missile);
    if (t0 < 0.0 || !std::isfinite(t0)) return out;

    if (earliestInfluence <= t0 + 1e-12) {
      firstTick = t0;
    } else {
      const double dt = earliestInfluence - t0;
      const double n = std::ceil(dt / period);
      if (!std::isfinite(n) || n < 0.0) {
        firstTick = t0;
      } else {
        firstTick = t0 + n * period;
      }
    }
  }

  out.timeToFirstSeekerTickSec = firstTick;

  // No useful scan before impact.
  if (firstTick >= tti - 1e-9) {
    return out;
  }

  out.firstReleaseDelaySec = std::max(0.0, firstTick - lead);
  out.valid = true;

  // Repeated bursts for discrete seekers: one per scan tick (up to maxBursts).
  out.repeatEverySec = 0.0;
  out.burstCount = 1;

  if (params.enableRepeatBursts && period > 1e-9) {
    const int maxB = std::max(1, params.maxBursts);
    const double remain = std::max(0.0, tti - firstTick);
    const int available = (int)std::floor(remain / period) + 1;
    out.burstCount = std::clamp(available, 1, maxB);
    if (out.burstCount > 1) out.repeatEverySec = period;
  }

  return out;
}

MissileThreatSummary nearestInboundMissile(const Missile* missiles,
                                          std::size_t missileCount,
                                          CombatTargetKind targetKind,
                                          core::u64 targetId,
                                          const math::Vec3d& targetPosKm,
                                          const math::Vec3d& targetVelKmS,
                                          const MissileThreatParams& params) {
  MissileThreatSummary best{};
  double bestTtiSec = std::numeric_limits<double>::infinity();

  if (!missiles || missileCount == 0) return best;

  const double minCos = std::clamp(params.minApproachCos, -1.0, 1.0);
  const double minClosing = std::max(0.0, params.minClosingKmS);
  const double maxDist = std::max(0.0, params.maxConsiderDistKm);

  for (std::size_t i = 0; i < missileCount; ++i) {
    const Missile& m = missiles[i];

    if (!m.hasTarget) continue;
    if (m.targetKind != targetKind) continue;
    if (m.targetId != targetId) continue;
    if (m.ttlSimSec <= 0.0) continue;

    const math::Vec3d toTarget = targetPosKm - m.posKm;
    const double distKm = toTarget.length();
    if (distKm > maxDist) continue;

    const math::Vec3d toDir = math::safeNormalized(toTarget, math::Vec3d{0, 0, 1}, 1e-12);
    const math::Vec3d mvDir = math::safeNormalized(m.velKmS, math::Vec3d{0, 0, 1}, 1e-12);
    const double approachCos = math::dot(mvDir, toDir);
    if (approachCos < minCos) continue;

    // Relative closing speed along the line-of-sight.
    const math::Vec3d relVel = m.velKmS - targetVelKmS;
    const double closingKmS = math::dot(relVel, toDir);
    if (closingKmS <= minClosing) continue;

    const double ttiSec = (distKm > 1e-9) ? (distKm / std::max(closingKmS, 1e-9)) : 0.0;
    if (ttiSec < bestTtiSec) {
      bestTtiSec = ttiSec;

      best.inbound = true;
      best.seeker = m.seeker;
      best.distKm = distKm;
      best.closingKmS = closingKmS;
      best.ttiSec = ttiSec;
      best.approachCos = approachCos;
      best.missileIndex = i;
      best.shooterId = m.shooterId;
      best.fromPlayer = m.fromPlayer;
    }
  }

  return best;
}


static double receivedJamming01ForDefense(const Missile& m, double distKm, double jammerPower) {
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

JammerResponsePlan planJammerResponse(const Missile& missile,
                                     const math::Vec3d& targetPosKm,
                                     const math::Vec3d& targetVelKmS,
                                     const JammerResponseParams& params) {
  JammerResponsePlan out{};

  if (!params.enable) return out;

  if (missile.seeker != MissileSeekerType::Radar) {
    out.reason = JammerResponseReason::NotRadarThreat;
    return out;
  }

  const math::Vec3d toTarget = targetPosKm - missile.posKm;
  const double distKm = toTarget.length();
  const math::Vec3d toDir = math::safeNormalized(toTarget, math::Vec3d{0, 0, 1}, 1e-12);

  const math::Vec3d relVel = missile.velKmS - targetVelKmS;
  const double closingKmS = math::dot(relVel, toDir);

  out.closingKmS = closingKmS;
  out.rangeRateKmS = math::dot(targetVelKmS - missile.velKmS, toDir);

  if (!(closingKmS > 1e-9) || !(distKm >= 0.0) || !std::isfinite(distKm)) {
    out.reason = JammerResponseReason::NotInbound;
    return out;
  }

  double tti = distKm / closingKmS;
  if (!std::isfinite(tti) || tti < 0.0) tti = 0.0;

  // Clamp by TTL (if known) to avoid recommending late jamming against a missile that will expire.
  if (missile.ttlSimSec > 0.0) {
    tti = std::min(tti, std::max(0.0, missile.ttlSimSec));
  }

  out.timeToImpactSec = tti;

  const double maxTti = params.maxTtiToJamSec;
  if (maxTti > 0.0 && tti > maxTti) {
    out.valid = true;
    out.reason = JammerResponseReason::TooFar;
    return out;
  }

  const double jp = std::max(0.0, params.jammerPower);
  out.receivedJamming01 = receivedJamming01ForDefense(missile, distKm, jp);

  const double minJ01 = std::clamp(params.minReceivedJamming01, 0.0, 1.0);
  if (out.receivedJamming01 < minJ01) {
    out.valid = true;
    out.reason = JammerResponseReason::TooWeak;
    return out;
  }

  // Determine whether the missile could actually engage HOJ given this jammer power.
  out.homeOnJamCapable = false;
  if (missile.homeOnJam && jp >= std::max(0.0, missile.homeOnJamMinJammerPower)) {
    const double minHojJ01 = std::clamp(missile.homeOnJamMinJamming01, 0.0, 1.0);
    if (out.receivedJamming01 >= minHojJ01 - 1e-12) {
      out.homeOnJamCapable = true;
    }
  }

  bool jammerOn = true;
  JammerResponseReason reason = JammerResponseReason::None;

  if (params.denyHomeOnJam && out.homeOnJamCapable) {
    const double denyClose = std::max(0.0, params.denyHomeOnJamTtiSec);
    if (denyClose > 0.0 && tti <= denyClose) {
      jammerOn = false;
      reason = JammerResponseReason::DenyHomeOnJamClose;
    }

    const double notch = std::max(0.0, missile.radarDopplerNotchKmS);
    const double mul = std::max(0.0, params.denyHomeOnJamNotchMultiple);

    if (jammerOn && notch > 0.0 && mul > 0.0) {
      if (std::fabs(out.rangeRateKmS) <= notch * mul) {
        jammerOn = false;
        reason = JammerResponseReason::DenyHomeOnJamNotch;
      }
    }
  }

  out.valid = true;
  out.jammerOn = jammerOn;
  out.jammerPower = jammerOn ? jp : 0.0;
  out.reason = reason;
  return out;
}


static math::Vec3d unitPerpFromSeed(const math::Vec3d& axisIn, core::u64 seed) {
  math::Vec3d axis = axisIn;
  if (axis.lengthSq() < 1e-12) axis = {0, 0, 1};
  axis = axis.normalized();

  math::Vec3d b = math::cross(axis, math::Vec3d{0, 1, 0});
  if (b.lengthSq() < 1e-12) b = math::cross(axis, math::Vec3d{1, 0, 0});
  if (b.lengthSq() < 1e-12) b = math::cross(axis, math::Vec3d{0, 0, 1});
  if (b.lengthSq() < 1e-12) return {1, 0, 0};
  b = b.normalized();

  math::Vec3d c = math::cross(axis, b);
  if (c.lengthSq() < 1e-12) return b;
  c = c.normalized();

  core::SplitMix64 rng(seed);
  const double ang = rng.range(0.0, 2.0 * std::numbers::pi);
  const double cs = std::cos(ang);
  const double sn = std::sin(ang);
  return math::safeNormalized(b * cs + c * sn, b, 1e-12);
}

MissileEvasionPlan planMissileEvasion(const Missile& missile,
                                     const math::Vec3d& targetPosKm,
                                     const math::Vec3d& targetVelKmS,
                                     core::u64 seed,
                                     const MissileEvasionParams& params) {
  MissileEvasionPlan out{};

  // Relative motion: missile in the target frame.
  const math::Vec3d r = missile.posKm - targetPosKm;
  const math::Vec3d v = missile.velKmS - targetVelKmS;

  const double distSq = r.lengthSq();
  const double dist = std::sqrt(std::max(0.0, distSq));
  const math::Vec3d los = math::safeNormalized(r, math::Vec3d{0, 0, 1}, 1e-12);

  const double vSq = v.lengthSq();
  const double vLen = std::sqrt(std::max(0.0, vSq));

  // Closing speed along LOS (positive when inbound).
  const double closing = -math::dot(v, los);
  out.closingKmS = closing;

  const double minRel = std::max(0.0, params.minRelSpeedKmS);
  if (!(vSq > 1e-18) || vLen < minRel || dist < 1e-9) {
    // Degenerate: pick any direction perpendicular to LOS.
    out.dirWorld = unitPerpFromSeed(los, seed);
    out.valid = (out.dirWorld.lengthSq() > 1e-12);
    return out;
  }

  if (closing <= 0.0) {
    // Not inbound (or grazing away).
    return out;
  }

  // Time of closest approach under constant relative velocity (clamped to now..inf).
  double tClosest = -math::dot(r, v) / vSq;
  if (!std::isfinite(tClosest)) tClosest = 0.0;
  if (tClosest < 0.0) tClosest = 0.0;
  out.tClosestSec = tClosest;

  // Predicted relative offset at closest approach.
  math::Vec3d missVec = r + v * tClosest;
  // Numerically ensure perpendicular to v (especially if tClosest was clamped).
  missVec = missVec - v * (math::dot(missVec, v) / vSq);

  const double missDist = missVec.length();
  out.missDistanceKm = missDist;

  math::Vec3d dir{0, 0, 0};
  if (missDist > std::max(0.0, params.minMissVecKm)) {
    // Move away from the predicted closest-approach point.
    dir = (-missVec) / missDist;
  } else {
    // Head-on: any direction perpendicular to relative velocity works; break ties via seed.
    dir = unitPerpFromSeed(v, seed);
  }

  // Optional: radar "beaming" bias.
  //
  // If the inbound missile is a radar seeker that models a doppler notch, bias the
  // jink direction toward rotating the target's velocity into the plane perpendicular
  // to the current line-of-sight. This increases LOS rate / reduces range-rate in a
  // deterministic way and synergizes with the notch mechanic.
  if (params.enableRadarBeaming && missile.seeker == MissileSeekerType::Radar) {
    const double notch = std::max(0.0, missile.radarDopplerNotchKmS);
    const double blend = std::clamp(params.radarBeamBlend, 0.0, 1.0);
    const double engageMul = std::max(0.0, params.radarBeamEngageNotchMultiple);
    if (notch > 0.0 && blend > 0.0) {
      // LOS from missile to target.
      const math::Vec3d toDir = -los;
      const double vrKmS = math::dot(targetVelKmS - missile.velKmS, toDir);
      const double absVr = std::fabs(vrKmS);

      if (absVr > notch * engageMul) {
        // Desired velocity direction: current velocity projected into the LOS-perpendicular plane.
        const double vT = targetVelKmS.length();

        math::Vec3d vPerp = targetVelKmS - toDir * math::dot(targetVelKmS, toDir);
        if (vPerp.lengthSq() < 1e-12) {
          // Degenerate (already LOS-aligned or nearly zero speed): pick any perpendicular direction.
          const math::Vec3d u = unitPerpFromSeed(toDir, seed ^ 0x6d4f2e31b7a3c0d5ull);
          vPerp = (vT > 1e-9) ? (u * vT) : u;
        } else {
          vPerp = vPerp.normalized() * std::max(0.0, vT);
        }

        const math::Vec3d deltaV = vPerp - targetVelKmS;
        const math::Vec3d fallback = unitPerpFromSeed(toDir, seed ^ 0x9e3779b97f4a7c15ull);
        const math::Vec3d beamDir = math::safeNormalized(deltaV, fallback, 1e-12);

        // Engagement weight ramps up as we exceed the notch threshold.
        const double w = std::clamp((absVr - notch) / std::max(absVr, 1e-9), 0.0, 1.0) * blend;

        dir = math::safeNormalized(dir * (1.0 - w) + beamDir * w, dir, 1e-12);
      }
    }
  }

  if (params.enforceLateralToLos) {
    // Project into the plane perpendicular to LOS (lateral jink).
    dir = dir - los * math::dot(dir, los);
    dir = math::safeNormalized(dir, unitPerpFromSeed(los, seed ^ 0x9e3779b97f4a7c15ull), 1e-12);
  } else {
    dir = math::safeNormalized(dir, unitPerpFromSeed(los, seed ^ 0x9e3779b97f4a7c15ull), 1e-12);
  }

  out.dirWorld = dir;
  out.valid = (dir.lengthSq() > 1e-12);
  return out;
}


static bool occludedByAsteroidSpheres(const math::Vec3d& fromKm,
                                     const math::Vec3d& toKm,
                                     const SphereTarget* targets,
                                     std::size_t targetCount,
                                     double radiusPadKm,
                                     core::u64 ignoreAsteroidId) {
  if (!targets || targetCount == 0) return false;

  const double pad = std::max(0.0, radiusPadKm);

  for (std::size_t i = 0; i < targetCount; ++i) {
    const SphereTarget& t = targets[i];
    if (t.kind != CombatTargetKind::Asteroid) continue;
    if (ignoreAsteroidId != 0 && t.id == ignoreAsteroidId) continue;

    const double r = std::max(0.0, t.radiusKm) + pad;
    if (r <= 0.0) continue;

    if (math::segmentHitsSphere(fromKm, toKm, t.centerKm, r)) {
      return true;
    }
  }

  return false;
}

MissileTerrainMaskPlan planMissileTerrainMasking(const Missile& missile,
                                                 const math::Vec3d& targetPosKm,
                                                 const math::Vec3d& targetVelKmS,
                                                 const SphereTarget* targets,
                                                 std::size_t targetCount,
                                                 core::u64 seed,
                                                 const MissileTerrainMaskParams& params) {
  MissileTerrainMaskPlan out{};

  if (!targets || targetCount == 0) return out;

  // If we're already occluded, we can keep our current plan (caller may re-plan later).
  if (params.ignoreIfAlreadyOccluded) {
    if (occludedByAsteroidSpheres(missile.posKm, targetPosKm, targets, targetCount,
                                 /*radiusPadKm=*/0.0, /*ignoreAsteroidId=*/0)) {
      return out;
    }
  }

  // Predict missile position slightly ahead to account for high missile speed.
  // If the relative geometry does not imply an inbound threat, fall back to no-lookahead.
  double lookahead = std::max(0.0, params.lookaheadSec);
  if (lookahead > 0.0) {
    const math::Vec3d toTarget = targetPosKm - missile.posKm;
    const double distKm = toTarget.length();
    const math::Vec3d toDir = math::safeNormalized(toTarget, math::Vec3d{0, 0, 1}, 1e-12);
    const math::Vec3d relVel = missile.velKmS - targetVelKmS;
    const double closingKmS = math::dot(relVel, toDir);
    if (!(closingKmS > 1e-9) || !(distKm >= 0.0) || !std::isfinite(distKm)) {
      lookahead = 0.0;
    } else {
      const double tti = distKm / closingKmS;
      if (std::isfinite(tti)) lookahead = std::min(lookahead, std::max(0.0, tti));
    }
  }

  const math::Vec3d missilePred = missile.posKm + missile.velKmS * lookahead;

  const double maxTravel = std::max(0.0, params.maxCoverTravelKm);
  const double coverPad = std::max(0.0, params.coverPadKm);

  double bestTravel = std::numeric_limits<double>::infinity();
  core::u64 bestId = 0;
  math::Vec3d bestPoint{0, 0, 0};

  for (std::size_t i = 0; i < targetCount; ++i) {
    const SphereTarget& a = targets[i];
    if (a.kind != CombatTargetKind::Asteroid) continue;

    const double r = std::max(0.0, a.radiusKm);
    if (r <= 0.0) continue;

    const math::Vec3d toAst = a.centerKm - missilePred;
    math::Vec3d u = math::safeNormalized(toAst, math::Vec3d{0, 0, 1}, 1e-12);

    // Cover point is on the far side of the asteroid relative to the missile.
    const math::Vec3d coverPoint = a.centerKm + u * (r + coverPad);

    const math::Vec3d d = coverPoint - targetPosKm;
    const double travel = d.length();
    if (!(travel >= 0.0) || !std::isfinite(travel)) continue;
    if (travel > maxTravel) continue;

    // Sanity: the cover point should actually be occluded by some asteroid from the
    // predicted missile position (usually the chosen one). Use a small pad so touches
    // count as occlusion.
    if (!occludedByAsteroidSpheres(missilePred, coverPoint, targets, targetCount,
                                   /*radiusPadKm=*/0.0, /*ignoreAsteroidId=*/0)) {
      continue;
    }

    // Prefer the nearest cover point (less travel).
    if (travel + 1e-9 < bestTravel) {
      bestTravel = travel;
      bestId = a.id;
      bestPoint = coverPoint;
    } else if (std::fabs(travel - bestTravel) <= 1e-9) {
      // Deterministic tie-break: prefer lower id; if ids are equal or missing, break via seed.
      if (bestId == 0 || (a.id != 0 && a.id < bestId)) {
        bestTravel = travel;
        bestId = a.id;
        bestPoint = coverPoint;
      } else if (a.id == bestId && bestId == 0) {
        // Extremely degenerate: random-but-deterministic pick.
        core::SplitMix64 rng(seed);
        if ((rng.nextU64() & 1ull) == 0ull) {
          bestTravel = travel;
          bestId = a.id;
          bestPoint = coverPoint;
        }
      }
    }
  }

  if (!std::isfinite(bestTravel) || bestTravel == std::numeric_limits<double>::infinity()) {
    return out;
  }

  math::Vec3d dir = bestPoint - targetPosKm;
  if (dir.lengthSq() < 1e-18) {
    // Already at (or extremely near) the cover point. Pick a lateral direction to avoid
    // feeding a zero vector into caller steering.
    dir = unitPerpFromSeed(missilePred - targetPosKm, seed ^ 0x6d4f2e31b7a3c0d5ull);
  } else {
    dir = dir.normalized();
  }

  out.valid = (dir.lengthSq() > 1e-12);
  out.dirWorld = dir;
  out.asteroidId = bestId;
  out.coverPointKm = bestPoint;
  out.travelKm = bestTravel;
  return out;
}


IntegratedDefensePlan planIntegratedDefense(const Missile* missiles,
                                            std::size_t missileCount,
                                            CombatTargetKind targetKind,
                                            core::u64 targetId,
                                            const math::Vec3d& targetPosKm,
                                            const math::Vec3d& targetVelKmS,
                                            const SphereTarget* worldTargets,
                                            std::size_t worldTargetCount,
                                            core::u64 seed,
                                            const IntegratedDefenseParams& params) {
  IntegratedDefensePlan out{};

  if (!missiles || missileCount == 0) return out;

  out.threat = nearestInboundMissile(missiles,
                                    missileCount,
                                    targetKind,
                                    targetId,
                                    targetPosKm,
                                    targetVelKmS,
                                    params.threat);

  if (!out.threat.inbound) return out;
  if (out.threat.missileIndex >= missileCount) return IntegratedDefensePlan{};

  const Missile& m = missiles[out.threat.missileIndex];
  out.valid = true;

  // Urgency mapping based on time-to-impact.
  const double scale = std::max(1e-6, params.urgencyTimeScaleSec);
  out.urgency01 = std::clamp(1.0 - out.threat.ttiSec / scale, 0.0, 1.0);

  // Always compute sub-plans for callers that want to introspect.
  out.countermeasures = planCountermeasureResponse(m, targetPosKm, targetVelKmS, params.countermeasures);
  out.jammer = planJammerResponse(m, targetPosKm, targetVelKmS, params.jammer);
  out.evasion = planMissileEvasion(m, targetPosKm, targetVelKmS, seed, params.evasion);

  // Terrain masking is only meaningful when the inbound missile depends on LOS
  // either during seeker-active tracking, or during pre-activation datalink updates.
  bool considerTerrain = false;
  if (params.enableTerrainMasking) {
    const double tActive = timeToSeekerActiveSec(m);
    considerTerrain = m.requireLineOfSight || (tActive > 1e-9 && m.datalinkRequireLineOfSight);
    if (considerTerrain) {
      out.terrain = planMissileTerrainMasking(m,
                                             targetPosKm,
                                             targetVelKmS,
                                             worldTargets,
                                             worldTargetCount,
                                             seed,
                                             params.terrain);
    }
  }

  const double minTti = std::max(0.0, params.minTtiForTerrainSec);
  const bool chooseTerrain = considerTerrain && out.terrain.valid && (out.threat.ttiSec > minTti);

  if (chooseTerrain) {
    out.maneuver = DefenseManeuverKind::TerrainMask;
    out.maneuverDirWorld = out.terrain.dirWorld;
  } else if (out.evasion.valid) {
    out.maneuver = DefenseManeuverKind::EvasionJink;
    out.maneuverDirWorld = out.evasion.dirWorld;
  } else {
    out.maneuver = DefenseManeuverKind::None;
    out.maneuverDirWorld = {0, 0, 0};
  }

  // If we are intentionally breaking LOS against a HOJ-capable radar missile,
  // keep the jammer off to avoid providing a strong HOJ cue on re-acquisition.
  if (chooseTerrain && out.jammer.valid && out.jammer.jammerOn && out.jammer.homeOnJamCapable) {
    out.jammer.jammerOn = false;
    out.jammer.jammerPower = 0.0;
    out.jammer.reason = JammerResponseReason::DenyHomeOnJamTerrainMask;
  }

  if (out.maneuver != DefenseManeuverKind::None) {
    if (out.maneuverDirWorld.lengthSq() < 1e-18) {
      out.maneuver = DefenseManeuverKind::None;
      out.maneuverDirWorld = {0, 0, 0};
    } else {
      out.maneuverDirWorld = out.maneuverDirWorld.normalized();
    }
  }

  return out;
}

IntegratedDefenseWithPointDefensePlan planIntegratedDefenseWithPointDefense(const Missile* missiles,
                                                                            std::size_t missileCount,
                                                                            CombatTargetKind targetKind,
                                                                            core::u64 targetId,
                                                                            const math::Vec3d& targetPosKm,
                                                                            const math::Vec3d& targetVelKmS,
                                                                            const math::Vec3d& turretAimDirWorld,
                                                                            const SphereTarget* worldTargets,
                                                                            std::size_t worldTargetCount,
                                                                            core::u64 seed,
                                                                            const IntegratedDefenseWithPointDefenseParams& params) {
  IntegratedDefenseWithPointDefensePlan out{};

  out.defense = planIntegratedDefense(missiles,
                                     missileCount,
                                     targetKind,
                                     targetId,
                                     targetPosKm,
                                     targetVelKmS,
                                     worldTargets,
                                     worldTargetCount,
                                     seed,
                                     params.defense);
  out.valid = out.defense.valid;

  if (!out.valid) return out;
  if (!params.enablePointDefense) return out;

  if (!missiles || missileCount == 0) return out;
  if (!out.defense.threat.inbound) return out;

  if (out.defense.threat.missileIndex >= missileCount) return IntegratedDefenseWithPointDefensePlan{};

  const double maxTti = std::max(0.0, params.maxTtiForPointDefenseSec);
  if (maxTti > 1e-9 && std::isfinite(out.defense.threat.ttiSec) && out.defense.threat.ttiSec > maxTti) {
    return out;
  }

  // Optionally clamp the PD planning horizon to the inbound time-to-impact estimate. This avoids
  // searching for shots that would occur after a likely impact.
  PointDefenseParams pd = params.pointDefense;
  if (params.clampPointDefenseHorizonToTti && std::isfinite(out.defense.threat.ttiSec) && out.defense.threat.ttiSec > 0.0) {
    const double tti = std::max(0.0, out.defense.threat.ttiSec);
    if (pd.maxPlanSimSec > 0.0) {
      pd.maxPlanSimSec = std::min(pd.maxPlanSimSec, tti);
    } else {
      pd.maxPlanSimSec = tti;
    }
  }

  const Missile& m = missiles[out.defense.threat.missileIndex];

  out.pointDefense = planPointDefenseShot(m,
                                         /*shooterPosKm=*/targetPosKm,
                                         /*shooterVelKmS=*/targetVelKmS,
                                         turretAimDirWorld,
                                         worldTargets,
                                         worldTargetCount,
                                         seed,
                                         pd);

  // Policy gate: don't recommend a shot that lands after the time-to-impact estimate.
  if (out.pointDefense.valid && std::isfinite(out.defense.threat.ttiSec)) {
    const double margin = std::max(0.0, params.interceptBeforeImpactMarginSec);
    const double latest = std::max(0.0, out.defense.threat.ttiSec - margin);
    if (out.pointDefense.interceptTimeSec > latest + 1e-9) {
      out.pointDefense = PointDefenseSolution{};
    }
  }

  return out;
}



// -----------------------------------------------------------------------------
// Multi-threat point defense engagement scheduling
// -----------------------------------------------------------------------------

static bool computeInboundThreatSummary(const Missile& m,
                                       std::size_t missileIndex,
                                       CombatTargetKind targetKind,
                                       core::u64 targetId,
                                       const math::Vec3d& targetPosKm,
                                       const math::Vec3d& targetVelKmS,
                                       const MissileThreatParams& params,
                                       MissileThreatSummary* out) {
  if (!out) return false;
  *out = MissileThreatSummary{};

  if (!m.hasTarget) return false;
  if (m.targetKind != targetKind) return false;
  if (m.targetId != targetId) return false;
  if (m.ttlSimSec <= 0.0) return false;

  const double minCos = std::clamp(params.minApproachCos, -1.0, 1.0);
  const double minClosing = std::max(0.0, params.minClosingKmS);
  const double maxDist = std::max(0.0, params.maxConsiderDistKm);

  const math::Vec3d toTarget = targetPosKm - m.posKm;
  const double distKm = toTarget.length();
  if (!(distKm >= 0.0) || !std::isfinite(distKm)) return false;
  if (distKm > maxDist) return false;

  const math::Vec3d toDir = math::safeNormalized(toTarget, math::Vec3d{0, 0, 1}, 1e-12);
  const math::Vec3d mvDir = math::safeNormalized(m.velKmS, math::Vec3d{0, 0, 1}, 1e-12);
  const double approachCos = math::dot(mvDir, toDir);
  if (approachCos < minCos) return false;

  const math::Vec3d relVel = m.velKmS - targetVelKmS;
  const double closingKmS = math::dot(relVel, toDir);
  if (!(closingKmS > minClosing)) return false;

  double ttiSec = (distKm > 1e-9) ? (distKm / std::max(closingKmS, 1e-9)) : 0.0;
  if (!std::isfinite(ttiSec) || ttiSec < 0.0) ttiSec = 0.0;

  // Clamp to missile TTL when provided.
  if (m.ttlSimSec > 0.0) {
    ttiSec = std::min(ttiSec, std::max(0.0, m.ttlSimSec));
  }

  out->inbound = true;
  out->seeker = m.seeker;
  out->distKm = distKm;
  out->closingKmS = closingKmS;
  out->ttiSec = ttiSec;
  out->approachCos = approachCos;
  out->missileIndex = missileIndex;
  out->shooterId = m.shooterId;
  out->fromPlayer = m.fromPlayer;
  return true;
}

static void advanceSphereTargets(std::vector<SphereTarget>& targets, double dtSim) {
  if (!(dtSim > 1e-9) || !std::isfinite(dtSim)) return;

  for (auto& t : targets) {
    t.centerKm = t.centerKm + t.velKmS * dtSim + t.accelKmS2 * (0.5 * dtSim * dtSim);
    t.velKmS = t.velKmS + t.accelKmS2 * dtSim;
  }
}

PointDefenseEngagementPlan planPointDefenseEngagement(const Missile* missiles,
                                                      std::size_t missileCount,
                                                      CombatTargetKind targetKind,
                                                      core::u64 targetId,
                                                      const math::Vec3d& targetPosKm,
                                                      const math::Vec3d& targetVelKmS,
                                                      const math::Vec3d& turretAimDirWorld,
                                                      const SphereTarget* worldTargets,
                                                      std::size_t worldTargetCount,
                                                      core::u64 seed,
                                                      const PointDefenseEngagementParams& params) {
  PointDefenseEngagementPlan out{};

  if (!missiles || missileCount == 0) return out;

  const int maxShots = std::max(0, params.maxShots);
  if (maxShots <= 0) return out;

  // Total planning horizon.
  double horizon = params.maxPlanSimSec;
  if (!(horizon > 1e-9) || !std::isfinite(horizon)) {
    horizon = params.pointDefense.maxPlanSimSec;
  }
  if (!(horizon > 1e-9) || !std::isfinite(horizon)) {
    horizon = 6.0;
  }
  horizon = std::max(0.0, horizon);

  const double cooldown = std::max(0.0, params.shotCooldownSec);
  const double margin = std::max(0.0, params.interceptBeforeImpactMarginSec);

  // Copy candidate missiles (so we can safely advance simulation time).
  std::vector<Missile> simMissiles;
  std::vector<std::size_t> indexMap;
  simMissiles.reserve(missileCount);
  indexMap.reserve(missileCount);

  for (std::size_t i = 0; i < missileCount; ++i) {
    const Missile& m = missiles[i];
    if (!m.hasTarget) continue;
    if (m.targetKind != targetKind) continue;
    if (m.targetId != targetId) continue;
    if (m.ttlSimSec <= 0.0) continue;

    simMissiles.push_back(m);
    indexMap.push_back(i);
  }

  if (simMissiles.empty()) return out;

  // Copy world targets so we can advance them alongside missiles.
  std::vector<SphereTarget> simTargets;
  if (worldTargets && worldTargetCount > 0) {
    simTargets.assign(worldTargets, worldTargets + worldTargetCount);
  }

  // If the defended target is present in the target array, use its (advanced)
  // kinematics for threat estimation + fire-control. This keeps the engagement
  // schedule consistent with the same target list used by the main sim.
  bool haveDefenderTarget = false;
  std::size_t defenderTargetIndex = 0;
  for (std::size_t ti = 0; ti < simTargets.size(); ++ti) {
    if (simTargets[ti].kind == targetKind && simTargets[ti].id == targetId) {
      haveDefenderTarget = true;
      defenderTargetIndex = ti;
      break;
    }
  }

  math::Vec3d turretDir = turretAimDirWorld;
  if (turretDir.lengthSq() < 1e-12) turretDir = {0, 0, 1};
  turretDir = turretDir.normalized();

  std::vector<MissileDetonation> dets;
  std::vector<MissileHit> hits;

  double tNow = 0.0;

  while ((int)out.shots.size() < maxShots && tNow <= horizon + 1e-9 && !simMissiles.empty()) {
    // Target (defended platform) kinematics at the current planning time.
    math::Vec3d targetPosNow = targetPosKm + targetVelKmS * tNow;
    math::Vec3d targetVelNow = targetVelKmS;
    if (haveDefenderTarget && defenderTargetIndex < simTargets.size()) {
      targetPosNow = simTargets[defenderTargetIndex].centerKm;
      targetVelNow = simTargets[defenderTargetIndex].velKmS;
    }

    // Build a list of inbound threats at this time.
    struct Entry {
      std::size_t simIndex;
      MissileThreatSummary threat;
    };

    std::vector<Entry> threats;
    threats.reserve(simMissiles.size());

    for (std::size_t si = 0; si < simMissiles.size(); ++si) {
      MissileThreatSummary s;
      if (computeInboundThreatSummary(simMissiles[si],
                                      indexMap[si],
                                      targetKind,
                                      targetId,
                                      targetPosNow,
                                      targetVelNow,
                                      params.threat,
                                      &s)) {
        threats.push_back({si, s});
      }
    }

    if (threats.empty()) break;

    std::sort(threats.begin(), threats.end(), [](const Entry& a, const Entry& b) {
      return a.threat.ttiSec < b.threat.ttiSec;
    });

    // Clamp per-shot planning horizon to the remaining schedule horizon.
    const double remainHorizon = std::max(0.0, horizon - tNow);
    if (!(remainHorizon > 1e-9)) break;

    bool scheduled = false;

    for (const Entry& e : threats) {
      if (e.simIndex >= simMissiles.size()) continue;
      const Missile& m = simMissiles[e.simIndex];

      // Configure the per-shot solver.
      PointDefenseParams pd = params.pointDefense;
      if (pd.maxPlanSimSec > 0.0) pd.maxPlanSimSec = std::min(pd.maxPlanSimSec, remainHorizon);
      else pd.maxPlanSimSec = remainHorizon;

      if (params.clampShotHorizonToThreatTti && std::isfinite(e.threat.ttiSec) && e.threat.ttiSec > 0.0) {
        pd.maxPlanSimSec = std::min(pd.maxPlanSimSec, std::max(0.0, e.threat.ttiSec));
      }

      const SphereTarget* tgtPtr = simTargets.empty() ? nullptr : simTargets.data();
      const std::size_t tgtCount = simTargets.size();

      // Seed is mixed with the original missile index and shot index to keep ordering stable.
      const core::u64 shotSeed = seed ^ (core::u64)indexMap[e.simIndex] ^ (core::u64)out.shots.size() * 0x9E3779B97F4A7C15ull;

      PointDefenseSolution sol = planPointDefenseShot(m,
                                                    /*shooterPosKm=*/targetPosNow,
                                                    /*shooterVelKmS=*/targetVelNow,
                                                    turretDir,
                                                    tgtPtr,
                                                    tgtCount,
                                                    shotSeed,
                                                    pd);

      // Policy gate: require intercept before estimated impact.
      if (sol.valid && std::isfinite(e.threat.ttiSec)) {
        const double latest = std::max(0.0, e.threat.ttiSec - margin);
        if (sol.interceptTimeSec > latest + 1e-9) {
          sol = PointDefenseSolution{};
        }
      }

      if (!sol.valid) continue;

      // Convert to absolute times within this plan.
      PointDefenseEngagementShot shot{};
      shot.valid = true;
      shot.threat = e.threat;
      shot.solution = sol;
      shot.solution.fireDelaySec = tNow + sol.fireDelaySec;
      shot.solution.interceptTimeSec = tNow + sol.interceptTimeSec;

      out.shots.push_back(shot);
      out.valid = true;

      // Update turret state: assume it is aligned to the shot aim direction after firing.
      if (sol.aimDirWorld.lengthSq() > 1e-18) {
        turretDir = sol.aimDirWorld.normalized();
      }

      // Remove the engaged missile (assume the shot succeeds).
      simMissiles.erase(simMissiles.begin() + (std::ptrdiff_t)e.simIndex);
      indexMap.erase(indexMap.begin() + (std::ptrdiff_t)e.simIndex);

      // Advance local simulation time to the next shot opportunity.
      double dtAdvance = sol.fireDelaySec + cooldown;

      // Ensure forward progress even with instant-fire, zero-cooldown setups.
      double dtMin = pd.simStepSec;
      if (!(dtMin > 1e-6) || !std::isfinite(dtMin)) dtMin = 0.02;
      dtMin = std::clamp(dtMin, 1e-4, 0.25);
      if (dtAdvance < dtMin) dtAdvance = dtMin;

      if (dtAdvance > 1e-9 && !simMissiles.empty()) {
        dets.clear();
        hits.clear();
        stepMissiles(simMissiles, dtAdvance, tgtPtr, tgtCount, dets, hits);
      }

      if (dtAdvance > 1e-9 && !simTargets.empty()) {
        advanceSphereTargets(simTargets, dtAdvance);
      }

      tNow += std::max(0.0, dtAdvance);

      // Cull missiles that detonated / expired during the advance.
      if (!simMissiles.empty()) {
        std::size_t w = 0;
        for (std::size_t r = 0; r < simMissiles.size(); ++r) {
          if (simMissiles[r].ttlSimSec > 0.0) {
            if (w != r) {
              simMissiles[w] = simMissiles[r];
              indexMap[w] = indexMap[r];
            }
            ++w;
          }
        }
        simMissiles.resize(w);
        indexMap.resize(w);
      }

      scheduled = true;
      break;
    }

    if (!scheduled) break;
  }

  return out;
}

PointDefenseSolution planPointDefenseShot(const Missile& missile,
                                          const math::Vec3d& shooterPosKm,
                                          const math::Vec3d& shooterVelKmS,
                                          const math::Vec3d& turretAimDirWorld,
                                          const SphereTarget* worldTargets,
                                          std::size_t worldTargetCount,
                                          core::u64 seed,
                                          const PointDefenseParams& params) {
  (void)seed;

  PointDefenseSolution best{};

  const double muzzle = std::max(0.0, params.muzzleSpeedKmS);
  if (!(muzzle > 1e-9) || !std::isfinite(muzzle)) {
    return best;
  }

  // Planning horizon (clamped by missile TTL when provided).
  double horizon = std::max(0.0, params.maxPlanSimSec);
  if (missile.ttlSimSec > 0.0) {
    horizon = std::min(horizon, std::max(0.0, missile.ttlSimSec));
  }
  if (!(horizon > 1e-6) || !std::isfinite(horizon)) {
    return best;
  }

  // Prediction step.
  double dt = params.simStepSec;
  if (!(dt > 1e-6) || !std::isfinite(dt)) dt = 0.02;
  dt = std::clamp(dt, 1e-4, 0.25);

  // Turret mount axis + current aim direction.
  math::Vec3d mountAxis = params.mountAxisWorld;
  if (mountAxis.lengthSq() < 1e-12) mountAxis = {0, 0, 1};
  mountAxis = mountAxis.normalized();

  math::Vec3d aimNow = turretAimDirWorld;
  if (aimNow.lengthSq() < 1e-12) aimNow = mountAxis;
  aimNow = aimNow.normalized();

  const double mountCos = std::clamp(params.mountConeCos, -1.0, 1.0);
  const double slewRate = std::max(0.0, params.turretSlewRateRadS);
  const bool useSlew = (slewRate > 1e-9);

  auto isBetter = [&](double fireDelay, double interceptTime) -> bool {
    if (!best.valid) return true;

    if (params.preferEarliestIntercept) {
      if (interceptTime + 1e-9 < best.interceptTimeSec) return true;
      if (std::fabs(interceptTime - best.interceptTimeSec) <= 1e-9 && fireDelay + 1e-9 < best.fireDelaySec) return true;
      return false;
    }

    // Prefer earliest fire.
    if (fireDelay + 1e-9 < best.fireDelaySec) return true;
    if (std::fabs(fireDelay - best.fireDelaySec) <= 1e-9 && interceptTime + 1e-9 < best.interceptTimeSec) return true;
    return false;
  };

  auto consider = [&](double tSec, const math::Vec3d& missilePosKm) {
    if (!(tSec > 1e-9) || !std::isfinite(tSec)) return;

    const math::Vec3d shooterAtT = shooterPosKm + shooterVelKmS * tSec;
    const math::Vec3d rel = missilePosKm - shooterAtT;
    const double dist = rel.length();
    if (!(dist >= 0.0) || !std::isfinite(dist)) return;

    const double travel = dist / muzzle;
    if (!std::isfinite(travel)) return;

    const double fireDelayRaw = tSec - travel;
    if (fireDelayRaw < -1e-9) return;

    const double fireDelay = std::max(0.0, fireDelayRaw);

    math::Vec3d aim = (dist > 1e-9) ? (rel / dist) : aimNow;
    if (aim.lengthSq() < 1e-12) return;
    aim = aim.normalized();

    // Mount cone check.
    if (mountCos > -0.5) {
      if (math::dot(mountAxis, aim) < mountCos) return;
    }

    // Slew check.
    const double d = std::clamp(math::dot(aimNow, aim), -1.0, 1.0);
    const double ang = std::acos(d);
    double requiredSlew = 0.0;
    if (useSlew) {
      requiredSlew = ang / slewRate;
      if (fireDelay + 1e-9 < requiredSlew) return;
    }

    // Optional LOS check (asteroids only).
    if (params.requireLineOfSight && worldTargets && worldTargetCount > 0) {
      const math::Vec3d shooterAtFire = shooterPosKm + shooterVelKmS * fireDelay;
      if (occludedByAsteroidSpheres(shooterAtFire,
                                    missilePosKm,
                                    worldTargets,
                                    worldTargetCount,
                                    params.lineOfSightOcclusionPadKm,
                                    /*ignoreAsteroidId=*/0)) {
        return;
      }
    }

    if (!isBetter(fireDelay, tSec)) return;

    best.valid = true;
    best.fireDelaySec = fireDelay;
    best.interceptTimeSec = tSec;
    best.aimDirWorld = aim;
    best.interceptPointKm = missilePosKm;
    best.aimAngleRad = ang;
    best.requiredSlewTimeSec = requiredSlew;
    best.rangeAtInterceptKm = dist;
  };

  // --- Predict missile motion and search for a feasible intercept window. ---

  if (params.useFullMissileSimulation) {
    std::vector<Missile> sim;
    sim.reserve(1);
    sim.push_back(missile);

    std::vector<MissileDetonation> dets;
    std::vector<MissileHit> hits;

    double t = 0.0;
    while (t <= horizon + 1e-9 && !sim.empty()) {
      consider(t, sim[0].posKm);

      const double step = std::min(dt, horizon - t);
      if (!(step > 1e-9)) break;

      stepMissiles(sim, step, worldTargets, worldTargetCount, dets, hits);
      t += step;
    }
  } else {
    // Inertial (constant-velocity) fallback with a minimal motor model.
    Missile m = missile;
    double t = 0.0;

    while (t <= horizon + 1e-9 && m.ttlSimSec > 0.0) {
      consider(t, m.posKm);

      const double step = std::min(dt, horizon - t);
      if (!(step > 1e-9)) break;

      // Simple boost/coast speed-up along current velocity.
      math::Vec3d v = m.velKmS;
      double speed = v.length();
      if (speed < 1e-9) {
        speed = 0.0;
        v = {0, 0, 1};
      }

      const double accel = std::max(0.0, m.thrustAccelKmS2);
      double burn = std::max(0.0, m.motorBurnRemainingSimSec);
      const double maxSpeed = std::max(0.0, m.maxSpeedKmS);

      if (accel > 0.0 && burn > 0.0 && maxSpeed > speed + 1e-12 && speed > 1e-9) {
        const double burnDt = std::min(step, burn);
        const double newSpeed = std::min(maxSpeed, speed + accel * burnDt);
        const math::Vec3d dir = v.normalized();
        m.velKmS = dir * newSpeed;
        m.motorBurnRemainingSimSec = burn - burnDt;
      }

      m.posKm = m.posKm + m.velKmS * step;
      m.ttlSimSec = std::max(0.0, m.ttlSimSec - step);
      m.ageSimSec = std::max(0.0, m.ageSimSec + step);

      t += step;
    }
  }

  if (best.valid && best.aimDirWorld.lengthSq() > 1e-18) {
    best.aimDirWorld = best.aimDirWorld.normalized();
  }

  return best;
}

}  // namespace stellar::sim
