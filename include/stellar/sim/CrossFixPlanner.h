#pragma once

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/BearingTrack.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// CrossFixPlanner — helper for bearing-only "cross-fix" workflows
// -----------------------------------------------------------------------------
//
// Many sensing/gameplay situations are bearing-only:
//   - weak/ghost radar contacts under jamming
//   - passive sensor bearings with no range
//   - mission clues where the player must triangulate a signal
//
// A single bearing defines a *ray*; to triangulate a unique point, the observer
// must move to create a baseline. This header provides a deterministic helper
// that suggests a lateral waypoint for creating that baseline.
//
// Design goals:
//   - deterministic + headless (safe for tests/tools/server sims)
//   - minimal state: caller supplies an optional seed to keep a stable side
//   - does not require the target to be known (uses BearingTrack estimate if any)

struct CrossFixParams {
  // --- Baseline selection -----------------------------------------------------
  // If the bearing track already has a fix (initialized), and sigmaKm is finite,
  // we can choose a baseline scaled by the uncertainty.
  bool baselineFromSigma{true};
  double baselineFromSigmaMult{2.0};

  // If no fix exists yet, callers often have an approximate range guess (e.g.
  // radar strength inversion). When enabled, baseline = rangeGuess * frac.
  bool baselineFromRange{true};
  double baselineFromRangeFrac{0.12};

  // Fallback baseline if neither sigma nor range guess are available.
  double baselineKm{8000.0};

  // Clamp recommended baselines.
  double minBaselineKm{500.0};
  double maxBaselineKm{250000.0};

  // Optional clamp for the waypoint distance from the observer.
  // 0 => no additional clamp.
  double maxWaypointDistKm{0.0};

  // --- Frame selection --------------------------------------------------------
  // Preferred world "up" axis for selecting a tangent direction.
  // This is only used to pick a stable perpendicular direction; the suggested
  // waypoint is always lateral to the bearing.
  math::Vec3d worldUp{0, 1, 0};

  // --- Side selection ---------------------------------------------------------
  // Preferred side sign in {-1, +1}. 0 => pick deterministically from seed.
  int preferredSideSign{0};

  // Optional seed used to break symmetry when preferredSideSign==0.
  // In a game context, pass a stable contact/signal id so the suggestion
  // doesn't flip every frame.
  core::u64 seed{0};

  // --- Predicted information gain (optional) ---------------------------------
  //
  // These options let tools/UI ask: "If I take ONE more bearing measurement from
  // the suggested waypoint, how much closer am I to a solvable cross-fix?"
  //
  // When enabled, the planner predicts the next bearing direction from each
  // candidate waypoint to the current aim point, updates the bearing track's
  // normal matrix A (information matrix) with that hypothetical measurement,
  // and evaluates solvability progress using BearingTrackSolveDiagnostics.

  // If true, populate CrossFixCandidate::predicted* fields.
  // (Automatically enabled when chooseSideByPrediction or optimizeBaselineByPrediction is on.)
  bool includePrediction{false};

  // If true, pick the preferred side (+/-) by predicted progress instead of
  // by seed/preferredSideSign.
  bool chooseSideByPrediction{false};

  // If true, search a set of baseline multipliers around the initially computed
  // baseline and pick the baseline/side that maximizes predicted progress.
  // (Implies chooseSideByPrediction.)
  bool optimizeBaselineByPrediction{false};

  // Assumed measurement weight for the *next* bearing sample used in prediction.
  // If <= 0, prediction is disabled.
  double predictedMeasurementWeight{1.0};

  // Baseline search configuration (only used when optimizeBaselineByPrediction is true).
  //
  // The search will evaluate up to `baselineSearchSteps` multipliers in the range
  // [baselineSearchMinMult, baselineSearchMaxMult], plus always considering 1.0.
  int baselineSearchSteps{7};
  double baselineSearchMinMult{0.5};
  double baselineSearchMaxMult{3.0};

  // Solvability gate parameters used for prediction (determinantEps, minEffectiveWeight).
  // Other fields are ignored.
  BearingTrackParams predictedSolveParams{};
};

struct CrossFixCandidate {
  // Suggested waypoint in world space.
  math::Vec3d waypointKm{0, 0, 0};

  // Unit direction from observer to waypoint.
  math::Vec3d moveDirUnit{0, 0, 0};

  // +1 or -1 relative to the planner's stable basis.
  int sideSign{0};

  // --- Optional prediction (when enabled) ------------------------------------
  bool hasPrediction{false};

  // Predicted solvability gates after taking ONE more measurement from this waypoint.
  BearingTrackSolveDiagnostics predictedSolve{};

  // Eigenvalues (ascending) of the predicted information matrix A' (for conditioning intuition).
  math::Vec3d predictedInfoEigenvalues{0, 0, 0};

  // Condition number estimate of A' (lambda_max / max(lambda_min, eps)).
  double predictedConditionNumber{0.0};
};

struct CrossFixRecommendation {
  bool valid{false};

  // True when the provided bearing track already has a triangulated position.
  bool trackHasFix{false};

  // Baseline used to place waypoints (km).
  double baselineKm{0.0};

  // Distance from observer to the aim point (estimated target location) if known.
  // 0 when unknown.
  double rangeToAimKm{0.0};

  // Expected change in bearing angle when moving to `preferred` (degrees).
  // 0 when unknown.
  double expectedBearingChangeDeg{0.0};

  // True when predicted metrics were computed.
  bool hasPrediction{false};

  // True when the planner used predicted metrics to choose side/baseline.
  bool usedPredictionForChoice{false};

  // The planner chooses a deterministic perpendicular direction and returns
  // a preferred + alternate option (left/right around the bearing).
  CrossFixCandidate preferred{};
  CrossFixCandidate alternate{};

  // Debug: the point the planner treated as the target aim point.
  // If the bearing track has a fix, this is track.posKm. Otherwise it may be
  // (observer + bearingDirUnit * rangeGuess) if rangeGuess was provided.
  math::Vec3d aimPointKm{0, 0, 0};

  // Debug: the normalized bearing direction used by the planner.
  math::Vec3d bearingDirUnit{0, 0, 0};
};

inline int clampSideSign(int s) {
  if (s > 0) return 1;
  if (s < 0) return -1;
  return 0;
}

inline core::u64 quantizeUnitComponentToU64(double x) {
  if (!std::isfinite(x)) x = 0.0;
  x = std::clamp(x, -1.0, 1.0);
  // Quantize to ~20 bits.
  const long long q = (long long)std::llround(x * 1048576.0);
  // Bias to be non-negative.
  return (core::u64)(q + 1048576ll);
}

inline math::Vec3d stablePerpUnit(const math::Vec3d& dirUnit, const math::Vec3d& preferredUp) {
  if (!math::isFinite(dirUnit)) return {0, 0, 0};

  math::Vec3d u = preferredUp;
  if (!math::isFinite(u) || u.lengthSq() < 1.0e-18) u = {0, 1, 0};
  u = math::safeNormalized(u, {0, 1, 0}, 1e-18);

  const auto chooseAlt = [&](const math::Vec3d& alt) {
    u = math::safeNormalized(alt, {0, 1, 0}, 1e-18);
  };

  const double c0 = std::abs(math::dot(dirUnit, u));
  if (c0 > 0.95) {
    chooseAlt({1, 0, 0});
    const double c1 = std::abs(math::dot(dirUnit, u));
    if (c1 > 0.95) {
      chooseAlt({0, 0, 1});
    }
  }

  // Perpendicular (sign depends on order; callers can flip).
  math::Vec3d perp = math::cross(dirUnit, u);
  perp = math::safeNormalized(perp, {0, 0, 0}, 1e-18);

  // If still degenerate (should be rare), fall back to a fixed axis.
  if (perp.lengthSq() <= 1.0e-18) {
    perp = math::cross(dirUnit, math::Vec3d{0, 0, 1});
    perp = math::safeNormalized(perp, {0, 0, 0}, 1e-18);
  }

  return perp;
}

inline double clampBaselineKm(double baselineKm, const CrossFixParams& params) {
  double b = baselineKm;
  if (!std::isfinite(b)) b = params.baselineKm;

  const double bMin = (std::isfinite(params.minBaselineKm) && params.minBaselineKm > 0.0)
                        ? params.minBaselineKm
                        : 0.0;
  const double bMax = params.maxBaselineKm;

  if (std::isfinite(bMax) && bMax > 0.0) {
    b = std::clamp(b, bMin, bMax);
  } else {
    b = std::max(b, bMin);
  }

  if (std::isfinite(params.maxWaypointDistKm) && params.maxWaypointDistKm > 0.0) {
    b = std::min(b, params.maxWaypointDistKm);
  }

  return b;
}

// Score used for predicted information-gain comparisons.
//
// - If the predicted geometry is not yet solvable, we prioritize progressing toward
//   solvability (0..1).
// - Once solvable, we continue to prefer geometries with more determinant "margin"
//   (loosely: better conditioning / less numerical fragility).
inline double predictedScore(const BearingTrackSolveDiagnostics& d) {
  if (!d.canSolve) return d.progress01;

  // Use determinant margin as a simple scalar that keeps increasing beyond the gate.
  const double denom = (d.detThresh > 1.0e-30) ? d.detThresh : 1.0;
  const double margin = d.detAbs / denom;

  // Log keeps values tame but still monotonic.
  return 1.0 + std::log(1.0 + std::max(0.0, margin));
}

// Recommend a cross-fix waypoint given an observer position and a bearing.
//
// track:
//   Optional. If provided and initialized, its position/sigma are used to scale
//   the baseline and compute expected bearing change.
//
// rangeGuessKm:
//   Optional approximate range to the target when no fix exists yet.
//   Pass <= 0 to disable.
inline CrossFixRecommendation recommendCrossFixWaypoint(const BearingTrack3d* track,
                                                        const math::Vec3d& observerPosKm,
                                                        const math::Vec3d& bearingDirWorld,
                                                        double rangeGuessKm,
                                                        const CrossFixParams& params = {}) {
  CrossFixRecommendation out{};

  if (!math::isFinite(observerPosKm) || !math::isFinite(bearingDirWorld)) return out;

  out.bearingDirUnit = math::safeNormalized(bearingDirWorld, {0, 0, 0}, 1e-18);
  if (out.bearingDirUnit.lengthSq() <= 1.0e-18) return out;

  // Determine if we have a triangulated aim point.
  out.trackHasFix = (track && track->initialized && math::isFinite(track->posKm));

  if (out.trackHasFix) {
    out.aimPointKm = track->posKm;
  } else if (std::isfinite(rangeGuessKm) && rangeGuessKm > 0.0) {
    out.aimPointKm = observerPosKm + out.bearingDirUnit * rangeGuessKm;
  } else {
    // Unknown range: pick a dummy point along the bearing for direction selection.
    const double d = std::max(1.0, params.baselineKm * 5.0);
    out.aimPointKm = observerPosKm + out.bearingDirUnit * d;
  }

  // Direction from observer to the aim point.
  const math::Vec3d relAim = out.aimPointKm - observerPosKm;
  const double relAimLenSq = relAim.lengthSq();
  const double relAimLen = (relAimLenSq > 0.0) ? std::sqrt(relAimLenSq) : 0.0;
  if (std::isfinite(relAimLen) && relAimLen > 1.0e-9) {
    out.rangeToAimKm = relAimLen;
  }

  const math::Vec3d rUnit = (relAimLenSq > 1.0e-18)
                              ? math::safeNormalized(relAim, out.bearingDirUnit, 1e-18)
                              : out.bearingDirUnit;

  math::Vec3d tangentUnit = stablePerpUnit(rUnit, params.worldUp);
  if (tangentUnit.lengthSq() <= 1.0e-18) return out;

  // --- Baseline selection -----------------------------------------------------
  double baseline0 = params.baselineKm;

  if (params.baselineFromSigma && out.trackHasFix) {
    const double sig = (track && std::isfinite(track->sigmaKm)) ? track->sigmaKm : 0.0;
    if (sig > 0.0 && std::isfinite(params.baselineFromSigmaMult) && params.baselineFromSigmaMult > 0.0) {
      baseline0 = sig * params.baselineFromSigmaMult;
    }
  } else if (params.baselineFromRange && std::isfinite(rangeGuessKm) && rangeGuessKm > 0.0) {
    if (std::isfinite(params.baselineFromRangeFrac) && params.baselineFromRangeFrac > 0.0) {
      baseline0 = rangeGuessKm * params.baselineFromRangeFrac;
    }
  }

  baseline0 = clampBaselineKm(baseline0, params);

  // --- Deterministic side selection (fallback/tie-break) ---------------------
  int deterministicSide = clampSideSign(params.preferredSideSign);

  if (deterministicSide == 0) {
    core::u64 h = params.seed;
    h = core::hashCombine(h, core::fnv1a64("CrossFixPlanner.side"));

    // Incorporate direction so distinct bearings don't all pick the same side when seed=0.
    h = core::hashCombine(h, quantizeUnitComponentToU64(rUnit.x));
    h = core::hashCombine(h, quantizeUnitComponentToU64(rUnit.y));
    h = core::hashCombine(h, quantizeUnitComponentToU64(rUnit.z));

    core::SplitMix64 rng(h);
    deterministicSide = rng.chance(0.5) ? 1 : -1;
  }

  // --- Prediction configuration ----------------------------------------------
  const bool wantPrediction = params.includePrediction || params.chooseSideByPrediction || params.optimizeBaselineByPrediction;

  double wPred = params.predictedMeasurementWeight;
  if (!std::isfinite(wPred) || wPred < 0.0) wPred = 0.0;

  const bool canPredict = wantPrediction && track && track->A.isFinite() && std::isfinite(track->weight) && (wPred > 0.0);

  const auto makeCandidates = [&](double baselineKm, CrossFixCandidate& plus, CrossFixCandidate& minus) {
    plus = CrossFixCandidate{};
    plus.sideSign = 1;
    plus.moveDirUnit = tangentUnit;
    plus.waypointKm = observerPosKm + tangentUnit * baselineKm;

    minus = CrossFixCandidate{};
    minus.sideSign = -1;
    minus.moveDirUnit = tangentUnit * -1.0;
    minus.waypointKm = observerPosKm - tangentUnit * baselineKm;
  };

  const auto predict = [&](CrossFixCandidate& c) {
    if (!canPredict) return;

    const math::Vec3d rel = out.aimPointKm - c.waypointKm;
    const math::Vec3d dUnit = math::safeNormalized(rel, {0, 0, 0}, 1e-18);
    if (dUnit.lengthSq() <= 1.0e-18) return;

    const math::Mat3d ddT = math::Mat3d::outerProduct(dUnit, dUnit);
    const math::Mat3d P = math::Mat3d::identity() - ddT;

    // Predicted information matrix after one more measurement.
    const math::Mat3d Ap = track->A + P * wPred;
    const double wp = track->weight + wPred;

    BearingTrack3d tmp{};
    tmp.A = Ap;
    tmp.weight = wp;

    c.predictedSolve = bearingTrackSolveDiagnostics(tmp, params.predictedSolveParams);

    const auto eig = math::symmetricEigenDecomposition3x3(Ap);
    if (eig.valid) {
      c.predictedInfoEigenvalues = eig.eigenvalues;

      const double l0 = std::max(1.0e-18, std::abs(eig.eigenvalues.x));
      const double l2 = std::max(1.0e-18, std::abs(eig.eigenvalues.z));
      c.predictedConditionNumber = l2 / l0;
      if (!std::isfinite(c.predictedConditionNumber)) c.predictedConditionNumber = 0.0;
    }

    c.hasPrediction = true;
  };

  // --- Baseline search (optional) --------------------------------------------
  // Always consider baseline0 first so we preserve legacy behavior when
  // prediction doesn't provide a strictly better option.
  static constexpr int kMaxBaselines = 16;
  std::array<double, kMaxBaselines> baselines{};
  int baselineCount = 0;

  baselines[baselineCount++] = baseline0;

  if (params.optimizeBaselineByPrediction && canPredict) {
    int steps = params.baselineSearchSteps;
    if (!std::isfinite((double)steps) || steps < 0) steps = 0;
    steps = std::clamp(steps, 0, kMaxBaselines);

    double mMin = params.baselineSearchMinMult;
    double mMax = params.baselineSearchMaxMult;
    if (!std::isfinite(mMin) || mMin <= 0.0) mMin = 0.5;
    if (!std::isfinite(mMax) || mMax <= 0.0) mMax = 3.0;
    if (mMax < mMin) std::swap(mMin, mMax);

    const auto addBaseline = [&](double b) {
      if (!(std::isfinite(b) && b > 0.0)) return;
      b = clampBaselineKm(b, params);

      for (int i = 0; i < baselineCount; ++i) {
        if (std::abs(baselines[i] - b) <= 1.0e-6) return;
      }

      if (baselineCount < kMaxBaselines) {
        baselines[baselineCount++] = b;
      }
    };

    // Always include 1.0 multiplier explicitly.
    addBaseline(baseline0 * 1.0);

    if (steps >= 2) {
      for (int i = 0; i < steps; ++i) {
        const double t = (steps == 1) ? 0.0 : (double)i / (double)(steps - 1);
        const double mult = lerp(mMin, mMax, t);
        addBaseline(baseline0 * mult);
      }
    }
  }

  // Evaluate candidates and pick baseline/side.
  double bestScore = -1.0e300;
  double bestBaseline = baseline0;
  int bestSide = deterministicSide;
  CrossFixCandidate bestPlus{};
  CrossFixCandidate bestMinus{};

  for (int i = 0; i < baselineCount; ++i) {
    const double b = baselines[i];

    CrossFixCandidate plus{};
    CrossFixCandidate minus{};
    makeCandidates(b, plus, minus);

    if (wantPrediction) {
      predict(plus);
      predict(minus);
    }

    // Choose side for this baseline.
    int localSide = deterministicSide;

    const bool usePredictionForChoice = params.chooseSideByPrediction || params.optimizeBaselineByPrediction;

    if (usePredictionForChoice && plus.hasPrediction && minus.hasPrediction) {
      const double sPlus = predictedScore(plus.predictedSolve);
      const double sMinus = predictedScore(minus.predictedSolve);

      if (sPlus > sMinus + 1.0e-12) {
        localSide = 1;
      } else if (sMinus > sPlus + 1.0e-12) {
        localSide = -1;
      } else {
        // Tie: keep deterministic seed-based side.
        localSide = deterministicSide;
      }
    }

    // Score this baseline by the preferred side's predicted score (if available).
    double localScore = 0.0;
    if (params.optimizeBaselineByPrediction && plus.hasPrediction && minus.hasPrediction) {
      // Optimization mode: score by the better of the two sides.
      const double sPlus = predictedScore(plus.predictedSolve);
      const double sMinus = predictedScore(minus.predictedSolve);
      localScore = std::max(sPlus, sMinus);

      // Use the maximizing side (tie -> deterministic).
      if (sPlus > sMinus + 1.0e-12) {
        localSide = 1;
      } else if (sMinus > sPlus + 1.0e-12) {
        localSide = -1;
      } else {
        localSide = deterministicSide;
      }
    } else if (usePredictionForChoice && plus.hasPrediction && minus.hasPrediction) {
      // Not optimizing baseline: still allow predicted choice to influence the score.
      localScore = (localSide > 0) ? predictedScore(plus.predictedSolve) : predictedScore(minus.predictedSolve);
    }

    // Baseline optimization: update if strictly better.
    if (params.optimizeBaselineByPrediction && canPredict) {
      if (localScore > bestScore + 1.0e-12) {
        bestScore = localScore;
        bestBaseline = b;
        bestSide = localSide;
        bestPlus = plus;
        bestMinus = minus;
      }
    } else {
      // No baseline optimization: keep baseline0 and break.
      bestBaseline = b;
      bestSide = localSide;
      bestPlus = plus;
      bestMinus = minus;
      break;
    }
  }

  out.baselineKm = bestBaseline;

  // Side selection result.
  out.usedPredictionForChoice = (canPredict && (params.chooseSideByPrediction || params.optimizeBaselineByPrediction));

  out.preferred = (bestSide > 0) ? bestPlus : bestMinus;
  out.alternate = (bestSide > 0) ? bestMinus : bestPlus;

  out.hasPrediction = out.preferred.hasPrediction || out.alternate.hasPrediction;

  // --- Expected bearing change ------------------------------------------------
  // Only meaningful if we have an aim point at a sensible distance.
  if (out.rangeToAimKm > 1.0e-6 && std::isfinite(out.rangeToAimKm)) {
    const math::Vec3d newRel = out.aimPointKm - out.preferred.waypointKm;
    const math::Vec3d newUnit = math::safeNormalized(newRel, {0, 0, 0}, 1e-18);
    if (newUnit.lengthSq() > 1.0e-18) {
      double c = math::dot(rUnit, newUnit);
      c = std::clamp(c, -1.0, 1.0);
      const double angRad = std::acos(c);
      if (std::isfinite(angRad)) {
        out.expectedBearingChangeDeg = angRad * (180.0 / 3.14159265358979323846);
      }
    }
  }

  out.valid = true;
  return out;
}

// Convenience overload when a track reference is available.
inline CrossFixRecommendation recommendCrossFixWaypoint(const BearingTrack3d& track,
                                                        const math::Vec3d& observerPosKm,
                                                        const math::Vec3d& bearingDirWorld,
                                                        double rangeGuessKm,
                                                        const CrossFixParams& params = {}) {
  return recommendCrossFixWaypoint(&track, observerPosKm, bearingDirWorld, rangeGuessKm, params);
}

} // namespace stellar::sim
