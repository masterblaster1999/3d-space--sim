#pragma once

#include "stellar/math/Mat3.h"
#include "stellar/math/SymmetricEigen3.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// BearingTrack — bearing-only 3D triangulation (headless, deterministic)
// -----------------------------------------------------------------------------
//
// Motivation:
//   A lot of interesting "space sim" sensing situations are bearing-only:
//     - passive EM/radiation bearings (no range)
//     - intermittent contacts (range drops out under jamming/occlusion)
//     - mission clues where the player must "cross-fix" a signal
//
//   Given multiple observer positions {o_i} and corresponding unit bearings
//   {d_i} pointing at the same target, we can estimate the 3D point p that
//   best fits all lines in a least-squares sense.
//
// Math:
//   We minimize sum_i || (I - d_i d_i^T) (p - o_i) ||^2.
//   This yields a symmetric linear system:
//     A p = b
//   where:
//     A = sum_i w_i (I - d_i d_i^T)
//     b = sum_i w_i (I - d_i d_i^T) o_i
//
// Design constraints:
//   - header-only
//   - deterministic
//   - no dynamic allocations
//
// Notes:
//   - With a single bearing, A is rank-2 (unsolvable for unique p). The track
//     remains invalid until enough geometric diversity exists.
//   - We support exponential forgetting (half-life) so the estimate adapts when
//     the target or observer moves.

struct BearingTrackParams {
  // Exponential forgetting for the normal-equation accumulators.
  // <= 0 disables decay (infinite memory).
  double observationHalfLifeSec{6.0};

  // Small diagonal regularization added when solving. This helps avoid numeric
  // blowups when bearings are near-parallel, but *does not* override the
  // determinant gate (we still require the unregularized system to be solvable).
  double solveRegularization{1.0e-6};

  // Relative determinant threshold used to gate solvability.
  // We scale this by (maxAbsElement^3) so it behaves sensibly under different
  // weight magnitudes.
  double determinantEps{1.0e-12};

  // Require at least this total (decayed) weight before attempting a solve.
  // This helps avoid "solving" from a single overweight measurement.
  double minEffectiveWeight{1.8};

  // Velocity estimation smoothing (half-life). <= 0 snaps to measured velocity.
  double velHalfLifeSec{0.65};

  // Clamp the velocity magnitude (safety against pathological solves).
  double maxSpeedKmS{25000.0};

  // Residual-derived 1-sigma clamp.
  double sigmaMinKm{0.0};
  double sigmaMaxKm{800000.0};
};

struct BearingTrack3d {
  // Normal equation accumulators.
  math::Mat3d A{};
  math::Vec3d b{0, 0, 0};
  double weight{0.0};

  // Residual accumulator (fit quality): sigma^2 ~= E[dist(line, p)^2].
  double residualWeight{0.0};
  double residualSqSum{0.0};

  // Current estimate (km / km/s).
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double sigmaKm{0.0};

  // Seconds since last used measurement.
  double ageSinceMeasSec{0.0};

  bool initialized{false};
};

struct BearingTrackResult {
  bool valid{false};
  bool hasMeasurement{false};

  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double sigmaKm{0.0};
  double ageSinceMeasSec{0.0};
};

inline double halfLifeDecay(double halfLifeSec, double dtSec) {
  if (!std::isfinite(dtSec) || dtSec <= 0.0) return 1.0;
  if (!std::isfinite(halfLifeSec) || halfLifeSec <= 0.0) return 1.0;
  return std::pow(0.5, dtSec / halfLifeSec);
}

inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

inline double lerp(double a, double b, double t) {
  return a + (b - a) * t;
}

inline math::Vec3d lerp(const math::Vec3d& a, const math::Vec3d& b, double t) {
  return a * (1.0 - t) + b * t;
}

inline double maxAbsElement(const math::Mat3d& M) {
  double m = 0.0;
  for (double v : M.m) m = std::max(m, std::abs(v));
  return m;
}


// Diagnostics helper for UI/tools: expose the same solvability gates used by
// updateBearingTrack() without mutating the track.
//
// This is useful for cross-fix workflows where the caller wants to show
// progress toward a solvable geometry before a unique position estimate exists.
struct BearingTrackSolveDiagnostics {
  double effectiveWeight{0.0};
  double minEffectiveWeight{0.0};

  // Determinant of the unregularized normal matrix A.
  double detAbs{0.0};
  double detThresh{0.0};
  double scale{0.0};

  bool weightOk{false};
  bool detOk{false};
  bool canSolve{false};

  // Normalized progress metrics in [0,1].
  double weightProgress01{0.0};
  double detProgress01{0.0};
  double progress01{0.0};
};

inline BearingTrackSolveDiagnostics bearingTrackSolveDiagnostics(const BearingTrack3d& track,
                                                                const BearingTrackParams& params = {}) {
  BearingTrackSolveDiagnostics out{};

  out.effectiveWeight = track.weight;
  out.minEffectiveWeight = params.minEffectiveWeight;

  const double wMin = (std::isfinite(params.minEffectiveWeight) && params.minEffectiveWeight > 0.0)
                        ? params.minEffectiveWeight
                        : 0.0;

  out.weightOk = (wMin <= 0.0) ? true : (track.weight >= wMin);
  out.weightProgress01 = (wMin <= 1.0e-12) ? 1.0 : clamp01(track.weight / wMin);

  const double scale = std::max(1.0, maxAbsElement(track.A));
  out.scale = scale;

  double detEps = params.determinantEps;
  if (!std::isfinite(detEps) || detEps < 0.0) detEps = 0.0;

  out.detThresh = detEps * scale * scale * scale;
  out.detAbs = std::abs(track.A.determinant());

  out.detOk = (out.detAbs > out.detThresh);
  if (out.detThresh <= 0.0) {
    // If the threshold is zero, treat any non-zero determinant as "full" progress.
    out.detProgress01 = (out.detAbs > 0.0) ? 1.0 : 0.0;
  } else {
    out.detProgress01 = clamp01(out.detAbs / out.detThresh);
  }

  out.canSolve = out.weightOk && out.detOk;
  out.progress01 = std::min(out.weightProgress01, out.detProgress01);
  return out;
}

// -----------------------------------------------------------------------------
// BearingTrack uncertainty ellipsoid (principal axes)
// -----------------------------------------------------------------------------
//
// This is a lightweight helper for tools/UI that want to visualize "how uncertain"
// a bearing-only fix is (or how directionally ill-conditioned the geometry is).
//
// If we interpret the normal matrix A as an information matrix, then a rough
// position covariance estimate is:
//   Cov ~= sigma^2 * A^{-1}
// where sigma is the (perpendicular) measurement noise scale.
//
// Because A is symmetric, eigenvectors of A are eigenvectors of Cov, and the
// corresponding covariance eigenvalues are sigma^2 / lambda_i.

struct BearingTrackUncertaintyEllipsoidParams {
  // If > 0, use this sigma scale (km).
  // If <= 0, use track.sigmaKm if it is finite and > 0, otherwise fall back to 1.
  double sigmaScaleKm{-1.0};

  // Minimum eigenvalue (information) used when inverting. This prevents
  // infinite axes for rank-deficient geometries (e.g. a single bearing).
  double minInfoEigenvalue{1.0e-12};

  // Clamp output axis lengths (km) for UI stability.
  double maxAxisKm{1.0e9};
};

struct BearingTrackUncertaintyEllipsoid {
  bool valid{false};

  // Principal axes (unit vectors) as columns.
  math::Mat3d axes{math::Mat3d::identity()};

  // Eigenvalues of the information matrix A (ascending).
  math::Vec3d infoEigenvalues{0, 0, 0};

  // Corresponding covariance eigenvalues (km^2) and 1-sigma semi-axis lengths (km).
  // These are aligned with `axes` / `infoEigenvalues` (same column order).
  math::Vec3d varianceAxisKm2{0, 0, 0};
  math::Vec3d sigmaAxisKm{0, 0, 0};

  // Condition number of the information matrix (lambdaMax / lambdaMinClamped).
  double conditionNumber{0.0};

  // Sigma scale actually used (km).
  double sigmaScaleKm{1.0};
};

inline BearingTrackUncertaintyEllipsoid bearingTrackUncertaintyEllipsoid(
  const BearingTrack3d& track,
  const BearingTrackUncertaintyEllipsoidParams& params = {}) {

  BearingTrackUncertaintyEllipsoid out{};
  if (!track.A.isFinite()) return out;

  const math::SymmetricEigen3Result eig = math::symmetricEigenDecomposition3x3(track.A);
  if (!eig.valid) return out;

  double sigma = params.sigmaScaleKm;
  if (!(std::isfinite(sigma) && sigma > 0.0)) {
    sigma = (std::isfinite(track.sigmaKm) && track.sigmaKm > 0.0) ? track.sigmaKm : 1.0;
  }

  double minEig = params.minInfoEigenvalue;
  if (!std::isfinite(minEig) || minEig <= 0.0) minEig = 1.0e-18;

  double maxAxis = params.maxAxisKm;
  if (!std::isfinite(maxAxis) || maxAxis <= 0.0) maxAxis = 1.0e12;

  out.axes = eig.eigenvectors;
  out.infoEigenvalues = eig.eigenvalues;
  out.sigmaScaleKm = sigma;

  const double sigma2 = sigma * sigma;

  auto invClamp = [&](double lambda) {
    double l = lambda;
    if (!std::isfinite(l) || l < minEig) l = minEig;
    return std::max(l, minEig);
  };

  const double l0 = invClamp(out.infoEigenvalues.x);
  const double l1 = invClamp(out.infoEigenvalues.y);
  const double l2 = invClamp(out.infoEigenvalues.z);

  const double v0 = sigma2 / l0;
  const double v1 = sigma2 / l1;
  const double v2 = sigma2 / l2;

  out.varianceAxisKm2 = {
    (std::isfinite(v0) && v0 >= 0.0) ? v0 : (maxAxis * maxAxis),
    (std::isfinite(v1) && v1 >= 0.0) ? v1 : (maxAxis * maxAxis),
    (std::isfinite(v2) && v2 >= 0.0) ? v2 : (maxAxis * maxAxis),
  };

  auto toSigma = [&](double var) {
    double s = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : maxAxis;
    if (!std::isfinite(s)) s = maxAxis;
    return std::clamp(s, 0.0, maxAxis);
  };

  out.sigmaAxisKm = {
    toSigma(out.varianceAxisKm2.x),
    toSigma(out.varianceAxisKm2.y),
    toSigma(out.varianceAxisKm2.z),
  };

  out.conditionNumber = (l0 > 0.0) ? (l2 / l0) : 0.0;
  if (!std::isfinite(out.conditionNumber)) out.conditionNumber = 0.0;

  out.valid = true;
  return out;
}

// Update the bearing-only track.
//
// Inputs:
//  - dtSec: time step (seconds)
//  - hasMeasurement: if false, the track simply coasts and ages/decays
//  - observerPosKm: observer position in world space (km)
//  - bearingDirWorld: (not necessarily normalized) direction in world space
//  - measWeight: relative confidence (>= 0). Typical range: 0..1.
inline BearingTrackResult updateBearingTrack(BearingTrack3d& track,
                                             double dtSec,
                                             bool hasMeasurement,
                                             const math::Vec3d& observerPosKm,
                                             const math::Vec3d& bearingDirWorld,
                                             double measWeight = 1.0,
                                             const BearingTrackParams& params = {}) {
  if (!std::isfinite(dtSec) || dtSec < 0.0) dtSec = 0.0;

  const math::Vec3d prevPos = track.posKm;

  // Predict/coast.
  if (track.initialized && dtSec > 0.0) {
    track.posKm += track.velKmS * dtSec;
  }

  // Age always advances.
  if (dtSec > 0.0) track.ageSinceMeasSec += dtSec;

  // Exponential decay of accumulators.
  const double decay = halfLifeDecay(params.observationHalfLifeSec, dtSec);
  if (decay != 1.0) {
    track.A *= decay;
    track.b *= decay;
    track.weight *= decay;
    track.residualWeight *= decay;
    track.residualSqSum *= decay;
  }

  bool usedMeas = false;
  math::Vec3d dUnit{0, 0, 0};
  double wUsed = 0.0;

  if (hasMeasurement && math::isFinite(observerPosKm) && math::isFinite(bearingDirWorld)) {
    dUnit = math::safeNormalized(bearingDirWorld, {0, 0, 0}, 1e-18);

    if (dUnit.lengthSq() > 1e-18) {
      wUsed = std::max(0.0, measWeight);
      if (std::isfinite(wUsed) && wUsed > 0.0) {
        usedMeas = true;

        // P = I - d d^T (projector onto the plane perpendicular to d).
        const math::Mat3d ddT = math::Mat3d::outerProduct(dUnit, dUnit);
        math::Mat3d P = math::Mat3d::identity() - ddT;

        track.A += P * wUsed;
        track.b += (P * observerPosKm) * wUsed;
        track.weight += wUsed;

        track.ageSinceMeasSec = 0.0;
      }
    }
  }

  // Determine if the system is geometrically solvable.
  bool canSolve = false;
  if (track.weight >= params.minEffectiveWeight) {
    const double scale = std::max(1.0, maxAbsElement(track.A));
    double detEps = params.determinantEps;
    if (!std::isfinite(detEps) || detEps < 0.0) detEps = 0.0;

    const double detThresh = detEps * scale * scale * scale;
    const double detAbs = std::abs(track.A.determinant());
    if (detAbs > detThresh) {
      canSolve = true;
    }
  }

  bool solved = false;
  math::Vec3d solvedPos = track.posKm;

  if (canSolve) {
    math::Mat3d Areg = track.A;

    // Add gentle diagonal regularization for numeric stability.
    double lambda = params.solveRegularization;
    if (!std::isfinite(lambda) || lambda < 0.0) lambda = 0.0;

    const double tr = std::abs(Areg.trace());
    lambda *= std::max(1.0, tr / 3.0);
    if (lambda > 0.0) Areg.addToDiagonal(lambda);

    // Use an absolute epsilon derived from the matrix scale.
    const double scale = std::max(1.0, maxAbsElement(Areg));
    double detEps = params.determinantEps;
    if (!std::isfinite(detEps) || detEps < 0.0) detEps = 0.0;
    const double detEpsAbs = std::max(1.0e-18, detEps * scale * scale * scale);

    math::Vec3d p{0, 0, 0};
    if (Areg.solve(p, track.b, detEpsAbs) && math::isFinite(p)) {
      solved = true;
      solvedPos = p;
    }
  }

  // Update state when solvable.
  if (solved) {
    // Velocity estimate from successive position solutions.
    if (track.initialized && dtSec > 1.0e-6) {
      const math::Vec3d measVel = (solvedPos - prevPos) / dtSec;

      const double velAlpha = 1.0 - halfLifeDecay(params.velHalfLifeSec, dtSec);
      const double a = clamp01(velAlpha);
      track.velKmS = lerp(track.velKmS, measVel, a);
    } else if (!track.initialized) {
      track.velKmS = {0, 0, 0};
    }

    if (params.maxSpeedKmS > 0.0 && std::isfinite(params.maxSpeedKmS)) {
      track.velKmS = math::clampMagnitude(track.velKmS, params.maxSpeedKmS);
    }

    track.posKm = solvedPos;
    track.initialized = true;
  }

  // Update residual fit quality whenever we have an estimate and a measurement.
  if (track.initialized && usedMeas) {
    const math::Vec3d rel = (track.posKm - observerPosKm);
    // Distance from point to line: |(p - o) x d| (d must be unit).
    const double dist = math::cross(rel, dUnit).length();
    if (std::isfinite(dist)) {
      track.residualSqSum += wUsed * dist * dist;
      track.residualWeight += wUsed;
    }
  }

  if (track.residualWeight > 1e-12) {
    const double sigma2 = std::max(0.0, track.residualSqSum / track.residualWeight);
    track.sigmaKm = std::sqrt(sigma2);
  }

  if (!std::isfinite(track.sigmaKm)) track.sigmaKm = params.sigmaMaxKm;
  track.sigmaKm = std::clamp(track.sigmaKm, params.sigmaMinKm, params.sigmaMaxKm);

  BearingTrackResult out{};
  out.valid = track.initialized;
  out.hasMeasurement = usedMeas;
  out.posKm = track.posKm;
  out.velKmS = track.velKmS;
  out.sigmaKm = track.sigmaKm;
  out.ageSinceMeasSec = track.ageSinceMeasSec;
  return out;
}

} // namespace stellar::sim
