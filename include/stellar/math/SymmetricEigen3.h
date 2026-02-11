#pragma once

#include "stellar/math/Mat3.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace stellar::math {

// -----------------------------------------------------------------------------
// Symmetric eigen decomposition for 3x3 matrices (Jacobi iterations)
// -----------------------------------------------------------------------------
//
// This is a small, dependency-free helper for tools/sim code that need
// principal-axis information (e.g. uncertainty ellipsoids).
//
// Notes:
//  - Input is treated as symmetric; we symmetrize (A + A^T)/2 for robustness.
//  - Eigenvalues are returned in ascending order.
//  - Eigenvectors are returned as columns of `eigenvectors`.
//  - Eigenvector signs are fixed deterministically (largest-magnitude component
//    is made positive) so results are stable across platforms.

struct SymmetricEigen3Result {
  bool valid{false};
  Vec3d eigenvalues{0, 0, 0};
  Mat3d eigenvectors{Mat3d::identity()};
};

namespace detail {

inline bool finite3x3(const double a[3][3]) {
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      if (!std::isfinite(a[r][c])) return false;
    }
  }
  return true;
}

inline void normalizeSafe(Vec3d& v) {
  const double lsq = v.lengthSq();
  if (!(lsq > 1e-30) || !std::isfinite(lsq)) {
    v = {0, 0, 0};
    return;
  }
  v *= 1.0 / std::sqrt(lsq);
}

inline void fixDeterministicSign(Vec3d& v) {
  const double ax = std::abs(v.x);
  const double ay = std::abs(v.y);
  const double az = std::abs(v.z);

  if (ax >= ay && ax >= az) {
    if (v.x < 0.0) v = -v;
  } else if (ay >= ax && ay >= az) {
    if (v.y < 0.0) v = -v;
  } else {
    if (v.z < 0.0) v = -v;
  }
}

inline void swapCols(double v[3][3], int a, int b) {
  for (int r = 0; r < 3; ++r) {
    std::swap(v[r][a], v[r][b]);
  }
}

} // namespace detail

inline SymmetricEigen3Result symmetricEigenDecomposition3x3(const Mat3d& A) {
  SymmetricEigen3Result out{};
  if (!A.isFinite()) return out;

  // Convert to row-major and symmetrize.
  const double a00 = A.m[0], a01 = A.m[3], a02 = A.m[6];
  const double a10 = A.m[1], a11 = A.m[4], a12 = A.m[7];
  const double a20 = A.m[2], a21 = A.m[5], a22 = A.m[8];

  double a[3][3] = {
    {a00, 0.5 * (a01 + a10), 0.5 * (a02 + a20)},
    {0.5 * (a10 + a01), a11, 0.5 * (a12 + a21)},
    {0.5 * (a20 + a02), 0.5 * (a21 + a12), a22},
  };

  if (!detail::finite3x3(a)) return out;

  // Eigenvectors accumulator (row-major). Start with identity.
  double v[3][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
  };

  auto offDiagMaxAbs = [&]() {
    return std::max({
      std::abs(a[0][1]),
      std::abs(a[0][2]),
      std::abs(a[1][2]),
    });
  };

  const double scale = std::max({
    1.0,
    std::abs(a[0][0]), std::abs(a[1][1]), std::abs(a[2][2]),
    std::abs(a[0][1]), std::abs(a[0][2]), std::abs(a[1][2]),
  });

  const double eps = 1e-14 * scale;

  // Jacobi iterations for 3x3 symmetric matrices.
  for (int iter = 0; iter < 32; ++iter) {
    // Find the largest off-diagonal element among (0,1), (0,2), (1,2).
    int p = 0;
    int q = 1;
    double maxVal = std::abs(a[0][1]);

    const double a02abs = std::abs(a[0][2]);
    if (a02abs > maxVal) {
      maxVal = a02abs;
      p = 0;
      q = 2;
    }

    const double a12abs = std::abs(a[1][2]);
    if (a12abs > maxVal) {
      maxVal = a12abs;
      p = 1;
      q = 2;
    }

    if (maxVal <= eps) break;

    const double app = a[p][p];
    const double aqq = a[q][q];
    const double apq = a[p][q];

    // Compute Jacobi rotation.
    const double phi = 0.5 * std::atan2(2.0 * apq, (aqq - app));
    const double c = std::cos(phi);
    const double s = std::sin(phi);

    // Apply rotation to A: A = J^T A J.
    for (int r = 0; r < 3; ++r) {
      const double arp = a[r][p];
      const double arq = a[r][q];
      a[r][p] = c * arp - s * arq;
      a[r][q] = s * arp + c * arq;
    }

    for (int ccol = 0; ccol < 3; ++ccol) {
      const double apc = a[p][ccol];
      const double aqc = a[q][ccol];
      a[p][ccol] = c * apc - s * aqc;
      a[q][ccol] = s * apc + c * aqc;
    }

    // Enforce exact symmetry for numerical stability.
    a[0][1] = a[1][0] = 0.5 * (a[0][1] + a[1][0]);
    a[0][2] = a[2][0] = 0.5 * (a[0][2] + a[2][0]);
    a[1][2] = a[2][1] = 0.5 * (a[1][2] + a[2][1]);

    // Apply the same rotation to eigenvectors accumulator V: V = V J.
    for (int r = 0; r < 3; ++r) {
      const double vrp = v[r][p];
      const double vrq = v[r][q];
      v[r][p] = c * vrp - s * vrq;
      v[r][q] = s * vrp + c * vrq;
    }
  }

  // Extract eigenvalues from the diagonal.
  double evals[3] = {a[0][0], a[1][1], a[2][2]};

  // Sort eigenvalues ascending and permute eigenvectors accordingly.
  int order[3] = {0, 1, 2};
  std::sort(order, order + 3, [&](int i, int j) {
    const double di = evals[i];
    const double dj = evals[j];
    if (std::abs(di - dj) > 1e-12) return di < dj;
    return i < j;
  });

  // Reorder evals and columns of v.
  double evalsSorted[3] = {evals[order[0]], evals[order[1]], evals[order[2]]};

  double vSorted[3][3] = {
    {v[0][order[0]], v[0][order[1]], v[0][order[2]]},
    {v[1][order[0]], v[1][order[1]], v[1][order[2]]},
    {v[2][order[0]], v[2][order[1]], v[2][order[2]]},
  };

  Vec3d c0{vSorted[0][0], vSorted[1][0], vSorted[2][0]};
  Vec3d c1{vSorted[0][1], vSorted[1][1], vSorted[2][1]};
  Vec3d c2{vSorted[0][2], vSorted[1][2], vSorted[2][2]};

  detail::normalizeSafe(c0);
  detail::normalizeSafe(c1);
  detail::normalizeSafe(c2);

  // Deterministic sign convention.
  detail::fixDeterministicSign(c0);
  detail::fixDeterministicSign(c1);
  detail::fixDeterministicSign(c2);

  out.eigenvalues = {evalsSorted[0], evalsSorted[1], evalsSorted[2]};
  out.eigenvectors = Mat3d::fromColumns(c0, c1, c2);
  out.valid = out.eigenvectors.isFinite() && isFinite(out.eigenvalues);
  return out;
}

} // namespace stellar::math
