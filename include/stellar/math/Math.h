#pragma once

#include <algorithm>
#include <cmath>

namespace stellar::math {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double twoPi = 2.0 * pi;

template <typename T>
constexpr T clamp(T v, T lo, T hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}

template <typename T>
constexpr T lerp(const T& a, const T& b, double t) {
  return static_cast<T>(a + (b - a) * t);
}

inline double degToRad(double deg) { return deg * (pi / 180.0); }
inline double radToDeg(double rad) { return rad * (180.0 / pi); }

inline double smoothstep(double t) {
  t = std::clamp(t, 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}

} // namespace stellar::math
