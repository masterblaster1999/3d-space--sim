#include "stellar/proc/Noise.h"

#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"

#include <cmath>

namespace stellar::proc {
namespace {
  stellar::core::u64 mix64(stellar::core::u64 z) {
    // SplitMix64 finalizer.
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }

  double u64ToUnitDouble(stellar::core::u64 x) {
    // 53 bits to [0,1)
    const stellar::core::u64 mantissa = x >> 11;
    return static_cast<double>(mantissa) * (1.0 / 9007199254740992.0);
  }

  double lerp(double a, double b, double t) {
    return a + (b - a) * t;
  }
} // namespace

double valueNoise2D(int x, int y, stellar::core::u64 seed) {
  using stellar::core::u64;

  u64 h = seed;
  h = stellar::core::hashCombine(h, static_cast<u64>(static_cast<std::uint32_t>(x)));
  h = stellar::core::hashCombine(h, static_cast<u64>(static_cast<std::uint32_t>(y)));
  h = mix64(h);
  return u64ToUnitDouble(h);
}

double valueNoise2D(double x, double y, stellar::core::u64 seed) {
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const int x1 = x0 + 1;
  const int y1 = y0 + 1;

  const double fx = x - static_cast<double>(x0);
  const double fy = y - static_cast<double>(y0);

  const double sx = stellar::math::smoothstep(fx);
  const double sy = stellar::math::smoothstep(fy);

  const double n00 = valueNoise2D(x0, y0, seed);
  const double n10 = valueNoise2D(x1, y0, seed);
  const double n01 = valueNoise2D(x0, y1, seed);
  const double n11 = valueNoise2D(x1, y1, seed);

  const double ix0 = lerp(n00, n10, sx);
  const double ix1 = lerp(n01, n11, sx);
  return lerp(ix0, ix1, sy);
}

double fbm2D(double x, double y, stellar::core::u64 seed, int octaves, double lacunarity, double gain) {
  double sum = 0.0;
  double amp = 0.5;
  double freq = 1.0;

  double norm = 0.0;

  for (int i = 0; i < octaves; ++i) {
    sum += amp * valueNoise2D(x * freq, y * freq, seed + static_cast<stellar::core::u64>(i) * 1013ull);
    norm += amp;
    amp *= gain;
    freq *= lacunarity;
  }

  if (norm > 0.0) sum /= norm;
  return sum;
}

} // namespace stellar::proc
