#pragma once

#include "stellar/core/Hash.h"
#include "stellar/core/Types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <string_view>
#include <type_traits>

namespace stellar::core {

// SplitMix64: fast, simple PRNG with 64-bit state.
// Great for procedural generation & seeding other RNGs.
// Not suitable for crypto.
class SplitMix64 {
public:
  explicit SplitMix64(u64 seed) : m_state(seed) {}

  u64 nextU64() {
    u64 z = (m_state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }

  u32 nextU32() { return static_cast<u32>(nextU64() >> 32); }

  // [0, 1)
  double nextDouble01() {
    // 53 random bits -> double in [0,1)
    const u64 r = nextU64();
    const u64 mantissa = r >> 11; // keep top 53 bits
    return static_cast<double>(mantissa) * (1.0 / 9007199254740992.0); // 2^53
  }

  // Uniform real in [min, max)
  double uniform(double min, double max) {
    return min + (max - min) * nextDouble01();
  }

  // Uniform integer in [minInclusive, maxInclusive]
  template <typename Int>
  Int range(Int minInclusive, Int maxInclusive) {
    static_assert(std::is_integral_v<Int>);
    if (maxInclusive < minInclusive) {
      std::swap(minInclusive, maxInclusive);
    }
    const u64 span = static_cast<u64>(maxInclusive) - static_cast<u64>(minInclusive) + 1ull;

    // Rejection sampling to avoid modulo bias.
    const u64 limit = (std::numeric_limits<u64>::max() / span) * span;
    u64 r = 0;
    do {
      r = nextU64();
    } while (r >= limit);

    return static_cast<Int>(static_cast<u64>(minInclusive) + (r % span));
  }

  bool chance(double probability01) {
    probability01 = std::clamp(probability01, 0.0, 1.0);
    return nextDouble01() < probability01;
  }

  template <typename T, std::size_t Extent>
  const T& pick(std::span<const T, Extent> items) {
    const auto idx = range<std::size_t>(0u, items.size() - 1u);
    return items[idx];
  }

  u64 state() const { return m_state; }

private:
  u64 m_state = 0;
};

// Convenience: stable seed from a string
inline u64 seedFromString(std::string_view s) { return fnv1a64(s); }

// Convenience: derive a child seed from a parent seed and a tag.
// Great for deterministic sub-generators (system -> planet -> terrain, etc.)
inline u64 deriveSeed(u64 parent, u64 tag) { return hashCombine(parent, tag); }
inline u64 deriveSeed(u64 parent, std::string_view tag) { return hashCombine(parent, fnv1a64(tag)); }

} // namespace stellar::core
