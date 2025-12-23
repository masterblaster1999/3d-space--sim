#include "stellar/core/Hash.h"

namespace stellar::core {

u64 fnv1a64(std::string_view text) {
  // FNV-1a 64-bit parameters
  constexpr u64 offsetBasis = 14695981039346656037ull;
  constexpr u64 prime       = 1099511628211ull;

  u64 hash = offsetBasis;
  for (char c : text) {
    hash ^= static_cast<u64>(static_cast<unsigned char>(c));
    hash *= prime;
  }
  return hash;
}

u64 hashCombine(u64 a, u64 b) {
  // A common 64-bit hash combine pattern.
  // Similar to boost::hash_combine but 64-bit.
  a ^= b + 0x9e3779b97f4a7c15ull + (a << 6) + (a >> 2);
  return a;
}

} // namespace stellar::core
