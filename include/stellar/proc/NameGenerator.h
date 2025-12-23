#pragma once

#include "stellar/core/Random.h"

#include <string>

namespace stellar::proc {

// A tiny deterministic syllable-based name generator.
// You can replace this later with Markov chains, real star catalogs, etc.
class NameGenerator {
public:
  explicit NameGenerator(stellar::core::u64 seed);

  std::string makeName(stellar::core::SplitMix64& rng, int minSyllables = 2, int maxSyllables = 4) const;

private:
  stellar::core::u64 m_seed = 0;
};

} // namespace stellar::proc
