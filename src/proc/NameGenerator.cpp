#include "stellar/proc/NameGenerator.h"

#include <algorithm>
#include <cctype>
#include <span>
#include <string>
#include <vector>

namespace stellar::proc {
namespace {
  constexpr const char* consonants[] = {
    "b","c","d","f","g","h","j","k","l","m","n","p","qu","r","s","t","v","w","x","y","z",
    "br","cr","dr","fr","gr","kr","pr","tr","vr",
    "ch","sh","th","st","sk","sl","sm","sn","sp","sw"
  };

  constexpr const char* vowels[] = {
    "a","e","i","o","u",
    "ae","ai","ao","au",
    "ea","ei","eo","eu",
    "ia","ie","io","iu",
    "oa","oe","oi","ou",
    "ua","ue","ui","uo"
  };

  constexpr const char* endings[] = {
    "","n","s","r","th","m","x","on","ar","is","us","or","ion","ara","ea","ium"
  };

  std::string capitalize(std::string s) {
    if (!s.empty()) {
      s[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
    }
    return s;
  }
} // namespace

NameGenerator::NameGenerator(stellar::core::u64 seed) : m_seed(seed) {}

std::string NameGenerator::makeName(stellar::core::SplitMix64& rng, int minSyllables, int maxSyllables) const {
  minSyllables = std::max(1, minSyllables);
  maxSyllables = std::max(minSyllables, maxSyllables);

  const int syllables = rng.range<int>(minSyllables, maxSyllables);

  std::string out;
  out.reserve(static_cast<std::size_t>(syllables) * 4u + 4u);

  for (int i = 0; i < syllables; ++i) {
    // Occasionally start with a vowel for variety.
    const bool startWithVowel = (i == 0) && rng.chance(0.2);

    if (!startWithVowel) {
      const auto* c = rng.pick(std::span(consonants));
      out += c;
    }

    const auto* v = rng.pick(std::span(vowels));
    out += v;

    // Small chance of consonant cluster mid-name
    if (rng.chance(0.15)) {
      const auto* c2 = rng.pick(std::span(consonants));
      out += c2;
    }
  }

  // Optional ending
  if (rng.chance(0.6)) {
    out += rng.pick(std::span(endings));
  }

  // Optional numeric suffix (rare)
  if (rng.chance(0.05)) {
    const int n = rng.range<int>(1, 999);
    out += "-";
    out += std::to_string(n);
  }

  // Use stored seed to lightly perturb (so two NameGenerator instances can differ).
  // (This is intentionally subtle.)
  if ((m_seed & 1ull) == 0ull && !out.empty()) {
    // Very rarely swap two characters
    if (out.size() > 4 && rng.chance(0.02)) {
      const auto a = rng.range<std::size_t>(1u, out.size() - 2u);
      const auto b = rng.range<std::size_t>(1u, out.size() - 2u);
      std::swap(out[a], out[b]);
    }
  }

  return capitalize(out);
}

} // namespace stellar::proc
