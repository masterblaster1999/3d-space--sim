#pragma once

#include "stellar/sim/SaveGame.h"

#include <algorithm>

namespace stellar::sim {

inline double clampReputation(double rep) {
  return std::clamp(rep, -100.0, 100.0);
}

inline double reputationNorm(double rep) {
  return clampReputation(rep) / 100.0;
}

// Positive rep reduces fees, negative rep increases fees.
// Returns a fee rate clamped to a reasonable range.
inline double applyReputationToFeeRate(double baseFeeRate, double rep) {
  const double kMaxEffect = 0.35; // +/- 35%
  const double eff = baseFeeRate * (1.0 - kMaxEffect * reputationNorm(rep));
  return std::clamp(eff, 0.0, 0.25);
}

inline double getReputation(const SaveGame& s, core::u32 factionId) {
  for (const auto& r : s.reputation) {
    if (r.factionId == factionId) return r.rep;
  }
  return 0.0;
}

inline void setReputation(SaveGame& s, core::u32 factionId, double rep) {
  rep = clampReputation(rep);
  for (auto& r : s.reputation) {
    if (r.factionId == factionId) {
      r.rep = rep;
      return;
    }
  }
  s.reputation.push_back(FactionReputation{factionId, rep});
}

inline void addReputation(SaveGame& s, core::u32 factionId, double delta) {
  setReputation(s, factionId, getReputation(s, factionId) + delta);
}

inline double getBounty(const SaveGame& s, core::u32 factionId) {
  for (const auto& b : s.bounties) {
    if (b.factionId == factionId) return b.bountyCr;
  }
  return 0.0;
}

inline void setBounty(SaveGame& s, core::u32 factionId, double bountyCr) {
  bountyCr = std::max(0.0, bountyCr);
  for (auto& b : s.bounties) {
    if (b.factionId == factionId) {
      b.bountyCr = bountyCr;
      return;
    }
  }
  if (bountyCr > 0.0) s.bounties.push_back(FactionBounty{factionId, bountyCr});
}

inline void addBounty(SaveGame& s, core::u32 factionId, double deltaCr) {
  setBounty(s, factionId, getBounty(s, factionId) + deltaCr);
}

} // namespace stellar::sim
