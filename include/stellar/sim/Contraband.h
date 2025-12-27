#pragma once

#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"

#include <array>
#include <cstddef>
#include <string>

namespace stellar::sim {

// Bitmask helpers for simple contraband/legality rules.
// Commodity count is intentionally small in the prototype, so a 32-bit mask is sufficient.
inline core::u32 commodityBit(econ::CommodityId cid) {
  return (core::u32)1u << (core::u32)cid;
}

// Deterministic per-faction illegality mask.
//
// This is used by both the SDL prototype and the headless modules (missions / tooling) so that
// "illegal here" is consistent across UI, tests, and any CLI utilities.
//
// Design: each non-zero faction bans 1 "vice" good and, with some probability,
// bans 1 "controlled tech" good.
inline core::u32 illegalCommodityMask(core::u64 universeSeed, core::u32 factionId) {
  if (factionId == 0) return 0u;

  // Large odd constant to improve bit diffusion. Keep stable for save/behavior continuity.
  core::SplitMix64 r(universeSeed ^ (core::u64)factionId * 0x9E3779B97F4A7C15ull);

  core::u32 mask = 0u;

  // Keep these lists small and "gamey" for now: it creates consistent but varied
  // contraband rules without needing a full law simulation.
  const std::array<econ::CommodityId, 3> vice{
      econ::CommodityId::Luxury,
      econ::CommodityId::Medicine,
      econ::CommodityId::Stimulants,
  };
  const std::array<econ::CommodityId, 3> tech{
      econ::CommodityId::Electronics,
      econ::CommodityId::Machinery,
      econ::CommodityId::Weapons,
  };

  mask |= commodityBit(vice[(std::size_t)(r.nextU32() % (core::u32)vice.size())]);

  if (r.nextUnit() < 0.70) {
    mask |= commodityBit(tech[(std::size_t)(r.nextU32() % (core::u32)tech.size())]);
  }

  return mask;
}

inline bool isIllegalCommodity(core::u64 universeSeed, core::u32 factionId, econ::CommodityId cid) {
  if (factionId == 0) return false;
  return (illegalCommodityMask(universeSeed, factionId) & commodityBit(cid)) != 0u;
}

inline std::string illegalCommodityListString(core::u64 universeSeed, core::u32 factionId) {
  const core::u32 mask = illegalCommodityMask(universeSeed, factionId);
  if (mask == 0u) return "None";

  std::string out;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    if ((mask & ((core::u32)1u << (core::u32)i)) == 0u) continue;
    if (!out.empty()) out += ", ";
    out += std::string(econ::commodityName((econ::CommodityId)i));
  }
  return out.empty() ? "None" : out;
}

inline bool hasIllegalCargo(core::u64 universeSeed,
                            core::u32 factionId,
                            const std::array<double, econ::kCommodityCount>& cargo) {
  const core::u32 mask = illegalCommodityMask(universeSeed, factionId);
  if (mask == 0u) return false;

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    if ((mask & ((core::u32)1u << (core::u32)i)) != 0u) {
      if (cargo[i] > 1e-6) return true;
    }
  }
  return false;
}

} // namespace stellar::sim
