#pragma once

#include "stellar/econ/Commodity.h"

#include <array>

namespace stellar::econ {

// Convenience helper used by gameplay, tools, and mission logic.
// Cargo is represented as an array of commodity *units*.
inline double cargoMassKg(const std::array<double, econ::kCommodityCount>& cargo) {
  double kg = 0.0;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const auto cid = static_cast<econ::CommodityId>(i);
    const double units = cargo[i];
    if (units <= 0.0) continue;
    kg += units * econ::commodityDef(cid).massKg;
  }
  return kg;
}

} // namespace stellar::econ
