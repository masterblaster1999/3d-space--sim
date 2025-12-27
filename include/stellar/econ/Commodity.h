#pragma once

#include "stellar/core/Types.h"

#include <array>
#include <string_view>

namespace stellar::econ {

enum class CommodityId : core::u16 {
  Food = 0,
  Water,
  Ore,
  Metals,
  Fuel,
  Machinery,
  Medicine,
  Electronics,
  Luxury,

  // Higher-value / specialty goods.
  Weapons,
  Stimulants,

  Count
};

constexpr std::size_t kCommodityCount = static_cast<std::size_t>(CommodityId::Count);

struct CommodityDef {
  CommodityId id{};
  const char* code{};     // short symbol for save files / UI
  const char* name{};     // display name
  double basePrice{};     // "credits" per unit (mid price baseline)
  double massKg{};        // kg per unit (for cargo capacity later)
};

const std::array<CommodityDef, kCommodityCount>& commodityTable();
const CommodityDef& commodityDef(CommodityId id);
std::string_view commodityName(CommodityId id);
std::string_view commodityCode(CommodityId id);

} // namespace stellar::econ
