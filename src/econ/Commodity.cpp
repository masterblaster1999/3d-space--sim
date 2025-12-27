#include "stellar/econ/Commodity.h"

namespace stellar::econ {

static const std::array<CommodityDef, kCommodityCount> kTable = {{
  {CommodityId::Food,        "FOOD", "Food",        12.0,  1.0},
  {CommodityId::Water,       "H2O",  "Water",        6.0,  1.0},
  {CommodityId::Ore,         "ORE",  "Ore",          9.0,  4.0},
  {CommodityId::Metals,      "MET",  "Metals",      22.0,  3.0},
  {CommodityId::Fuel,        "FUEL", "Fuel",        35.0,  1.0},
  {CommodityId::Machinery,   "MACH", "Machinery",   55.0,  2.0},
  {CommodityId::Medicine,    "MED",  "Medicine",    75.0,  0.5},
  {CommodityId::Electronics, "ELEC", "Electronics", 85.0,  0.5},
  {CommodityId::Luxury,      "LUX",  "Luxury",     120.0,  0.2},

  {CommodityId::Weapons,     "ARMS", "Weapons",    150.0,  1.6},
  {CommodityId::Stimulants,  "STIM", "Stimulants", 165.0,  0.2},
}};

const std::array<CommodityDef, kCommodityCount>& commodityTable() { return kTable; }

const CommodityDef& commodityDef(CommodityId id) {
  return kTable[static_cast<std::size_t>(id)];
}

std::string_view commodityName(CommodityId id) {
  return commodityDef(id).name;
}

std::string_view commodityCode(CommodityId id) {
  return commodityDef(id).code;
}

} // namespace stellar::econ
