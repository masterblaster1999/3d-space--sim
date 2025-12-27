#include "stellar/econ/RoutePlanner.h"

#include "stellar/econ/Commodity.h"

#include <algorithm>
#include <cmath>

namespace stellar::econ {

static constexpr std::size_t idx(CommodityId id) { return static_cast<std::size_t>(id); }

std::vector<RouteOpportunity> bestRoutes(const StationEconomyState& fromState,
                                        const StationEconomyModel& fromModel,
                                        const StationEconomyState& toState,
                                        const StationEconomyModel& toModel,
                                        double bidAskSpread,
                                        std::size_t maxResults) {
  std::vector<RouteOpportunity> out;
  out.reserve(kCommodityCount);

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const CommodityId cid = static_cast<CommodityId>(i);
    const auto qFrom = quote(fromState, fromModel, cid, bidAskSpread);
    const auto qTo   = quote(toState, toModel, cid, bidAskSpread);

    const double profit = qTo.bid - qFrom.ask;

    // Feasibility: don't recommend routes where you cannot actually buy/sell.
    const double unitsFrom = std::max(0.0, qFrom.inventory);
    const double unitsToSpace = std::max(0.0, toModel.capacity[i] - toState.inventory[i]);
    const double unitsPossible = std::min(unitsFrom, unitsToSpace);

    if (profit > 0.0 && unitsPossible > 1e-9) {
      const double massKg = std::max(1e-6, commodityDef(cid).massKg);

      RouteOpportunity r{};
      r.commodity = cid;
      r.profitPerUnit = profit;
      r.buyPrice = qFrom.ask;
      r.sellPrice = qTo.bid;

      r.unitsFrom = unitsFrom;
      r.unitsToSpace = unitsToSpace;
      r.unitsPossible = unitsPossible;

      r.unitMassKg = massKg;
      r.profitPerKg = profit / massKg;
      r.profitTotal = profit * unitsPossible;

      // Defaults: no fees applied.
      r.netProfitPerUnit = r.profitPerUnit;
      r.netProfitTotal = r.profitTotal;

      out.push_back(r);
    }
  }

  std::sort(out.begin(), out.end(), [](const RouteOpportunity& a, const RouteOpportunity& b) {
    return a.profitPerUnit > b.profitPerUnit;
  });

  if (out.size() > maxResults) out.resize(maxResults);
  return out;
}

std::vector<RouteOpportunity> bestRoutesForCargo(const StationEconomyState& fromState,
                                                const StationEconomyModel& fromModel,
                                                const StationEconomyState& toState,
                                                const StationEconomyModel& toModel,
                                                double cargoCapacityKg,
                                                double fromFeeRate,
                                                double toFeeRate,
                                                double bidAskSpread,
                                                std::size_t maxResults) {
  std::vector<RouteOpportunity> out;
  out.reserve(kCommodityCount);

  cargoCapacityKg = std::max(0.0, cargoCapacityKg);
  fromFeeRate = std::clamp(fromFeeRate, 0.0, 1.0);
  toFeeRate = std::clamp(toFeeRate, 0.0, 1.0);

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const CommodityId cid = static_cast<CommodityId>(i);
    const auto qFrom = quote(fromState, fromModel, cid, bidAskSpread);
    const auto qTo   = quote(toState, toModel, cid, bidAskSpread);

    const double rawProfit = qTo.bid - qFrom.ask;
    if (rawProfit <= 0.0) continue;

    // Station feasibility limits.
    const double unitsFrom = std::max(0.0, qFrom.inventory);
    const double unitsToSpace = std::max(0.0, toModel.capacity[i] - toState.inventory[i]);
    double unitsPossible = std::min(unitsFrom, unitsToSpace);
    if (unitsPossible <= 1e-9) continue;

    const double massKg = std::max(1e-6, commodityDef(cid).massKg);

    // Cargo feasibility limit.
    if (cargoCapacityKg > 0.0) {
      const double cargoUnits = std::floor(cargoCapacityKg / massKg + 1e-9);
      unitsPossible = std::min(unitsPossible, std::max(0.0, cargoUnits));
      if (unitsPossible <= 1e-9) continue;
    }

    // Net-of-fees profit.
    const double netBuy = qFrom.ask * (1.0 + fromFeeRate);
    const double netSell = qTo.bid * (1.0 - toFeeRate);
    const double netProfit = netSell - netBuy;
    if (netProfit <= 0.0) continue;

    RouteOpportunity r{};
    r.commodity = cid;
    r.profitPerUnit = rawProfit;
    r.buyPrice = qFrom.ask;
    r.sellPrice = qTo.bid;
    r.unitsFrom = unitsFrom;
    r.unitsToSpace = unitsToSpace;
    r.unitsPossible = unitsPossible;
    r.unitMassKg = massKg;
    r.profitPerKg = rawProfit / massKg;
    r.profitTotal = rawProfit * unitsPossible;

    r.feeFrom = fromFeeRate;
    r.feeTo = toFeeRate;
    r.netProfitPerUnit = netProfit;
    r.netProfitTotal = netProfit * unitsPossible;

    out.push_back(r);
  }

  std::sort(out.begin(), out.end(), [](const RouteOpportunity& a, const RouteOpportunity& b) {
    // Primary: trip profit.
    if (a.netProfitTotal != b.netProfitTotal) return a.netProfitTotal > b.netProfitTotal;
    // Secondary: per-unit profit.
    return a.netProfitPerUnit > b.netProfitPerUnit;
  });

  if (out.size() > maxResults) out.resize(maxResults);
  return out;
}

} // namespace stellar::econ
