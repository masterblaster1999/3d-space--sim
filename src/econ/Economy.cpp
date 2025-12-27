#include "stellar/econ/Economy.h"

#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::econ {

static constexpr std::size_t idx(CommodityId id) { return static_cast<std::size_t>(id); }

void StationEconomyState::clampToCapacity(const StationEconomyModel& model) {
  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    inventory[i] = std::max(0.0, inventory[i]);
    inventory[i] = std::min(inventory[i], model.capacity[i]);
  }
}

static void setAll(std::array<double, kCommodityCount>& a, double v) {
  for (double& x : a) x = v;
}

StationEconomyModel makeEconomyModel(StationType type, double bias) {
  StationEconomyModel m{};
  m.type = type;

  // Keep station inventories intentionally a bit tight so that trading/missions can
  // actually create shortages and price spikes. (This is a gameplay-oriented model,
  // not a realistic macro-economy.)
  const double kStockScale = 0.65;

  setAll(m.productionPerDay, 0.0);
  setAll(m.consumptionPerDay, 0.0);
  setAll(m.desiredStock, 0.0);
  setAll(m.capacity, 0.0);

  auto D = [&](CommodityId id, double desired) {
    desired *= kStockScale;
    m.desiredStock[idx(id)] = desired;
    m.capacity[idx(id)] = desired * 1.8;
  };
  auto P = [&](CommodityId id, double prod) { m.productionPerDay[idx(id)] = prod; };
  auto C = [&](CommodityId id, double cons) { m.consumptionPerDay[idx(id)] = cons; };

  // Common baselines
  D(CommodityId::Food, 200);
  D(CommodityId::Water, 200);
  D(CommodityId::Fuel, 120);
  D(CommodityId::Medicine, 40);
  D(CommodityId::Electronics, 60);
  D(CommodityId::Machinery, 60);
  D(CommodityId::Metals, 140);
  D(CommodityId::Ore, 180);
  D(CommodityId::Luxury, 25);
  D(CommodityId::Weapons, 35);
  D(CommodityId::Stimulants, 18);

  switch (type) {
    case StationType::Outpost:
      C(CommodityId::Food, 10);
      C(CommodityId::Water, 10);
      C(CommodityId::Fuel, 4);
      C(CommodityId::Medicine, 1);
      break;

    case StationType::Agricultural:
      // produces food; consumes machinery/fuel.
      P(CommodityId::Food, 90 + 30 * std::max(0.0, -bias));
      P(CommodityId::Water, 12);
      // Small pharma/agro-chemical output.
      P(CommodityId::Stimulants, 6.0 + 3.0 * std::max(0.0, -bias));
      C(CommodityId::Fuel, 6);
      C(CommodityId::Machinery, 2);
      C(CommodityId::Electronics, 1);
      C(CommodityId::Medicine, 0.5);
      D(CommodityId::Food, 700);
      D(CommodityId::Water, 350);
      D(CommodityId::Stimulants, 35);
      break;

    case StationType::Mining:
      P(CommodityId::Ore, 80);
      C(CommodityId::Fuel, 10);
      C(CommodityId::Machinery, 3);
      C(CommodityId::Food, 6);
      C(CommodityId::Water, 6);
      D(CommodityId::Ore, 600);
      break;

    case StationType::Refinery:
      C(CommodityId::Ore, 70);
      C(CommodityId::Fuel, 4);
      P(CommodityId::Metals, 40);
      P(CommodityId::Fuel, 12);
      D(CommodityId::Ore, 500);
      D(CommodityId::Metals, 500);
      D(CommodityId::Fuel, 300);
      break;

    case StationType::Industrial:
      C(CommodityId::Metals, 35);
      C(CommodityId::Fuel, 8);
      P(CommodityId::Machinery, 18 + 8 * std::max(0.0, bias));
      P(CommodityId::Electronics, 10 + 6 * std::max(0.0, bias));
      P(CommodityId::Weapons, 5.5 + 2.5 * std::max(0.0, bias));
      C(CommodityId::Stimulants, 0.6);
      D(CommodityId::Machinery, 350);
      D(CommodityId::Electronics, 300);
      D(CommodityId::Metals, 450);
      D(CommodityId::Weapons, 140);
      break;

    case StationType::Research:
      C(CommodityId::Electronics, 5);
      C(CommodityId::Medicine, 2);
      C(CommodityId::Luxury, 1);
      C(CommodityId::Stimulants, 0.8);
      C(CommodityId::Fuel, 3);
      P(CommodityId::Medicine, 4);
      D(CommodityId::Medicine, 250);
      D(CommodityId::Electronics, 250);
      D(CommodityId::Luxury, 120);
      D(CommodityId::Stimulants, 60);
      break;

    case StationType::TradeHub:
      // small "balancing" production/consumption; mostly demand-driven.
      C(CommodityId::Food, 6);
      C(CommodityId::Water, 6);
      C(CommodityId::Fuel, 4);
      C(CommodityId::Medicine, 1);
      C(CommodityId::Electronics, 2);
      C(CommodityId::Machinery, 2);
      C(CommodityId::Metals, 4);
      C(CommodityId::Ore, 4);
      C(CommodityId::Weapons, 1.2);
      C(CommodityId::Stimulants, 1.0);
      D(CommodityId::Luxury, 80);
      D(CommodityId::Weapons, 120);
      D(CommodityId::Stimulants, 90);
      break;

    case StationType::Shipyard:
      C(CommodityId::Metals, 60);
      C(CommodityId::Machinery, 30);
      C(CommodityId::Electronics, 18);
      C(CommodityId::Fuel, 10);
      C(CommodityId::Weapons, 3.0);
      D(CommodityId::Metals, 900);
      D(CommodityId::Machinery, 700);
      D(CommodityId::Electronics, 600);
      D(CommodityId::Weapons, 200);
      break;

    case StationType::Count:
      break;
  }

  // Slightly damp volatility for big hubs (more stable markets)
  if (type == StationType::TradeHub) m.priceVolatility = 0.65;
  if (type == StationType::Outpost)  m.priceVolatility = 1.10;

  return m;
}

StationEconomyState makeInitialState(const StationEconomyModel& model, core::SplitMix64& rng) {
  StationEconomyState s{};
  s.lastUpdateDay = 0.0;
  s.lastSampleDay = 0.0;
  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const double desired = model.desiredStock[i];
    const double base = desired * (0.35 + 0.55 * rng.nextDouble());
    s.inventory[i] = std::min(model.capacity[i], std::max(0.0, base));
  }
  return s;
}

double midPrice(const StationEconomyState& state, const StationEconomyModel& model, CommodityId id) {
  const auto& def = commodityDef(id);
  const std::size_t i = idx(id);

  const double desired = model.desiredStock[i];
  const double inv = state.inventory[i];

  if (desired <= 1e-9) return def.basePrice;

  const double ratio = inv / desired; // 1.0 means "normal"
  double factor = 1.0 + model.priceVolatility * (1.0 - ratio);

  // Scarcity spike if nearly empty
  if (inv < desired * 0.05) factor *= 1.4;

  factor = stellar::math::clamp(factor, 0.2, 5.0);
  return def.basePrice * factor;
}

void updateEconomyTo(StationEconomyState& state,
                     const StationEconomyModel& model,
                     double timeDays,
                     core::SplitMix64& rng,
                     double sampleIntervalDays) {
  if (timeDays <= state.lastUpdateDay) return;

  // Defensive clamp: avoid weird sampling behavior / infinite loops.
  sampleIntervalDays = std::max(1e-6, sampleIntervalDays);

  // Economy shocks need to be deterministic and *independent of update cadence*.
  // We do this by deriving a stable per-station base seed, then using an integer
  // day index to generate a fixed shock value for each commodity.
  //
  // NOTE: Universe::stationEconomy() currently re-seeds the RNG each call, so we
  // must not rely on RNG state continuity across calls.
  const core::u64 baseSeed = rng.nextU64();
  const double shockVol = std::max(0.0, model.shockVolatility);

  double t = state.lastUpdateDay;

  // Step at most one day at a time so our day-indexed shocks remain stable.
  while (t < timeDays) {
    const double day = std::floor(t);
    const double nextDay = day + 1.0;
    const double tNext = std::min(timeDays, nextDay);
    const double step = std::max(0.0, tNext - t);

    // Pre-mix a day seed once per loop.
    const core::u64 daySeed = core::hashCombine(baseSeed, static_cast<core::u64>(day));

    for (std::size_t i = 0; i < kCommodityCount; ++i) {
      double net = model.productionPerDay[i] - model.consumptionPerDay[i];

      if (shockVol > 0.0) {
        // Deterministic daily shock in [-1, 1].
        // Scaled by desired stock so all commodities get comparable *relative* drift.
        core::SplitMix64 srng(core::hashCombine(daySeed, static_cast<core::u64>(i)));
        const double noise = srng.nextDouble() * 2.0 - 1.0;
        const double desired = std::max(1e-9, model.desiredStock[i]);
        const double shockUnitsPerDay = noise * shockVol * desired;
        net += shockUnitsPerDay;
      }

      state.inventory[i] += net * step;
      state.inventory[i] = std::max(0.0, std::min(state.inventory[i], model.capacity[i]));
    }

    t = tNext;

    if (t - state.lastSampleDay >= sampleIntervalDays) {
      for (std::size_t i = 0; i < kCommodityCount; ++i) {
        const CommodityId id = static_cast<CommodityId>(i);
        const double p = midPrice(state, model, id);

        auto& hist = state.history[i];
        hist.push_back(PricePoint{t, p});
        if (hist.size() > 128) hist.erase(hist.begin(), hist.begin() + (hist.size() - 128));
      }
      state.lastSampleDay = t;
    }
  }

  state.lastUpdateDay = timeDays;
}

} // namespace stellar::econ
