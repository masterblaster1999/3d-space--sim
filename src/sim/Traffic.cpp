#include "stellar/sim/Traffic.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/econ/Market.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace stellar::sim {
namespace {

static constexpr std::size_t idx(econ::CommodityId id) { return static_cast<std::size_t>(id); }

struct BestCommodity {
  econ::CommodityId id{econ::CommodityId::Food};
  double score{0.0};
};

BestCommodity pickCommodity(const Station& src,
                           const econ::StationEconomyState& srcState,
                           const Station& dst,
                           const econ::StationEconomyState& dstState,
                           core::SplitMix64& rng) {
  std::array<BestCommodity, 3> top{};
  for (auto& t : top) { t.score = 0.0; }

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const econ::CommodityId cid = static_cast<econ::CommodityId>(i);

    const double prodS = src.economyModel.productionPerDay[i];
    const double consS = src.economyModel.consumptionPerDay[i];
    const double prodD = dst.economyModel.productionPerDay[i];
    const double consD = dst.economyModel.consumptionPerDay[i];

    const double netSrc = prodS - consS;
    const double netNeed = consD - prodD;
    if (netSrc <= 0.0 || netNeed <= 0.0) continue;

    const double invS = srcState.inventory[i];
    const double invD = dstState.inventory[i];
    const double desiredS = std::max(1e-9, src.economyModel.desiredStock[i]);
    const double desiredD = std::max(1e-9, dst.economyModel.desiredStock[i]);

    // If source is extremely dry, don't ship it.
    if (invS < desiredS * 0.22) continue;

    const double surplus = std::max(0.0, (invS - desiredS) / desiredS);
    const double shortage = std::max(0.0, (desiredD - invD) / desiredD);

    // Score favors natural producer->consumer flows, but is modulated by current
    // surplus/shortage so inventories can influence routes.
    double score = (netSrc * netNeed);
    score *= (0.55 + 1.25 * std::min(2.0, surplus));
    score *= (0.55 + 1.25 * std::min(2.0, shortage));

    // Small randomness so it doesn't collapse to a single lane.
    score *= (0.92 + 0.16 * rng.nextDouble());

    // Insert into a tiny "top N" list.
    if (score > top[0].score) {
      top[2] = top[1];
      top[1] = top[0];
      top[0] = {cid, score};
    } else if (score > top[1].score) {
      top[2] = top[1];
      top[1] = {cid, score};
    } else if (score > top[2].score) {
      top[2] = {cid, score};
    }
  }

  // If we found nothing producer->consumer, fall back to "anything with inventory".
  if (top[0].score <= 1e-9) {
    // Try a handful of random commodities.
    for (int tries = 0; tries < 6; ++tries) {
      const econ::CommodityId cid = static_cast<econ::CommodityId>(rng.range<int>(0, (int)econ::kCommodityCount - 1));
      const std::size_t i = idx(cid);
      const double invS = srcState.inventory[i];
      if (invS <= 1e-6) continue;
      top[0] = {cid, 1.0};
      break;
    }
  }

  // Choose among the top 3 with weighted probability (keeps some variety).
  const double s0 = top[0].score;
  const double s1 = top[1].score;
  const double s2 = top[2].score;
  const double sum = s0 + s1 + s2;
  if (sum > 1e-9) {
    const double r = rng.nextDouble() * sum;
    if (r < s0) return top[0];
    if (r < s0 + s1) return top[1];
    return top[2];
  }

  return top[0];
}

} // namespace

void simulateNpcTradeTraffic(Universe& universe,
                             const StarSystem& system,
                             double timeDays,
                             std::unordered_map<SystemId, int>& lastTrafficDayBySystem,
                             int kMaxBackfillDays) {
  if (system.stations.size() < 2) return;
  if (timeDays < 0.0) return;

  const SystemId sysId = system.stub.id;
  const int currentDay = (int)std::floor(timeDays);

  auto it = lastTrafficDayBySystem.find(sysId);
  if (it == lastTrafficDayBySystem.end()) {
    // First time we see this system: don't backfill the whole universe.
    // If we're already far into the timeline, do a tiny "warm start" (simulate one day)
    // so markets aren't always perfectly "fresh" on first arrival.
    lastTrafficDayBySystem.emplace(sysId, std::max(-1, currentDay - 1));
    it = lastTrafficDayBySystem.find(sysId);
  }

  int lastDay = it->second;
  if (currentDay <= lastDay) {
    // Still advance station economies to the current time so callers see up-to-date state.
    for (const auto& st : system.stations) (void)universe.stationEconomy(st, timeDays);
    return;
  }

  kMaxBackfillDays = std::max(1, kMaxBackfillDays);
  int startDay = lastDay + 1;
  if (currentDay - startDay + 1 > kMaxBackfillDays) {
    startDay = currentDay - kMaxBackfillDays + 1;
  }

  // Build stable station pointers once.
  std::vector<const Station*> stations;
  stations.reserve(system.stations.size());
  for (const auto& st : system.stations) stations.push_back(&st);

  for (int day = startDay; day <= currentDay; ++day) {
    // Seed per (universe, system, day)
    core::u64 s = core::hashCombine(universe.seed(), static_cast<core::u64>(sysId));
    s = core::hashCombine(s, core::seedFromText("traffic"));
    s = core::hashCombine(s, static_cast<core::u64>(day));
    core::SplitMix64 rng(s);

    // Update station economies to a stable time within the day.
    const double t = (double)day + 0.5;
    std::vector<econ::StationEconomyState*> states;
    states.reserve(stations.size());
    for (const Station* st : stations) {
      states.push_back(&universe.stationEconomy(*st, t));
    }

    const int n = (int)stations.size();
    const int baseRuns = std::clamp(1 + n, 2, 12);
    const int runs = rng.range<int>(std::max(1, baseRuns - 1), std::min(16, baseRuns + 2));

    for (int r = 0; r < runs; ++r) {
      int srcIdx = rng.range<int>(0, n - 1);
      int dstIdx = rng.range<int>(0, n - 1);
      if (n > 1) {
        int guard = 0;
        while (dstIdx == srcIdx && guard++ < 6) dstIdx = rng.range<int>(0, n - 1);
      }
      if (dstIdx == srcIdx) continue;

      const Station& src = *stations[srcIdx];
      const Station& dst = *stations[dstIdx];
      econ::StationEconomyState& srcState = *states[srcIdx];
      econ::StationEconomyState& dstState = *states[dstIdx];

      const BestCommodity pick = pickCommodity(src, srcState, dst, dstState, rng);
      const econ::CommodityId cid = pick.id;
      const std::size_t i = idx(cid);

      const double invS = srcState.inventory[i];
      const double invD = dstState.inventory[i];
      const double desiredS = std::max(1e-9, src.economyModel.desiredStock[i]);
      const double desiredD = std::max(1e-9, dst.economyModel.desiredStock[i]);
      const double capD = std::max(1e-9, dst.economyModel.capacity[i]);

      // Compute a "reasonable" shipment size.
      const double surplus = std::max(0.0, invS - desiredS * 0.25);
      const double need = std::max(0.0, (desiredD * 1.05) - invD);
      if (surplus <= 1e-6 || need <= 1e-6) continue;

      double units = std::min(surplus, need);
      units = std::min(units, desiredS * 0.45);
      units = std::min(units, capD * 0.30);
      units *= rng.range<double>(0.25, 0.85);

      if (units <= 1e-4) continue;

      const double taken = econ::takeInventory(srcState, src.economyModel, cid, units);
      if (taken <= 1e-6) continue;
      (void)econ::addInventory(dstState, dst.economyModel, cid, taken);
    }
  }

  // Advance all affected stations to the actual current time.
  for (const auto& st : system.stations) (void)universe.stationEconomy(st, timeDays);

  it->second = currentDay;
}

} // namespace stellar::sim
