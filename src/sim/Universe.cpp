#include "stellar/sim/Universe.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/proc/NameGenerator.h"
#include "stellar/proc/SystemGenerator.h"
#include "stellar/econ/Economy.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace stellar::sim {

static proc::SectorCoord decodeSector(SystemId id, core::u32& outLocalIndex) {
  const auto unbias = [](core::u16 b) -> core::i32 { return static_cast<core::i32>(b) - 32768; };

  const core::u16 bx = static_cast<core::u16>((id >> 48) & 0xFFFFu);
  const core::u16 by = static_cast<core::u16>((id >> 32) & 0xFFFFu);
  const core::u16 bz = static_cast<core::u16>((id >> 16) & 0xFFFFu);
  const core::u16 bi = static_cast<core::u16>((id >> 0)  & 0xFFFFu);

  outLocalIndex = static_cast<core::u32>(bi);
  return proc::SectorCoord{ unbias(bx), unbias(by), unbias(bz) };
}

static StarClass pickStarClass(core::SplitMix64& rng) {
  const double r = rng.nextDouble();
  // Very rough main-sequence-ish distribution.
  if (r < 0.0003) return StarClass::O;
  if (r < 0.0016) return StarClass::B;
  if (r < 0.006)  return StarClass::A;
  if (r < 0.03)   return StarClass::F;
  if (r < 0.10)   return StarClass::G;
  if (r < 0.30)   return StarClass::K;
  return StarClass::M;
}

Universe::Universe(core::u64 seed, proc::GalaxyParams params)
: seed_(seed),
  galaxyParams_(params),
  galaxyGen_(seed, params) {
  factions_ = generateFactions(seed_, 8);
}

void Universe::setCacheCaps(std::size_t sectorCap, std::size_t systemCap, std::size_t stationCap) {
  sectorCache_.setCapacity(sectorCap);
  systemCache_.setCapacity(systemCap);
  stationEconomyCache_.setCapacity(stationCap);
}

const proc::Sector& Universe::sector(const proc::SectorCoord& coord) {
  if (auto* cached = sectorCache_.get(coord)) return *cached;
  auto sec = galaxyGen_.generateSector(coord, factions_);
  return sectorCache_.put(coord, std::move(sec));
}

std::vector<SystemStub> Universe::queryNearby(const math::Vec3d& posLy,
                                              double radiusLy,
                                              std::size_t maxResults) {
  std::vector<SystemStub> out;
  if (radiusLy <= 0.0) return out;

  const double r2 = radiusLy * radiusLy;
  const double s = galaxyParams_.sectorSizeLy;

  const proc::SectorCoord minC{
    static_cast<core::i32>(std::floor((posLy.x - radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.y - radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.z - radiusLy) / s)),
  };
  const proc::SectorCoord maxC{
    static_cast<core::i32>(std::floor((posLy.x + radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.y + radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.z + radiusLy) / s)),
  };

  struct Item { SystemStub stub; double d2; };
  std::vector<Item> items;

  for (core::i32 x = minC.x; x <= maxC.x; ++x) {
    for (core::i32 y = minC.y; y <= maxC.y; ++y) {
      for (core::i32 z = minC.z; z <= maxC.z; ++z) {
        const proc::SectorCoord c{x,y,z};
        const proc::Sector& sec = sector(c);

        for (const auto& stub : sec.systems) {
          const math::Vec3d d = stub.posLy - posLy;
          const double dd = d.lengthSq();
          if (dd <= r2) items.push_back(Item{stub, dd});
        }
      }
    }
  }

  std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
    if (a.d2 != b.d2) return a.d2 < b.d2;
    return a.stub.id < b.stub.id;
  });

  if (items.size() > maxResults) items.resize(maxResults);

  out.reserve(items.size());
  for (auto& it : items) out.push_back(std::move(it.stub));
  return out;
}

std::optional<SystemStub> Universe::findClosestSystem(const math::Vec3d& posLy, double maxRadiusLy) {
  auto list = queryNearby(posLy, maxRadiusLy, 64);
  if (list.empty()) return std::nullopt;

  // queryNearby is distance-sorted already
  return list.front();
}

const StarSystem& Universe::getSystem(SystemId id, const SystemStub* hintStub) {
  if (auto* cached = systemCache_.get(id)) return *cached;

  SystemStub stub{};
  bool haveStub = false;

  if (hintStub && hintStub->id == id) {
    stub = *hintStub;
    haveStub = true;
  } else {
    // Decode sector from id and search it.
    core::u32 localIndex = 0;
    const proc::SectorCoord c = decodeSector(id, localIndex);

    const proc::Sector& sec = sector(c);

    auto it = std::lower_bound(sec.systems.begin(), sec.systems.end(), id,
                               [](const SystemStub& s, SystemId idv) { return s.id < idv; });
    if (it != sec.systems.end() && it->id == id) {
      stub = *it;
      haveStub = true;
    } else {
      std::ostringstream oss;
      oss << "Universe::getSystem: stub not found for id=" << id
          << " (sector " << c.x << "," << c.y << "," << c.z
          << " localIndex=" << localIndex << "). Generating fallback stub.";
      stellar::core::log(stellar::core::LogLevel::Warn, oss.str());

      // Fallback: deterministic from id + seed.
      //
      // We also synthesize a plausible disc position so legacy/unresolvable ids don't
      // all collapse to the origin (which can make debugging/nav UI confusing).
      stub.id = id;
      stub.seed = core::hashCombine(seed_, static_cast<core::u64>(id));
      proc::NameGenerator ng(stub.seed);
      stub.name = ng.systemName();

      core::SplitMix64 frng(core::hashCombine(stub.seed, core::seedFromText("fallback_stub")));

      const double r = galaxyParams_.radiusLy * std::sqrt(frng.nextDouble());
      const double a = frng.nextDouble() * 6.283185307179586;
      const double z = (frng.nextDouble() - 0.5) * galaxyParams_.thicknessLy;
      stub.posLy = { r * std::cos(a), r * std::sin(a), z };

      stub.primaryClass = pickStarClass(frng);
      stub.planetCount = frng.range(0, 12);
      stub.stationCount = std::max(1, frng.range(0, 3));
      stub.factionId = 0;
      haveStub = true;


    }
  }

  (void)haveStub;

  StarSystem sys = proc::generateSystem(stub, factions_);
  return systemCache_.put(id, std::move(sys));
}

econ::StationEconomyState& Universe::stationEconomy(const Station& station, double timeDays) {
  econ::StationEconomyState* st = stationEconomyCache_.get(station.id);
  if (!st) {
    core::SplitMix64 rng(core::hashCombine(seed_, static_cast<core::u64>(station.id)));
    econ::StationEconomyState init = econ::makeInitialState(station.economyModel, rng);
    st = &stationEconomyCache_.put(station.id, std::move(init));
  }

  core::SplitMix64 rng(core::hashCombine(seed_, static_cast<core::u64>(station.id)));
  econ::updateEconomyTo(*st, station.economyModel, timeDays, rng);
  return *st;
}

std::vector<StationEconomyOverride> Universe::exportStationOverrides() const {
  std::vector<StationEconomyOverride> out;

  auto snap = stationEconomyCache_.snapshot();
  out.reserve(snap.size());
  for (auto& kv : snap) {
    StationEconomyOverride ov{};
    ov.stationId = kv.first;
    ov.state = std::move(kv.second);
    out.push_back(std::move(ov));
  }

  return out;
}

void Universe::importStationOverrides(const std::vector<StationEconomyOverride>& overrides) {
  for (const auto& ov : overrides) {
    stationEconomyCache_.put(ov.stationId, ov.state);
  }
}

} // namespace stellar::sim
