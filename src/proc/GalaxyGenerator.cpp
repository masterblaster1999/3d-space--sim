#include "stellar/proc/GalaxyGenerator.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/proc/NameGenerator.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

std::size_t SectorCoordHash::operator()(const SectorCoord& c) const noexcept {
  core::u64 h = 0xcbf29ce484222325ull;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(c.x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(c.y)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(c.z)));
  return static_cast<std::size_t>(h);
}

GalaxyGenerator::GalaxyGenerator(core::u64 seed, GalaxyParams params)
: seed_(seed), params_(params) {}

SectorCoord GalaxyGenerator::sectorOf(const math::Vec3d& posLy) const {
  const double s = params_.sectorSizeLy;
  return {
    static_cast<core::i32>(std::floor(posLy.x / s)),
    static_cast<core::i32>(std::floor(posLy.y / s)),
    static_cast<core::i32>(std::floor(posLy.z / s)),
  };
}

sim::SystemId GalaxyGenerator::makeSystemId(const SectorCoord& coord, core::u32 localIndex) const {
  // Pack the sector coordinate + local index into 64 bits:
  //   [x:16][y:16][z:16][i:16]
  //
  // This allows Universe::getSystem(id) to decode the sector directly from the id
  // and stream/generate the correct sector on-demand (without requiring a hint stub).
  //
  // NOTE: Galaxy radius/sectorSize keep coords comfortably within +/- 32k sectors
  // (default: 50k ly / 10 ly = 5k sectors), so a 16-bit biased encoding is safe.
  if (coord.x < -32768 || coord.x > 32767 ||
      coord.y < -32768 || coord.y > 32767 ||
      coord.z < -32768 || coord.z > 32767) {
    // Extremely out-of-range query; fall back to a hashed id (won't be decodable).
    core::u64 h = seed_;
    h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.x)));
    h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.y)));
    h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.z)));
    h = core::hashCombine(h, static_cast<core::u64>(localIndex));
    return static_cast<sim::SystemId>(h);
  }

  const auto bias = [](core::i32 v) -> core::u16 {
    return static_cast<core::u16>(static_cast<core::i32>(v) + 32768);
  };

  const core::u64 bx = static_cast<core::u64>(bias(coord.x));
  const core::u64 by = static_cast<core::u64>(bias(coord.y));
  const core::u64 bz = static_cast<core::u64>(bias(coord.z));
  const core::u64 bi = static_cast<core::u64>(static_cast<core::u16>(localIndex & 0xFFFFu));

  const core::u64 id = (bx << 48) | (by << 32) | (bz << 16) | (bi << 0);
  return static_cast<sim::SystemId>(id);
}

static int poisson(core::SplitMix64& rng, double mean) {
  if (mean <= 0.0) return 0;
  if (mean < 30.0) {
    // Knuth
    const double L = std::exp(-mean);
    int k = 0;
    double p = 1.0;
    do {
      ++k;
      p *= rng.nextDouble();
    } while (p > L);
    return k - 1;
  }

  // Normal approximation for larger means (Box-Muller)
  const double u1 = std::max(1e-12, rng.nextDouble());
  const double u2 = rng.nextDouble();
  const double z0 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * stellar::math::kPi * u2);
  const double k = mean + z0 * std::sqrt(mean);
  return std::max(0, static_cast<int>(std::round(k)));
}

static sim::StarClass pickStarClass(core::SplitMix64& rng) {
  const double r = rng.nextDouble();
  // Very rough main-sequence-ish distribution (not astrophysically accurate).
  if (r < 0.0003) return sim::StarClass::O;
  if (r < 0.0016) return sim::StarClass::B;
  if (r < 0.006)  return sim::StarClass::A;
  if (r < 0.03)   return sim::StarClass::F;
  if (r < 0.10)   return sim::StarClass::G;
  if (r < 0.30)   return sim::StarClass::K;
  return sim::StarClass::M;
}

static core::u32 pickFaction(const math::Vec3d& posLy, const std::vector<sim::Faction>& factions) {
  // 0 is Independent
  core::u32 best = 0;
  double bestD = 1e30;

  for (const auto& f : factions) {
    if (f.id == 0) continue;
    const double d = (posLy - f.homePosLy).length();
    if (d < f.influenceRadiusLy && d < bestD) {
      bestD = d;
      best = f.id;
    }
  }
  return best;
}

Sector GalaxyGenerator::generateSector(const SectorCoord& coord, const std::vector<sim::Faction>& factions) const {
  // Sector seed
  core::u64 h = seed_;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.y)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.z)));

  core::SplitMix64 rng(h);

  Sector sector{};
  sector.coord = coord;

  // Density based on distance from galactic center + vertical falloff.
  const double s = params_.sectorSizeLy;
  const math::Vec3d centerLy{
    (coord.x + 0.5) * s,
    (coord.y + 0.5) * s,
    (coord.z + 0.5) * s
  };
  const double rxy = std::sqrt(centerLy.x*centerLy.x + centerLy.y*centerLy.y);
  const double z = std::abs(centerLy.z);

  const double radial = std::exp(-rxy / std::max(1.0, params_.radialScaleLengthLy));
  const double vertical = std::exp(-z / std::max(1.0, params_.thicknessLy * 0.5));
  const double mean = params_.baseMeanSystemsPerSector * radial * vertical;

  const int n = poisson(rng, mean);

  sector.systems.reserve(static_cast<std::size_t>(std::max(0, n)));

  NameGenerator ng{};
  const double halfThickness = params_.thicknessLy * 0.5;

  for (int i = 0; i < n; ++i) {
    // Try a few candidates to keep within disc bounds.
    bool placed = false;
    math::Vec3d pos{};
    for (int tries = 0; tries < 8 && !placed; ++tries) {
      const double ux = rng.nextDouble();
      const double uy = rng.nextDouble();
      const double uz = rng.nextDouble();
      pos = {
        (coord.x + ux) * s,
        (coord.y + uy) * s,
        (coord.z + uz) * s
      };

      const double rr = std::sqrt(pos.x*pos.x + pos.y*pos.y);
      if (rr > params_.radiusLy) continue;
      if (std::abs(pos.z) > halfThickness) continue;
      placed = true;
    }
    if (!placed) continue;

    sim::SystemStub stub{};
    stub.id = makeSystemId(coord, static_cast<core::u32>(i));
    stub.seed = core::hashCombine(seed_, static_cast<core::u64>(stub.id));

    ng.reseed(stub.seed);
    stub.name = ng.systemName();
    stub.posLy = pos;
    stub.primaryClass = pickStarClass(rng);
    stub.planetCount = rng.range(0, 12);
    stub.stationCount = std::max(1, rng.range(0, 3)); // ensure at least 1 station for gameplay
    stub.factionId = pickFaction(stub.posLy, factions);

    sector.systems.push_back(std::move(stub));
  }

  // Stable order for deterministic query results.
  std::sort(sector.systems.begin(), sector.systems.end(), [](const sim::SystemStub& a, const sim::SystemStub& b) {
    return a.id < b.id;
  });

  return sector;
}

} // namespace stellar::proc
