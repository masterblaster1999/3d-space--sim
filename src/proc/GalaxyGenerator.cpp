#include "stellar/proc/GalaxyGenerator.h"

#include "stellar/math/Math.h"

#include <cmath>

namespace stellar::proc {

std::vector<GalaxySystemStub> GalaxyGenerator::generate() const {
  std::vector<GalaxySystemStub> out;
  out.reserve(m_cfg.systemCount);

  for (std::size_t i = 0; i < m_cfg.systemCount; ++i) {
    const auto id = stellar::core::deriveSeed(m_cfg.seed, static_cast<stellar::core::u64>(i));

    // Per-index RNG: adding more systems won't change existing ones.
    stellar::core::SplitMix64 rng(id);

    const double u = rng.nextDouble01();
    const double v = rng.nextDouble01();

    const double r = m_cfg.radiusLy * std::sqrt(u); // sqrt for uniform disc density
    const double theta = stellar::math::twoPi * v;

    const double x = r * std::cos(theta);
    const double z = r * std::sin(theta);

    const double y = rng.uniform(-0.5, 0.5) * m_cfg.thicknessLy;

    out.push_back(GalaxySystemStub{
      .id = id,
      .positionLy = {x, y, z},
    });
  }

  return out;
}

} // namespace stellar::proc
