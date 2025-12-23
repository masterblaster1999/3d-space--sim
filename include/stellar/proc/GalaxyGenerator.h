#pragma once

#include "stellar/core/Random.h"
#include "stellar/math/Vec3.h"

#include <cstddef>
#include <vector>

namespace stellar::proc {

struct GalaxySystemStub {
  stellar::core::u64 id = 0;
  stellar::math::Vec3d positionLy{};
};

struct GalaxyGenConfig {
  stellar::core::u64 seed = 1;
  std::size_t systemCount = 1000;

  double radiusLy = 50000.0;
  double thicknessLy = 2000.0;
};

class GalaxyGenerator {
public:
  explicit GalaxyGenerator(GalaxyGenConfig cfg) : m_cfg(cfg) {}

  std::vector<GalaxySystemStub> generate() const;

private:
  GalaxyGenConfig m_cfg{};
};

} // namespace stellar::proc
