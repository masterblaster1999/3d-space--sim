#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Celestial.h"

#include <vector>

namespace stellar::sim {

struct UniverseConfig {
  stellar::core::u64 seed = 1;

  // How many star systems to generate for this sample universe.
  // (You can later replace this with on-demand streaming.)
  std::size_t systemCount = 1000;

  // Simple disc galaxy shape in light-years
  double radiusLy = 50000.0;
  double thicknessLy = 2000.0;
};

class Universe {
public:
  explicit Universe(UniverseConfig cfg);

  // Generates the full galaxy + all systems.
  void generate();

  const UniverseConfig& config() const { return m_cfg; }

  std::size_t size() const { return m_systems.size(); }
  const std::vector<StarSystem>& systems() const { return m_systems; }
  const StarSystem& system(std::size_t index) const { return m_systems.at(index); }

private:
  UniverseConfig m_cfg{};
  std::vector<StarSystem> m_systems;
};

} // namespace stellar::sim
