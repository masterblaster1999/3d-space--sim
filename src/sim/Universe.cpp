#include "stellar/sim/Universe.h"

#include "stellar/core/Log.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/proc/SystemGenerator.h"

#include <sstream>

namespace stellar::sim {

Universe::Universe(UniverseConfig cfg) : m_cfg(cfg) {}

void Universe::generate() {
  std::ostringstream oss;
  oss << "Generating universe: seed=" << m_cfg.seed
      << ", systems=" << m_cfg.systemCount
      << ", radiusLy=" << m_cfg.radiusLy
      << ", thicknessLy=" << m_cfg.thicknessLy;
  STELLAR_LOG_INFO(oss.str());

  proc::GalaxyGenerator galaxy(proc::GalaxyGenConfig{
    .seed = m_cfg.seed,
    .systemCount = m_cfg.systemCount,
    .radiusLy = m_cfg.radiusLy,
    .thicknessLy = m_cfg.thicknessLy,
  });

  const auto stubs = galaxy.generate();

  proc::SystemGenerator sysgen(proc::SystemGenConfig{ .galaxySeed = m_cfg.seed });

  m_systems.clear();
  m_systems.reserve(stubs.size());

  for (const auto& stub : stubs) {
    m_systems.push_back(sysgen.generate(stub));
  }

  STELLAR_LOG_INFO("Universe generation complete.");
}

} // namespace stellar::sim
