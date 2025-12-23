#pragma once

#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/proc/NameGenerator.h"
#include "stellar/sim/Celestial.h"

namespace stellar::proc {

struct SystemGenConfig {
  stellar::core::u64 galaxySeed = 1;
};

class SystemGenerator {
public:
  explicit SystemGenerator(SystemGenConfig cfg);

  stellar::sim::StarSystem generate(const GalaxySystemStub& stub) const;

private:
  SystemGenConfig m_cfg{};
  NameGenerator m_names;
};

} // namespace stellar::proc
