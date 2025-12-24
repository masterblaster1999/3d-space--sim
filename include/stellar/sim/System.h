#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/sim/Celestial.h"

#include <vector>

namespace stellar::sim {

struct Station {
  StationId id{0};
  std::string name;
  econ::StationType type{econ::StationType::Outpost};
  core::u32 factionId{0};       // 0 = independent
  double feeRate{0.0};          // market fee rate (0..1)
  econ::StationEconomyModel economyModel{};

  // Physical placement (around primary star)
  OrbitElements orbit{};

  // Size (km) - used for rendering scale and collision/docking volumes.
  double radiusKm{6000.0};

  // Docking / approach ("mail-slot" style)
  double commsRangeKm{120000.0};        // range to request clearance (km)
  double approachLengthKm{80000.0};     // approach corridor length (km), extending outward from slot
  double approachRadiusKm{12000.0};     // cylindrical approach radius (km)
  double maxApproachSpeedKmS{0.20};     // max relative speed in corridor (km/s)

  // Slot/tunnel dimensions in station-local space (km).
  // Station-local +Z points outward from the slot.
  double slotWidthKm{6000.0};
  double slotHeightKm{2500.0};
  double slotDepthKm{9000.0};           // depth of the "tunnel" from entrance into the hangar
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
