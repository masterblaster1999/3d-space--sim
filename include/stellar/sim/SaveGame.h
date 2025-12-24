#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"

#include <array>
#include <string>
#include <vector>

namespace stellar::sim {

struct StationEconomyOverride {
  StationId stationId{0};
  econ::StationEconomyState state{};
};

// Simple reputation record (per faction id).
// Stored in the save file so mission progression / market fees can persist.
struct FactionReputation {
  core::u32 factionId{0};
  double rep{0.0}; // roughly [-100,+100]
};

// Lightweight "gameplay" mission representation.
// Stored in the save file so early progression loops (cargo delivery/courier/bounties)
// persist across runs.
enum class MissionType : core::u8 {
  Courier = 0,
  Delivery,
  BountyScan,
  BountyKill,
  MultiDelivery,
};

struct Mission {
  core::u64 id{0};
  MissionType type{MissionType::Courier};

  // Issuing faction (used for reputation effects).
  core::u32 factionId{0};

  SystemId fromSystem{0};
  StationId fromStation{0};

  SystemId toSystem{0};
  StationId toStation{0};

  // Delivery missions use a commodity + units.
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // Multi-hop deliveries: optional intermediate stop.
  // If viaSystem/viaStation are non-zero, the player must first deliver to the via stop,
  // then proceed to the final destination (toSystem/toStation).
  SystemId viaSystem{0};
  StationId viaStation{0};
  core::u8 leg{0}; // 0 = heading to via (if any) else final; 1 = heading to final

  // Bounty scan missions can attach to a specific NPC id (spawned in the target system).
  // 0 means "no specific target".
  core::u64 targetNpcId{0};

  // Bounty scan progress flag.
  bool scanned{false};

  double reward{0.0};
  double deadlineDay{0.0};

  bool completed{false};
  bool failed{false};

  // If true, the station provided the cargo at acceptance time (taking from station inventory).
  // If false, the player must source the cargo themselves.
  bool cargoProvided{true};
};

struct SaveGame {
  int version{3};

  core::u64 seed{0};
  double timeDays{0.0};

  SystemId currentSystem{0};
  StationId dockedStation{0};

  // Player ship
  math::Vec3d shipPosKm{0,0,0};
  math::Vec3d shipVelKmS{0,0,0};
  math::Quatd shipOrient{1,0,0,0};
  math::Vec3d shipAngVelRadS{0,0,0};

  // Economy
  double credits{1000.0};
  std::array<double, econ::kCommodityCount> cargo{}; // units

  // Ship meta/progression
  double fuel{45.0};
  double fuelMax{45.0};
  double fsdRangeLy{18.0};
  double hull{1.0}; // 0..1
  double cargoCapacityKg{420.0};
  double fsdReadyDay{0.0}; // timeDays when the next hyperspace jump is allowed

  // Missions
  core::u64 nextMissionId{1};
  std::vector<Mission> missions{};

  // Reputation
  std::vector<FactionReputation> reputation{};

  std::vector<StationEconomyOverride> stationOverrides{};
};

bool saveToFile(const SaveGame& s, const std::string& path);
bool loadFromFile(const std::string& path, SaveGame& out);

} // namespace stellar::sim
