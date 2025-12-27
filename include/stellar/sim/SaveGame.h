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

// Simple bounty record (per faction id).
// If non-zero, local police will treat the player as WANTED in that faction's space.
struct FactionBounty {
  core::u32 factionId{0};
  double bountyCr{0.0};
};

// Persistent day-stamp for deterministic "background" traffic simulation.
// Only systems the player has visited need to be tracked.
struct SystemTrafficStamp {
  SystemId systemId{0};
  int dayStamp{-1};
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
  Passenger,
  Smuggle,
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

  // Delivery / smuggling missions use a commodity + units.
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
  int version{12};

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

  // Exploration
  double explorationDataCr{0.0};
  std::vector<core::u64> scannedKeys{};

  // Ship meta/progression
  double fuel{45.0};
  double fuelMax{45.0};
  double fsdRangeLy{18.0};
  double hull{1.0}; // 0..1
  double shield{1.0}; // 0..1
  double heat{0.0}; // gameplay heat (0..~120)
  double cargoCapacityKg{420.0};

  // Passenger capacity (seats) for cabin-style missions.
  // Kept intentionally simple for the prototype loop.
  int passengerSeats{2};
  double fsdReadyDay{0.0}; // timeDays when the next hyperspace jump is allowed

  // Loadout / progression (kept simple for now: small ints, interpreted by gameplay code).
  // These are *not* physics-critical; they tune HUD/combat feel and basic progression loops.
  core::u8 shipHull{0};        // 0 = Scout (starter), 1 = Hauler, 2 = Fighter
  core::u8 thrusterMk{1};      // 1..3
  core::u8 shieldMk{1};        // 1..3
  core::u8 distributorMk{1};   // 1..3
  core::u8 weaponPrimary{0};   // enum in gameplay (0=beam, 1=pulse, 2=cannon, 3=rail)
  core::u8 weaponSecondary{2}; // default cannon

  // Smuggling / stealth: reduces chance of cargo scans when carrying contraband.
  core::u8 smuggleHoldMk{0}; // 0..3

  // Missions
  core::u64 nextMissionId{1};
  std::vector<Mission> missions{};

  // Mission board (cached offers) - persisted so boards don't reroll on UI refresh / reload.
  StationId missionOffersStationId{0};
  int missionOffersDayStamp{-1};
  std::vector<Mission> missionOffers{};

  // Reputation
  std::vector<FactionReputation> reputation{};

  // Law / bounties
  std::vector<FactionBounty> bounties{};
  // Bounty vouchers earned for destroying criminals (redeem at stations).
  std::vector<FactionBounty> bountyVouchers{};

  // Background NPC traffic simulation (for markets). Stores last simulated day per system.
  std::vector<SystemTrafficStamp> trafficStamps{};

  std::vector<StationEconomyOverride> stationOverrides{};
};

bool saveToFile(const SaveGame& s, const std::string& path);
bool loadFromFile(const std::string& path, SaveGame& out);

} // namespace stellar::sim
