#pragma once

#include "stellar/sim/SaveGame.h"

#include <cstddef>
#include <string>

namespace stellar::sim {

class Universe;
struct StarSystem;
struct Station;

struct MissionBoardParams {
  double searchRadiusLy{160.0};
  std::size_t maxCandidates{128};

  int baseOfferCount{6};
  int minOfferCount{4};
  int maxOfferCount{9};

  // Rep thresholds that add one extra offer each.
  double repThreshold1{25.0};
  double repThreshold2{50.0};

  // Mission type weights (must sum <= 1.0; remainder goes to BountyKill).
  double wCourier{0.31};
  double wDelivery{0.31};
  double wMultiDelivery{0.13};
  double wPassenger{0.11};
  double wSmuggle{0.04};
  double wBountyScan{0.06};

  // Reward scaling (positive rep only): reward *= 1 + (rep/100)*repRewardBonus
  double repRewardBonus{0.10};
};

struct MissionTickResult {
  int failed{0};
};

struct MissionDockResult {
  int completed{0};
  int progressedMultiLeg{0};
};

// Ensure SaveGame::missionOffers is populated for the given docked station and day.
// Offers are deterministic per (station id, dayStamp, faction id) and persisted in the save.
void refreshMissionOffers(Universe& universe,
                          const StarSystem& currentSystem,
                          const Station& dockedStation,
                          double timeDays,
                          double playerRep,
                          SaveGame& ioSave,
                          const MissionBoardParams& params = {});

// Accept an offer from SaveGame::missionOffers.
// On success:
//  - consumes station inventory / player credits as required
//  - adds cargo to player if needed
//  - moves offer into SaveGame::missions and removes it from missionOffers
bool acceptMissionOffer(Universe& universe,
                        const Station& dockedStation,
                        double timeDays,
                        double playerRep,
                        SaveGame& ioSave,
                        std::size_t offerIndex,
                        std::string* outError = nullptr);

// Mark missions as failed when deadlines pass (applies a small rep hit).
MissionTickResult tickMissionDeadlines(SaveGame& ioSave,
                                       double timeDays,
                                       double repPenaltyOnFail = -4.0);

// When docked, attempt to complete Courier/Delivery/MultiDelivery missions.
// Applies rewards, rep gains, and returns delivered cargo to the destination market.
MissionDockResult tryCompleteMissionsAtDock(Universe& universe,
                                            const StarSystem& currentSystem,
                                            const Station& dockedStation,
                                            double timeDays,
                                            SaveGame& ioSave,
                                            double repRewardOnComplete = +2.0);

} // namespace stellar::sim
