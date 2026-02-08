#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Combat.h"
#include "stellar/sim/Countermeasures.h"

#include <cstddef>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// MissileDefense (headless helpers)
// -----------------------------------------------------------------------------
//
// Utilities for detecting inbound missiles in a deterministic, frame-rate
// agnostic way. This is used by the game app to give NPCs basic defensive
// behaviors (countermeasures + evasive thrust jinks) without baking the logic
// directly into the big game loop.

struct MissileThreatParams {
  // Cosine of the maximum "off boresight" angle for considering a missile
  // inbound. 1.0 means perfectly aligned, 0.0 means 90 degrees.
  double minApproachCos{0.2};

  // Ignore cases where relative closing speed is extremely low.
  double minClosingKmS{0.02};

  // Ignore missiles farther than this distance.
  double maxConsiderDistKm{250000.0};
};

struct MissileThreatSummary {
  bool inbound{false};
  MissileSeekerType seeker{MissileSeekerType::Heat};

  double distKm{0.0};
  double closingKmS{0.0};
  double ttiSec{0.0};
  double approachCos{0.0};

  std::size_t missileIndex{0};
  core::u64 shooterId{0};
  bool fromPlayer{false};
};

// Recommended evasion direction for an inbound missile (deterministic).
//
// The goal is to provide a *direction only*; the game can decide how strongly
// to apply thrust based on pilot skill, ship performance, and time-to-impact.
struct MissileEvasionPlan {
  bool valid{false};

  // Unit vector (world space) pointing in a good direction for applying
  // lateral thrust to increase the predicted miss distance.
  math::Vec3d dirWorld{0, 0, 0};

  // Time of closest approach under constant relative velocity (seconds).
  double tClosestSec{0.0};

  // Predicted miss distance at closest approach if the target does nothing (km).
  double missDistanceKm{0.0};

  // Relative closing speed along the line of sight (km/s). Positive means inbound.
  double closingKmS{0.0};
};

struct MissileEvasionParams {
  // Ignore cases where relative speed is extremely low (fallback direction is used).
  double minRelSpeedKmS{0.02};

  // Treat the closest-approach offset as degenerate below this threshold (km).
  double minMissVecKm{1e-3};

  // If true, project the output onto the plane perpendicular to the current
  // line-of-sight (produces a more "lateral jink" that tends to increase LOS rate).
  bool enforceLateralToLos{true};

  // --- Radar-specific evasion tuning ---
  //
  // When a radar seeker models a doppler notch (see Missile::radarDopplerNotchKmS),
  // targets can sometimes break lock by driving the geometry toward low range-rate.
  //
  // This struct optionally enables a simple deterministic "beaming" bias for NPCs:
  // blend the classic closest-approach "jink" direction with a direction that tends
  // to rotate the ship's velocity toward being perpendicular to the current LOS.
  //
  // NOTE:
  //   This is not a full aircraft energy model; it is a gameplay-oriented heuristic
  //   that interacts with the seeker notch mechanic in a stable way.

  // Enable the radar beaming bias when the inbound missile is a radar seeker.
  bool enableRadarBeaming{true};

  // Blend factor in [0,1] that controls how strongly the beaming direction is mixed
  // into the classic miss-distance jink.
  //
  // 0.0 -> pure jink (original behavior)
  // 1.0 -> pure beaming direction
  double radarBeamBlend{0.65};

  // Only engage the beaming bias when the absolute radial velocity exceeds
  //   radarBeamEngageNotchMultiple * missile.radarDopplerNotchKmS.
  //
  // This avoids over-steering when the target is already nearly notched.
  double radarBeamEngageNotchMultiple{1.10};
};


// -----------------------------------------------------------------------------
// Terrain masking planner (deterministic)
// -----------------------------------------------------------------------------
//
// Some seekers require a clear line-of-sight to maintain track (see
// Missile::requireLineOfSight and Missile::datalinkRequireLineOfSight).
//
// This helper gives NPCs a simple way to exploit that mechanic: pick a nearby
// asteroid and recommend a direction that moves the target toward a cover point
// "behind" the asteroid relative to the inbound missile.
//
// The output is intentionally only a *direction* — the game decides how much
// thrust to apply based on ship capability and pilot skill.

struct MissileTerrainMaskParams {
  // Only consider asteroids whose resulting cover point is within this distance
  // of the target (km).
  double maxCoverTravelKm{35000.0};

  // Predict the missile position this far into the future (seconds) when
  // computing the cover point direction (helps when missiles are fast).
  double lookaheadSec{0.65};

  // Extra padding added outside the asteroid radius when generating the cover
  // point (km). This effectively selects a point slightly outside the shadow
  // boundary so numerical LOS checks have margin.
  double coverPadKm{0.01};

  // If true and the target is already occluded from the missile by an asteroid,
  // return an invalid plan (NPC can keep its current maneuver).
  bool ignoreIfAlreadyOccluded{true};
};

struct MissileTerrainMaskPlan {
  bool valid{false};

  // Unit vector (world space) pointing toward the selected cover point.
  math::Vec3d dirWorld{0, 0, 0};

  // Chosen asteroid (0 when not applicable).
  core::u64 asteroidId{0};

  // A point "behind" the asteroid relative to the inbound missile, with
  // coverPadKm applied.
  math::Vec3d coverPointKm{0, 0, 0};

  // Distance from the target to the cover point (km).
  double travelKm{0.0};
};

// Compute a terrain-masking direction against a specific inbound missile.
//
// targets/targetCount should include asteroid SphereTargets (kind Asteroid).
// seed is used only for deterministic tie-breaking when multiple asteroids are
// equally good candidates.
MissileTerrainMaskPlan planMissileTerrainMasking(const Missile& missile,
                                                 const math::Vec3d& targetPosKm,
                                                 const math::Vec3d& targetVelKmS,
                                                 const SphereTarget* targets,
                                                 std::size_t targetCount,
                                                 core::u64 seed,
                                                 const MissileTerrainMaskParams& params = {});


// -----------------------------------------------------------------------------
// Countermeasure response planning (deterministic)
// -----------------------------------------------------------------------------
//
// A lightweight helper intended for NPCs:
//  - choose an appropriate countermeasure type based on seeker (flare vs chaff)
//  - schedule release timing to coincide with seeker measurement updates
//
// When missiles model discrete seeker updates (Missile::seekerUpdatePeriodSimSec),
// releasing countermeasures just before a measurement tick can materially improve
// effectiveness (the missile "sees" the decoy on the next scan).

struct CountermeasureResponseParams {
  // Aim to start influencing the missile no earlier than this window before
  // predicted impact (seconds). This avoids wasting countermeasures at long range.
  //
  // Set <= 0 to disable (release as soon as a useful seeker update exists).
  double decoyWindowSec{6.0};

  // Preferred lead time before a discrete seeker update tick (seconds).
  // The planner schedules the first burst at (tick - lead).
  double releaseLeadSec{0.05};

  // If time-to-impact is below this threshold, release immediately.
  double panicTtiSec{0.35};

  // If enabled and the missile uses discrete seeker updates, recommend
  // repeated bursts spaced by the update period.
  bool enableRepeatBursts{true};

  // Maximum number of bursts recommended before the caller should re-plan.
  int maxBursts{3};
};

struct CountermeasureResponsePlan {
  bool valid{false};

  CountermeasureType type{CountermeasureType::Flare};

  // Seconds from "now" to release the first burst.
  double firstReleaseDelaySec{0.0};

  // If > 0, caller may repeat bursts every this many seconds.
  double repeatEverySec{0.0};

  // Suggested number of bursts (including the first).
  int burstCount{0};

  // Debug/telemetry fields.
  double timeToImpactSec{0.0};
  double timeToSeekerActiveSec{0.0};
  double timeToFirstSeekerTickSec{0.0};
};

// Compute a deterministic countermeasure plan against a specific inbound missile.
//
// The plan is purely kinematic: it estimates time-to-impact from current relative
// velocity and aligns the first countermeasure burst to the next seeker update tick
// (or to seeker activation when the seeker is not yet active).
CountermeasureResponsePlan planCountermeasureResponse(const Missile& missile,
                                                      const math::Vec3d& targetPosKm,
                                                      const math::Vec3d& targetVelKmS,
                                                      const CountermeasureResponseParams& params = {});

// Compute a suggested evasion direction against a specific missile/target pair.
//
// seed:
//   Any stable seed (hash(universeSeed, npcId, timePhase, etc.)). Only used to break
//   ties in degenerate head-on cases.
MissileEvasionPlan planMissileEvasion(const Missile& missile,
                                     const math::Vec3d& targetPosKm,
                                     const math::Vec3d& targetVelKmS,
                                     core::u64 seed,
                                     const MissileEvasionParams& params = {});

// Find the nearest (minimum time-to-impact) inbound missile tracking the given
// target, according to a simple relative kinematic estimate.
MissileThreatSummary nearestInboundMissile(const Missile* missiles,
                                          std::size_t missileCount,
                                          CombatTargetKind targetKind,
                                          core::u64 targetId,
                                          const math::Vec3d& targetPosKm,
                                          const math::Vec3d& targetVelKmS,
                                          const MissileThreatParams& params = {});




// -----------------------------------------------------------------------------
// Jammer response planning (radar seekers)
// -----------------------------------------------------------------------------
//
// A lightweight deterministic heuristic for when a defender should emit a noise
// jammer against an inbound radar seeker.
//
// Why this exists:
//   - Noise jamming can reduce seeker track quality (harder guidance).
//   - However, some missiles support Home-On-Jam (HOJ). If HOJ would engage,
//     keeping the jammer on while attempting a doppler notch can be counter-productive,
//     because HOJ bypasses notch-based lock loss.
//   - The strategy here is: *jam when it matters*, but *shut down to deny HOJ*
//     when geometry is already near-notched or extremely close.
//
// The caller decides how to map `jammerPower` to an actual ship system.

enum class JammerResponseReason {
  None = 0,
  NotRadarThreat,
  NotInbound,
  TooFar,
  TooWeak,
  DenyHomeOnJamNotch,
  DenyHomeOnJamClose,
  DenyHomeOnJamTerrainMask,
};

struct JammerResponseParams {
  bool enable{true};

  // Jammer power level to use when enabled (same units as SphereTarget::jammerPower).
  double jammerPower{1.0};

  // Do not bother jamming if the estimated time-to-impact exceeds this window (seconds).
  // Set <= 0 to always consider jamming.
  double maxTtiToJamSec{60.0};

  // Require at least this received jamming level (0..1) before recommending the jammer.
  double minReceivedJamming01{0.05};

  // If the inbound radar missile supports HOJ and would receive enough jamming to engage it,
  // disable the jammer when the geometry is already near the doppler notch.
  bool denyHomeOnJam{true};

  // Disable when abs(range-rate) <= notch * multiple (only applies when notch > 0).
  double denyHomeOnJamNotchMultiple{1.2};

  // Also disable the jammer extremely close to impact (seconds) when HOJ would engage.
  double denyHomeOnJamTtiSec{2.0};
};

struct JammerResponsePlan {
  bool valid{false};

  bool jammerOn{false};
  double jammerPower{0.0};

  // Estimated received jamming level at the missile (0..1) if the jammer were on.
  double receivedJamming01{0.0};

  // True if the threat missile supports HOJ and would be able to use it with the suggested jammer power.
  bool homeOnJamCapable{false};

  JammerResponseReason reason{JammerResponseReason::None};

  // Debug/telemetry fields.
  double timeToImpactSec{0.0};
  double closingKmS{0.0};

  // Doppler range-rate used by notch/HOJ logic:
  //   rangeRate = dot(targetVel - missileVel, LOS_dir)  (km/s)
  double rangeRateKmS{0.0};
};

JammerResponsePlan planJammerResponse(const Missile& missile,
                                     const math::Vec3d& targetPosKm,
                                     const math::Vec3d& targetVelKmS,
                                     const JammerResponseParams& params = {});

// -----------------------------------------------------------------------------
// Integrated defense planning (deterministic)
// -----------------------------------------------------------------------------
//
// Convenience helper for gameplay/NPC code: combines threat detection with the
// existing per-threat planners (evasion, terrain masking, countermeasures) and
// returns a single coherent recommendation.
//
// This keeps "what should I do right now?" logic out of the main game loop,
// while remaining fully deterministic and headless.

enum class DefenseManeuverKind : core::u8 {
  None = 0,
  EvasionJink = 1,
  TerrainMask = 2,
};

struct IntegratedDefenseParams {
  MissileThreatParams threat{};
  MissileEvasionParams evasion{};
  MissileTerrainMaskParams terrain{};
  CountermeasureResponseParams countermeasures{};
  JammerResponseParams jammer{};

  // If enabled, terrain masking can be selected when the inbound missile depends
  // on line-of-sight (active seeker) and/or midcourse datalink LOS.
  bool enableTerrainMasking{true};

  // Minimum time-to-impact required before terrain masking will be selected
  // (even if a valid cover point exists).
  double minTtiForTerrainSec{1.0};

  // Urgency mapping:
  //   urgency01 = clamp01(1 - tti / urgencyTimeScaleSec)
  double urgencyTimeScaleSec{6.0};
};

struct IntegratedDefensePlan {
  bool valid{false};

  MissileThreatSummary threat{};

  DefenseManeuverKind maneuver{DefenseManeuverKind::None};

  // Unit vector (world space) for the recommended thrust direction.
  // When maneuver == None, this is {0,0,0}.
  math::Vec3d maneuverDirWorld{0, 0, 0};

  // Rough [0,1] urgency scalar derived from time-to-impact.
  double urgency01{0.0};

  // Sub-plans for introspection/telemetry.
  MissileEvasionPlan evasion{};
  MissileTerrainMaskPlan terrain{};
  CountermeasureResponsePlan countermeasures{};
  JammerResponsePlan jammer{};
};

// Plan defense against the nearest inbound missile (min time-to-impact).
//
// worldTargets/worldTargetCount should include asteroid SphereTargets (kind Asteroid)
// if terrain masking is desired.
IntegratedDefensePlan planIntegratedDefense(const Missile* missiles,
                                            std::size_t missileCount,
                                            CombatTargetKind targetKind,
                                            core::u64 targetId,
                                            const math::Vec3d& targetPosKm,
                                            const math::Vec3d& targetVelKmS,
                                            const SphereTarget* worldTargets,
                                            std::size_t worldTargetCount,
                                            core::u64 seed,
                                            const IntegratedDefenseParams& params = {});


// -----------------------------------------------------------------------------
// Point defense fire control (deterministic)
// -----------------------------------------------------------------------------
//
// A lightweight helper intended for NPCs / gameplay code:
//  - pick a ballistic lead direction to intercept an inbound missile
//  - respect simple turret constraints (slew rate + mount cone)
//  - optionally require asteroid-free line-of-sight to the intercept point
//
// The output is intentionally a *single-shot* recommendation. Callers can re-plan
// each frame or at a lower rate to create continuous tracking behavior.

struct PointDefenseParams {
  // Projectile muzzle speed relative to the firing platform (km/s).
  // Must be > 0.
  double muzzleSpeedKmS{0.0};

  // Maximum future time horizon to search for an intercept (seconds).
  // The effective horizon is also clamped by missile.ttlSimSec when > 0.
  double maxPlanSimSec{6.0};

  // Prediction step used for searching the intercept window (seconds).
  // If <= 0, a conservative default (~0.02s) is used.
  double simStepSec{0.02};

  // If true (default), use a small forward simulation via `stepMissiles()` to
  // predict missile motion. When false, fall back to simple inertial motion.
  bool useFullMissileSimulation{true};

  // Turret mount axis (unit) used for optional aim-cone gating.
  // When mountConeCos >= -0.5, a candidate aim direction must satisfy:
  //   dot(mountAxisWorld, aimDirWorld) >= mountConeCos
  math::Vec3d mountAxisWorld{0, 0, 1};
  double mountConeCos{-1.0};

  // Optional turret slew limit (rad/s).
  // When > 0, the turret can only rotate at this rate from its current direction.
  double turretSlewRateRadS{0.0};

  // Optional: require that the firing platform has a clear line-of-sight to the
  // intercept point at fire time (not occluded by asteroid SphereTargets).
  bool requireLineOfSight{false};
  double lineOfSightOcclusionPadKm{0.0};

  // Selection policy:
  //  - true: pick the earliest intercept time
  //  - false: pick the earliest fire time
  bool preferEarliestIntercept{true};
};

struct PointDefenseSolution {
  bool valid{false};

  // Seconds from "now" to fire.
  double fireDelaySec{0.0};

  // Seconds from "now" to intercept.
  double interceptTimeSec{0.0};

  // Unit vector (world space) for the firing direction.
  math::Vec3d aimDirWorld{0, 0, 0};

  // Predicted intercept point (world space, km).
  math::Vec3d interceptPointKm{0, 0, 0};

  // Debug/telemetry.
  double aimAngleRad{0.0};
  double requiredSlewTimeSec{0.0};
  double rangeAtInterceptKm{0.0};
};

// Plan a single projectile shot intended to intercept a specific missile.
//
// shooterPosKm/shooterVelKmS represent the firing platform kinematics.
// turretAimDirWorld is the current turret pointing direction (unit).
//
// worldTargets/worldTargetCount are only required when:
//  - params.useFullMissileSimulation is true (for guidance/LOS effects), and/or
//  - params.requireLineOfSight is true (for asteroid occlusion checks).
PointDefenseSolution planPointDefenseShot(const Missile& missile,
                                          const math::Vec3d& shooterPosKm,
                                          const math::Vec3d& shooterVelKmS,
                                          const math::Vec3d& turretAimDirWorld,
                                          const SphereTarget* worldTargets,
                                          std::size_t worldTargetCount,
                                          core::u64 seed,
                                          const PointDefenseParams& params = {});


// -----------------------------------------------------------------------------
// Multi-threat point defense engagement scheduling (deterministic)
// -----------------------------------------------------------------------------
//
// This helper extends the single-shot `planPointDefenseShot()` planner to handle
// multiple inbound missiles in a short horizon. It produces a conservative
// (turret-slew + cooldown aware) sequence of shots by repeatedly:
//   - selecting the nearest inbound threat (min time-to-impact)
//   - planning a feasible intercept
//   - advancing a local simulation forward to the next shot opportunity
//
// It is intentionally a deterministic heuristic, not an optimizer. Callers can
// re-plan each frame (or at a lower rate) to adapt to maneuvers and attrition.

struct PointDefenseEngagementParams {
  // Threat filtering / inbound detection.
  MissileThreatParams threat{};

  // Per-shot point-defense solver parameters.
  PointDefenseParams pointDefense{};

  // Minimum time between consecutive shots (seconds).
  // This is a simple reload / ROF constraint for scheduling.
  double shotCooldownSec{0.10};

  // Maximum number of shots to schedule.
  int maxShots{3};

  // Total planning horizon for the engagement schedule (seconds).
  // When <= 0, the horizon is derived from pointDefense.maxPlanSimSec.
  double maxPlanSimSec{0.0};

  // Require that planned intercept occurs at least this many seconds before the
  // threat's estimated time-to-impact.
  double interceptBeforeImpactMarginSec{0.05};

  // If true (default), clamp each per-shot planning horizon to the selected
  // threat's time-to-impact estimate.
  bool clampShotHorizonToThreatTti{true};
};

struct PointDefenseEngagementShot {
  bool valid{false};

  // Threat metadata at the time the shot was planned.
  MissileThreatSummary threat{};

  // Fire-control solution. Times are expressed as seconds-from-*initial-now*
  // (i.e., absolute within this plan).
  PointDefenseSolution solution{};
};

struct PointDefenseEngagementPlan {
  bool valid{false};

  // Ordered shot schedule.
  std::vector<PointDefenseEngagementShot> shots{};
};

// Plan a conservative schedule of point-defense shots against multiple inbound missiles.
//
// turretAimDirWorld is the current turret pointing direction (unit).
//
// worldTargets/worldTargetCount should match the target array used by the main sim
// when pointDefense.useFullMissileSimulation is enabled (for accurate guidance/LOS).
PointDefenseEngagementPlan planPointDefenseEngagement(const Missile* missiles,
                                                      std::size_t missileCount,
                                                      CombatTargetKind targetKind,
                                                      core::u64 targetId,
                                                      const math::Vec3d& targetPosKm,
                                                      const math::Vec3d& targetVelKmS,
                                                      const math::Vec3d& turretAimDirWorld,
                                                      const SphereTarget* worldTargets,
                                                      std::size_t worldTargetCount,
                                                      core::u64 seed,
                                                      const PointDefenseEngagementParams& params = {});

// -----------------------------------------------------------------------------
// Integrated defense + point defense convenience (deterministic)
// -----------------------------------------------------------------------------
//
// Convenience helper for gameplay/NPC code: runs `planIntegratedDefense()` and,
// if enabled, also computes a point-defense fire-control solution against the
// selected nearest inbound threat.
//
// This is intended to keep glue logic ("pick threat -> maneuver + decoys + shoot")
// out of higher level game code while remaining fully deterministic and headless.

struct IntegratedDefenseWithPointDefenseParams {
  IntegratedDefenseParams defense{};
  PointDefenseParams pointDefense{};

  // Enable/disable point defense planning (muzzleSpeedKmS must still be > 0).
  bool enablePointDefense{true};

  // Only attempt point defense when the inbound missile time-to-impact estimate
  // is within this window (seconds). Prevents firing at extreme range.
  double maxTtiForPointDefenseSec{6.0};

  // Require that the predicted intercept occurs at least this many seconds
  // before the time-to-impact estimate.
  double interceptBeforeImpactMarginSec{0.05};

  // If true (default), clamp pointDefense.maxPlanSimSec to the threat time-to-impact.
  bool clampPointDefenseHorizonToTti{true};
};

struct IntegratedDefenseWithPointDefensePlan {
  bool valid{false};

  IntegratedDefensePlan defense{};
  PointDefenseSolution pointDefense{};
};

// Plan integrated defense AND an optional point-defense shot against the nearest inbound missile.
//
// turretAimDirWorld is the current turret pointing direction (unit).
//
// worldTargets/worldTargetCount:
//   Passed through to both terrain masking and point-defense planning. For highest fidelity
//   with `pointDefense.useFullMissileSimulation`, provide the same target array you pass
//   to your main combat step (ships/player/decoys/asteroids).
IntegratedDefenseWithPointDefensePlan planIntegratedDefenseWithPointDefense(const Missile* missiles,
                                                                            std::size_t missileCount,
                                                                            CombatTargetKind targetKind,
                                                                            core::u64 targetId,
                                                                            const math::Vec3d& targetPosKm,
                                                                            const math::Vec3d& targetVelKmS,
                                                                            const math::Vec3d& turretAimDirWorld,
                                                                            const SphereTarget* worldTargets,
                                                                            std::size_t worldTargetCount,
                                                                            core::u64 seed,
                                                                            const IntegratedDefenseWithPointDefenseParams& params = {});


}  // namespace stellar::sim
