- **Round N:** Radar seeker upgrades: added `SphereTarget::radarSignature` (RCS-like weighting), plus optional distance-aware jammer reception (HOJ gating + radar track-quality suppression) with regression tests.

- **Round N:** Seeker gimbal slew now modulates decoy attraction + autonomous acquisition scoring (large off-axis pulls are suppressed unless the seeker can re-point quickly), with a regression test.

- **Round N-1:** Missile guidance now consumes `SphereTarget::accelKmS2` for acceleration-aware lead pursuit and better inertial target-memory prediction; APN feed-forward prefers sensor accel when available.

- **Round 2:** Windows build fixes (UI loop/save/menu). `NOMINMAX` + `WIN32_LEAN_AND_MEAN`, Time Trial ghost restore compile fix, and UI updates for `SystemConditionsSnapshot` + `CStrView` conversions.
- **Round 3:** Orbit Analyzer maneuver planner now wires cleanly to the in-game RTN maneuver node (dv along/normal/radial + ref-body choice) and compiles (bindings + function signature aligned).
- `apps/stellar_game/OrbitAnalyzerWindow.cpp`: fixed gravity body kind enum usage and corrected impact/closest-approach readouts (impact altitude computed from distance-to-body).
- `apps/stellar_game/main.cpp`: fixed Orbit Analyzer call wiring, GalNet digest -> Integration Hub event push, and removed undefined supercruise/hyperspace guards.
- `apps/stellar_game/GalNetInboxWindow.cpp`: updated to current `CommsLog` API (items/markRead/markPinned) and resolves system names via `Universe`.

- **Round 4:** Heat-seeker IR aspect weighting (tail-chase hotter, head-on cooler) added to missile target scoring (incl. decoy competition), with tunable parameters and regression tests.
- `include/stellar/sim/Combat.h`: added per-missile heat-aspect tuning fields (front/rear factors + speed blend).
- `src/sim/Combat.cpp`: applies aspect-weighted heat signature scoring for locked/memory/reacquire logic; sets default aspect tuning for heat-seeker missiles.
- `tests/test_combat.cpp`: added heat-aspect tests validating tail-on lock is harder to decoy than head-on.


Verify:
- Build the `game-release` preset.
- In-game: open Orbit Analyzer, click a planner action (e.g. Circularize) and confirm maneuver node + trajectory preview update; check Forecast impact/approach tables.
- In-game: open GalNet Inbox, select bulletins, mark read/unread/pin, and request an update; enable watchlist digest and confirm an event appears in Integration Hub.
## Round 5
- `include/stellar/sim/Combat.h`: added `Missile::seekerUpdatePeriodSimSec` to model discrete seeker scan/update ticks (0 = legacy continuous).
- `src/sim/Combat.cpp`: active seeker can now be rate-limited; between scan ticks missiles guide on inertial target memory and defer decoy + auto-acquire decisions.
- `tests/test_combat.cpp`: added deterministic test covering delayed decoy commitment for rate-limited seekers vs immediate commitment for continuous seekers.

## Round 6
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: added deterministic countermeasure response planner that recommends flare/chaff timing aligned to seeker scan ticks.
- `tests/test_countermeasure_response.cpp`: regression tests for scan-aligned scheduling, seeker activation delay, and invalid edge case.

## Round 7
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: added deterministic terrain-masking planner (`planMissileTerrainMasking`) that recommends a cover direction behind nearby asteroids relative to an inbound missile.
- `tests/test_terrain_masking.cpp`: regression tests for cover-point selection + the "already occluded" early-out.

## Round 8
- `include/stellar/sim/Combat.h` / `src/sim/Combat.cpp`: added projectile-driven missile interception via `stepMissilesWithProjectileInterception` + detonation metadata (`MissileDetonation::intercepted`, interceptor ids).
- `tests/test_projectile_missile_interception.cpp`: regression test for deterministic mid-course shootdown that preempts a would-be ship impact.

## Round 9
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: added `planIntegratedDefense()` to fuse threat selection + evasion + terrain masking (LOS-aware) + countermeasure timing into one deterministic NPC-friendly output.
- `tests/test_integrated_defense.cpp`: regression tests covering TerrainMask selection only when LOS-breaking is relevant; otherwise falls back to EvasionJink.

## Round 10
- `include/stellar/sim/Combat.h`: added `Missile::datalinkUpdatePeriodSimSec` + `Missile::datalinkLatencySimSec` for discrete pre-activation midcourse updates.
- `src/sim/Combat.cpp`: midcourse datalink can now be rate-limited + latency-shifted; between ticks guidance relies on inertial target memory (with an implied memory bridge when period>0).
- `tests/test_datalink_midcourse.cpp`: regression test ensuring target velocity changes are only ingested on the next (latency-shifted) datalink tick, even when `targetMemorySimSec` is left at 0.

## Round 11
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: added deterministic point-defense fire-control planner `planPointDefenseShot()` (projectile lead solution) with turret slew + mount-cone constraints and optional asteroid line-of-sight gating.
- `tests/test_point_defense_planner.cpp`: regression tests covering immediate-fire intercept, slew-limited delayed shot, and invalid case when the turret cannot rotate in time.


## Round 12
- `include/stellar/sim/Combat.h` / `src/sim/Combat.cpp`: added deterministic terminal weave (endgame jink) for missiles (`terminalWeave*` params) and enabled a small weave on RadarMissile defaults to complicate point defense.
- `tests/test_terminal_weave.cpp`: regression test validating deterministic weave-induced lateral oscillation.


## Round 13
- `include/stellar/sim/Combat.h` / `src/sim/Combat.cpp`: added optional deterministic track-error (sensor noise) that scales with seeker `trackQuality` (`trackErrorMaxRad`, `trackErrorFrequencyHz`) and is applied to non-decoy guidance while the active seeker runs.
- `tests/test_track_error.cpp`: regression test validating the phase-stable lateral bias when track quality is poor, no bias at `trackQuality==1`, and determinism.


## Round 14
- `include/stellar/sim/Combat.h` / `src/sim/Combat.cpp`: added optional seeker received-signal falloff (`seekerSignalHalfRangeKm`, `seekerSignalMinDistKm`) that modulates `trackQuality` (and therefore decoy susceptibility + track error) with range.
- `tests/test_seeker_signal_falloff.cpp`: regression test validating long-range signal falloff triggers decoy commitment while short-range keeps lock, plus determinism.


## Round 15
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: added `planIntegratedDefenseWithPointDefense()` meta-planner that combines the existing integrated defense output with an optional point-defense shot recommendation (with a time-to-impact clamp).
- `tests/test_integrated_defense_point_defense.cpp`: regression tests covering a valid PD recommendation when feasible, plus invalid cases for turret-slew limits and long-range TTI gating.


## Round 16
- `include/stellar/sim/Combat.h` / `src/sim/Combat.cpp`: added optional countermeasure "clutter" track-quality suppression (`Missile::decoyClutterTrackSuppressionGain`) so strong in-FOV decoys can degrade seeker measurement quality even before a full lock swap.
- `src/sim/Combat.cpp`: enabled conservative clutter suppression defaults for both RadarMissile (chaff) and Heat missiles (flares).
- `tests/test_decoy_clutter.cpp`: regression tests validating deterministic trackQuality suppression when gain>0 and legacy behavior when gain==0.


## Round 17
- `include/stellar/sim/Combat.h` / `src/sim/Combat.cpp`: added optional seeker-lag guidance bias (`Missile::seekerLagGuidanceGain`) so slow gimbal slew can produce real guidance error (not only lower track quality) during close-range high-LOS-rate fights.
- `tests/test_seeker_lag_guidance.cpp`: regression test validating that enabling lag bias changes the commanded steering direction deterministically.


## Round 18
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: added deterministic jammer-response planner (`planJammerResponse`) for radar threats, including Home-On-Jam (HOJ) denial heuristics (near-notch, close-range, and terrain-mask override).
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: integrated jammer recommendations into `planIntegratedDefense()` as an additional sub-plan for NPC decision making.
- `tests/test_jammer_response.cpp`: regression tests for non-HOJ jamming, HOJ denial (notch + close), far-range gating, and the terrain-mask override.


## Round 19
- `include/stellar/sim/MissileDefense.h` / `src/sim/MissileDefense.cpp`: added deterministic multi-threat point-defense engagement scheduler `planPointDefenseEngagement()` (shot sequence with cooldown + per-shot intercept feasibility).
- `tests/test_point_defense_engagement.cpp`: regression tests for two-shot scheduling and a cooldown-limited horizon case.


## Round 20
- `include/stellar/sim/BearingTrack.h`: added `BearingTrackSolveDiagnostics` + `bearingTrackSolveDiagnostics()` helper exposing the same weight/determinant solvability gates used by `updateBearingTrack()`.
- `include/stellar/sim/CrossFixPlanner.h`: added a deterministic cross-fix waypoint recommender for bearing-only tracks (preferred + alternate side) with baseline selection from sigma or range guess.
- `tests/test_bearing_track_diagnostics.cpp`: regression tests for solvability gate + progress behavior.
- `tests/test_cross_fix_planner.cpp`: regression tests for baseline selection, perpendicular move directions, expected bearing change, and deterministic side selection.

## Round 21
- Added BearingTrack uncertainty ellipsoid helper (principal axes + sigma lengths) for visualizing bearing-only fix conditioning.
- Upgraded CrossFixPlanner with optional predicted information-gain mode: can compute predicted solvability/conditioning per candidate, choose side, and optionally optimize baseline.
- Added tests covering uncertainty ellipsoid behavior and cross-fix prediction/baseline optimization.

## Round 22
- `include/stellar/sim/NavRoute.h` / `src/sim/NavRoute.cpp`: exposed hard routing constraints for A* planners (ban nodes + undirected edges) to support avoid-lists and detour analysis.
- `include/stellar/sim/NavRouteResilience.h`: added replacement-path style edge-detour resilience report ("best route if a jump fails").
- `tests/test_nav_route_resilience.cpp`: regression test covering deterministic detour selection and critical-edge detection.

## Round 23
- `include/stellar/sim/NavRoute.h` / `src/sim/NavRoute.cpp`: added MMR-style diversified K-route selection (`selectDiversifiedKRoutesMMR`) plus a convenience wrapper (`plotKRoutesAStarCostDiversified`).
- `tests/test_nav_route_diversified_k.cpp`: regression tests validating the cost-vs-diversity tradeoff and deterministic behavior.

## Round 24
- `include/stellar/sim/NavRouteConnectivity.h`: added deterministic edge-disjoint routing analysis (unit-capacity max-flow) with a minimum-cut report. Supports the same banned node/edge constraints as constrained A* routing.
- `tests/test_nav_route_edge_disjoint.cpp`: regression tests for 3-way edge-disjoint corridors, banned-edge behavior, and deterministic min-cut extraction.
