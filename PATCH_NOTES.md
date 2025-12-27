
## 2025-12-27 (Patch) - Alert decay + Lay Low station option
- **Bugfix/behavior:** local security `policeHeat` now decays over time (exponential decay), so "Alert" cools down naturally as time passes.
- **New dock service:** **Lay low (3h)** in **Security Services**: pay a fee to advance time by 3 hours and reduce local alert/response.

# Patch Notes (Dec 27, 2025 - Police Alert Fix + Port Security Services)

- **Bugfix:** contraband enforcement, police bribes, and black-market contraband sales now raise **local security alert** (`policeHeat`) instead of ship thermal heat.
- **UI clarity:** status labels now distinguish **Alert** (local security) from **Ship heat** (thermal/overheat).
- **New station option:** **Bribe port authority (reduce alert)** lowers local alert for a fee (does not clear bounties).
- **Black market consequence:** selling illegal goods applies a small reputation penalty and raises alert over time.

---

# Patch Notes (Dec 27, 2025 - New Commodities + Save Compatibility)

## Economy: 2 new commodities (new)
- Added two higher-value commodities:
  - **Weapons** (`ARMS`)
  - **Stimulants** (`STIM`)
- Updated station economy models so the new goods are actually produced/consumed:
  - Agricultural: produces small amounts of Stimulants
  - Industrial: produces Weapons
  - TradeHub/Research/Shipyard: increased demand for these goods

## Smuggling / contraband: more variety (upgrade)
- Contraband rules now include the new goods.
  - Factions may ban **Stimulants** as a vice.
  - Factions may ban **Weapons** as controlled tech.

## SaveGame: commodity-count forward compatibility (important)
- SaveGame version bumped to **12**.
- **Older saves remain loadable** even when commodity counts change:
  - `cargo` lines and station override `inventory` lines now parse the remainder of their line and
    default missing commodity slots to **0**.

---

# Patch Notes (Dec 27, 2025 - Bribes & Contraband Scans)

## New: Bribe & comply window for contraband scans (police)
- If a **police** cargo scan completes and detects contraband, they may offer a short "bribe or comply" window.
  - **C**: Pay bribe to keep your cargo.
  - **I**: Comply: cargo is confiscated and you are fined (unpaid fines become bounty).
  - **Run**: Leave scan range to keep cargo but receive a bounty and reputation hit.
- The window is visible in **Ship / Status** with a countdown timer.
- Police ships will generally **hold fire** during the bribe/compliance window (unless already hostile).

---

## Shipyard: Smuggling compartments (new)
- Added **Smuggling compartments Mk1–Mk3** as a Shipyard upgrade.
  - Reduces the *chance* of cargo scans when you are carrying contraband.
  - Makes scans take slightly longer (more time to break range if you try to run).
- Persisted via SaveGame (`smuggleHoldMk`). Save version bumped to **11** (older saves default to Mk0).

## Smuggling: meaningful failure case (new)
- If a cargo scan **confiscates contraband** that is tied to an active **Smuggle** mission, that mission now **fails immediately**.
  - Applies the usual mission-failure reputation hit with the issuing faction.

## Missions UI QoL
- Active cargo missions now show a quick **Cargo: have / need** line (Delivery / Multi-hop / Smuggle).

---

# Patch Notes (Dec 27, 2025 - Hotfix)

## Smuggle missions: completion fix (SDL prototype)
- Fixed the SDL prototype's dock-completion logic to include **Smuggle** missions.
  - Completing a Smuggle mission now consumes cargo and pays out as expected.
  - Smuggle deliveries still **do not** feed the destination's *legal* market inventory.

## Market QoL: reserve mission cargo
- The Market sell flow now reserves cargo required by active **Delivery**, **Multi-hop**, and **Smuggle** missions.
  - The Cargo column shows reserved units.
  - If all units are reserved, the Sell button is disabled with a tooltip explaining why.

---

# Patch Notes (Dec 27, 2025)

## Missions: Smuggling (new)
- Added **Smuggle** missions (new `MissionType::Smuggle`).
  - Smuggling jobs grant the contraband cargo from the contact on acceptance (no legal market inventory/credit checks).
  - Completion consumes the cargo and pays out normally, but **does not** feed the destination's legal market inventory.

## Simulation: shared contraband rules (refactor)
- Introduced `sim::Contraband` (`stellar/sim/Contraband.h`) to centralize deterministic contraband rules:
  - `illegalCommodityMask`, `isIllegalCommodity`, `illegalCommodityListString`, `hasIllegalCargo`
- The SDL prototype now uses the shared contraband helpers so illegality rules stay consistent across UI + headless modules.

## Tooling + tests
- `stellar_sandbox` now recognizes the **Smuggle** mission type name.
- Added a unit test to cover Smuggle mission acceptance + completion invariants.

---

# Patch Notes (Dec 26, 2025)

## Tooling: stellar_sandbox JSON + mission board (new)
- `stellar_sandbox` gained a small CLI parser and JSON emitter for scripting/balancing:
  - `--json` (stdout) and `--out <path>` (file) outputs nearby systems + optional trade ideas
- New headless mission workflow helpers:
  - `--missions` to print deterministic offers
  - `--acceptOffer <idx>` to simulate acceptance with cargo/credits/seats constraints
  - `--autoComplete` (+ `--advanceDays`) to test completion logic without the renderer
  - `--load/--save` to round-trip SaveGame files from tooling

## Navigation: jump route planner (new)
- Moved the game's A* hop-count jump planner into core: `sim::plotRouteAStarHops` (`stellar/sim/NavRoute.h`).
- `stellar_sandbox` gained `--route` mode to plot jump paths headlessly (supports `--json` output).
- Added unit test `test_nav` to cover route finding, unreachable cases, and determinism.

## Core utilities (new)
- Added `stellar/core/Args.h` (tiny arg parser) and `stellar/core/JsonWriter.h` (minimal JSON writer)


## Encounters: squads + pursuit escalation (new)
- Pirates can now spawn in **packs** (leader + wingmen) instead of always solo.
- Police can spawn with **wingmen**, and their response scales with your **local bounty** and a soft pursuit "heat".
- Pirates may **break off and flee** when badly damaged (adds a bit of morale/variety to fights).
- Contacts UI now shows squad role tags: `[LEAD]` / `[WING]`, plus local `Heat`.

## Passenger missions + cabins (new)
- Added a new mission type: **Passenger** (party size stored in `Mission::units`).
- Added `SaveGame::passengerSeats` with Shipyard upgrades to increase capacity.
- Added seat-capacity validation on mission acceptance + completion support.
- SaveGame version bumped to **10** (backwards-compatible: old saves default to 2 seats).

## Mission system scaffolding (new)
- Added a reusable, **headless** mission module: `sim::MissionLogic`.
  - Deterministic mission-board generation persisted in `SaveGame`.
  - Shared accept logic (cargo checks, inventory/credit handling, station fees).
  - Shared completion/deadline helpers for Courier/Delivery/Multi-hop.
  - Game UI now delegates mission-offer generation + acceptance to the shared module.
- Added lightweight helpers:
  - `econ::cargoMassKg(...)` (shared cargo mass calc)
  - `sim::Reputation` helpers (rep + fee shaping)

## Tooling / CI
- Added `CMakePresets.json` for one-command headless builds.
- Updated GitHub Actions to build via CMake+Ninja and run `ctest`.

## Tests
- Added `test_missions` (mission board determinism + accept/complete path).

This patch adds a low-cost **ambient NPC trade traffic** layer that moves commodities between stations, nudging inventories and prices over time (even while you fly). It also persists these "traffic day stamps" in the save file so markets remain deterministic across reloads.

## Build fixes
- Added missing asteroid mining fields used by gameplay code:
  - `AsteroidNode::yield`
  - `AsteroidNode::chunkAccumulator`
- Added missing signal-site one-shot flag:
  - `SignalSource::fieldSpawned`
- Added missing NPC combat stat fields referenced by encounter logic:
  - `Contact::shieldMax`, `Contact::shieldRegenPerSec`, `Contact::hullMax`
- Fixed commodity loop using the correct enum/constant (`econ::kCommodityCount` instead of `CommodityId::COUNT`).
- Updated `spawnCargoPod` to accept a 5th parameter for scatter velocity (fixes "term does not evaluate to a function taking 5 arguments").
- Fixed string concatenation with `const char*` that caused pointer arithmetic ("Scooped" toast).
- Fixed invalid `.c_str()` calls on `CommodityDef::name` (it is `const char*`).

## Small gameplay/quality tweaks
- Ensured spawned pirates/traders/police initialize `shieldMax/hullMax` consistently.
- Added a simple NPC shield regen loop (uses `Contact::shieldRegenPerSec`).

## Economy + tooling upgrades
- Route planner now filters out *infeasible* trades (no source inventory or no destination capacity).
- Added `bestRoutesForCargo(...)` which considers cargo mass limits and market fees, and returns net trip profit.
- `stellar_sandbox` gained a headless `--trade` scan mode so you can quickly find profitable routes without running the full SDL/OpenGL prototype.

## Tests + CI
- Restored a procedural determinism test (`test_proc`) against the current streaming `Universe`.
- Added unit tests for:
  - route planner feasibility and cargo/fee math (`test_route_planner`)
  - SaveGame round-trip serialization (`test_savegame`)
- GitHub Actions now runs a headless `make distcheck` build (render/ImGui disabled) to keep CI fast and reliable.

## Ambient NPC trade traffic (new)
- Added `sim::simulateNpcTradeTraffic(...)` (core library) to model background station-to-station trade flows.
- Traffic simulation is:
  - deterministic per (universe seed, system id, day stamp)
  - bounded in CPU cost (backfills a limited number of days on large time jumps)
- Added save/load support for traffic stamps (`SaveGame` version bumped to **10**).
- Updated `test_savegame` to cover the new serialization.

## Streaming fix: SystemId encodes sector (new)
- Fixed `GalaxyGenerator::makeSystemId(...)` to **pack sector coords + local index** into the 64-bit id.
  - This restores the intended behavior where `Universe::getSystem(id)` can decode the sector and fetch the correct stub **without requiring a hint**.
- Added a regression check to `test_streaming` to ensure `getSystem(id)` returns the correct stub position/fields on a fresh Universe.
- Improved `Universe::getSystem(...)` fallback stub generation to create a **plausible galaxy-disc position** (instead of collapsing to {0,0,0}).