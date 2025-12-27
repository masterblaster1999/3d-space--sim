# StellarForge — Procedural Space Sim Prototype (C++20)

A small C++20 procedural space-sim framework + an SDL2/OpenGL **3D prototype** you can evolve toward a Pioneer / Elite-style loop.

This repo contains:
- Deterministic **procedural generation** (galaxy → systems → planets → stations)
- A streaming **Universe** with LRU caches
- A real-time **flight prototype** with ship physics, docking, markets, and combat-lite threats
- "Ambient" **NPC trade traffic** that nudges station inventories (and therefore prices) over time
- **Squad encounters**: pirate packs + police wingmen, with basic morale/pursuit escalation

## What’s playable right now (prototype “feel” pass)

### Stations + docking
- Stations orbit planets and broadcast a **comms range**.
- You can request docking clearance (**L**) and dock (**G**) when in range.
- While docked you can use basic services:
  - **Repair** (cost scales with damage)
  - **Refuel** (buys Fuel from the station and fills your tank)
  - **Shipyard upgrades** (at Shipyards): cargo rack + fuel tank upgrades
  - Passenger cabins (at Shipyards): increases **passenger seat** capacity for Passenger missions

### Economy / market prototype
- Each station has its own inventory and prices (simple supply/demand model).
- Buy & sell commodities.
- Cargo is now **mass-limited (kg)** instead of infinite.
- Background NPC trade traffic moves commodities between stations so markets drift even while you fly.

### In-system flight
- 6DOF-ish flight with dampers / boost / brake.
- Autopilot can fly to the current target.
- Planet/star collision + low-altitude heating (simple danger close to gravity wells).
- Projectile lead indicator for kinetic weapons (HUD).

### Galaxy streaming
- A simple galaxy map that streams nearby systems from a procedural seed.
- System IDs encode their sector so `Universe::getSystem(id)` works without extra context (useful for saves + nav routes).

### FSD jump (early hyperspace loop)
- Select a destination system in the **Galaxy** window and press **J** (or click the button).
- Requires being undocked, not mass-locked near a station, and having enough fuel.
- Has a short **charge** and **cooldown** (balance tuning is placeholder).

### Missions (early)
- Dock at a station and open the **Missions** window (**F4**).
- Accept **Courier**, **Delivery**, **Multi-hop Delivery**, **Smuggle**, and **Passenger** jobs.
- Some deliveries require you to source the cargo; others provide it up-front (requires cargo space).
- **Smuggle** jobs carry contraband (illegal cargo) and may attract police attention.
- Passenger jobs require enough **Passenger Seats** (buy cabin upgrades at Shipyards).
- Accept early **Bounty Scan** / **Bounty Kill** jobs (scan targets with the scanner, or destroy marked criminals).
- Complete by docking at the destination station before the deadline.
## Controls (quick)
- Translate: **W/S** (forward/back), **A/D** (strafe), **R/F** (up/down)
- Rotate: **Arrow keys** (pitch/yaw), **Q/E** (roll)
- Mouse steering: **M** toggles Mouse Steer mode (captures cursor); move mouse (pitch/yaw)
- Boost: **LShift**   Brake: **X**
- Dampers: **Z** on / **C** off
- Targets: **T** station, **B** planet, **N** contact, **U** star | clear target: **Y**
- Request clearance: **L**
- Dock/Undock: **G**
- Autopilot approach: **P**
- Fire primary/secondary: **Left Mouse / Right Mouse**
- Supercruise: **H** engage/drop (assist auto-drops in SAFE window)
- Interdiction: align to escape vector (HUD) | **H** submit
- FSD jump: **J** (to selected system in Galaxy)
- UI: **TAB** Galaxy, **F1** Flight, **F2** Market, **F3** Contacts, **F4** Missions, **F6** Scanner
- Save/Load: **F5 / F9**
- Pause: **Space**
- Time scale: **[ / ]** cycle simulation speed

## Build

This is a CMake project.

### CMake presets (recommended)
If you have CMake 3.20+, you can use the included `CMakePresets.json`:

```bash
cmake --preset headless-release
cmake --build --preset headless-release
ctest --preset headless-release
```

### Standard build (game + tests)
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/apps/stellar_game/stellar_game
```

### Headless build (core + tests only)
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTELLAR_ENABLE_RENDER=OFF -DSTELLAR_ENABLE_IMGUI=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
```

### Sandbox tools (headless)
`stellar_sandbox` can print nearby systems and (optionally) scan for profitable trade opportunities without running the SDL/OpenGL prototype.

```bash
# Print nearby systems, then scan for the best trades from the chosen origin station.
./build/apps/stellar_sandbox/stellar_sandbox --seed 1337 --radius 120 --trade --cargoKg 240 --fromSys 0 --fromStation 0

# Plan a jump route (A* hops) inside the queried node set.
# Tip: increase --radius/--limit if a route cannot be found.
./build/apps/stellar_sandbox/stellar_sandbox --seed 1337 --radius 200 --route --fromSys 0 --toSys 12 --jr 18

# Emit JSON for scripting / balancing spreadsheets.
./build/apps/stellar_sandbox/stellar_sandbox --seed 1337 --radius 120 --trade --cargoKg 240 --json --out trade_scan.json

# Print the deterministic mission board at a station (and optionally simulate acceptance/completion).
./build/apps/stellar_sandbox/stellar_sandbox --seed 1337 --radius 120 --missions --rep 25 --acceptOffer 0 --autoComplete

# Load an existing save, then write the resulting state back out (tooling).
./build/apps/stellar_sandbox/stellar_sandbox --seed 1337 --radius 120 --missions --load mysave.sav --acceptOffer 1 --save mysave_updated.sav
```

### GitHub Actions / Autotools templates
Some CI templates run `./configure && make`. This repo now includes a **compat** `configure` script + `Makefile` that delegate to CMake.

## Project structure
- `include/stellar/...` core types, math, sim, econ, render
- `src/...` implementations
- `apps/stellar_game` SDL2/OpenGL prototype
- `tests` unit tests (headless)

## Next “feel” upgrades worth doing
- Supercruise drop/approach tuning (more natural exit + approach corridor guidance)
- Passenger cabins + passenger missions
- A real "traffic layer": convoys you can interdict/protect (hook ambient trade into gameplay)
- More encounter variety: convoy escorts, distress rescues, salvage ambushes
- Deeper crime loop: fines vs bounties, bribes, faction-specific contraband lists per station