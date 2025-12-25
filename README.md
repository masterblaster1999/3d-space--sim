# StellarForge — Procedural Space Sim Prototype (C++20)

A small C++20 procedural space-sim framework + an SDL2/OpenGL **3D prototype** you can evolve toward a Pioneer / Elite-style loop.

This repo contains:
- Deterministic **procedural generation** (galaxy → systems → planets → stations)
- A streaming **Universe** with LRU caches
- A real-time **flight prototype** with ship physics, docking, markets, and combat-lite threats

## What’s playable right now (prototype “feel” pass)

### Stations + docking
- Stations orbit planets and broadcast a **comms range**.
- You can request docking clearance (**L**) and dock (**G**) when in range.
- While docked you can use basic services:
  - **Repair** (cost scales with damage)
  - **Refuel** (buys Fuel from the station and fills your tank)
  - **Shipyard upgrades** (at Shipyards): cargo rack + fuel tank upgrades

### Economy / market prototype
- Each station has its own inventory and prices (simple supply/demand model).
- Buy & sell commodities.
- Cargo is now **mass-limited (kg)** instead of infinite.

### In-system flight
- 6DOF-ish flight with dampers / boost / brake.
- Autopilot can fly to the current target.

### Galaxy streaming
- A simple galaxy map that streams nearby systems from a procedural seed.

### FSD jump (early hyperspace loop)
- Select a destination system in the **Galaxy** window and press **J** (or click the button).
- Requires being undocked, not mass-locked near a station, and having enough fuel.
- Has a short **charge** and **cooldown** (balance tuning is placeholder).

### Missions (early)
- Dock at a station and open the **Missions** window (**F4**).
- Accept **Courier** and **Delivery** jobs.
- Delivery jobs may require you to buy cargo, or may provide cargo up-front (requires cargo space).
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
- Fire laser: **Left Mouse**
- FSD jump: **J** (to selected system in Galaxy)
- UI: **TAB** Galaxy, **F1** Flight, **F2** Market, **F3** Contacts, **F4** Missions, **F6** Scanner
- Save/Load: **F5 / F9**
- Pause: **Space**
- Time scale: **[ / ]** cycle simulation speed

## Build

This is a CMake project.

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

### GitHub Actions / Autotools templates
Some CI templates run `./configure && make`. This repo now includes a **compat** `configure` script + `Makefile` that delegate to CMake.

## Project structure
- `include/stellar/...` core types, math, sim, econ, render
- `src/...` implementations
- `apps/stellar_game` SDL2/OpenGL prototype
- `tests` unit tests (headless)

## Next “feel” upgrades worth doing
- Supercruise drop/approach tuning (more natural exit + approach corridor guidance)
- Mission variety (bounty / passenger / multi-hop) + reputation
- Mission board persistence + better failure/penalty handling
- NPC traffic that actually moves station inventories/prices
