# StellarForge — Procedural Space Sim Prototype (C++20)

A small C++20 procedural space-sim framework + an SDL2/OpenGL **3D prototype** you can evolve toward a Pioneer / Elite-style loop.

This repo contains:
- Deterministic **procedural generation** (galaxy → systems → planets → stations)
- A streaming **Universe** with LRU caches
- A real-time **flight prototype** with ship physics, docking, markets, and combat-lite threats

## What’s playable right now (prototype “feel” pass)

### Stations that exist in-world (with orientation)
- Stations are placed on **orbits** in the current system.
- Stations have **visible geometry** (ring + central body + a highlighted docking frame).
- Stations slowly **spin**, so orientation matters when approaching.

### “Mail-slot” style docking (with clearance)
- Target a station: **T** cycles stations.
- Request docking clearance: **L** (must be in comms range).
- Fly through the **mail-slot** (stay aligned, stay under speed limit).
- Press **G** to dock / undock.

### Combat-lite loop (pirates + basic weapon)
- Pirates periodically spawn as contacts and will **pursue and fire** on the player.
- Player weapon: **laser (Left Mouse)**.
- If you die: you **respawn** and lose cargo + 10% credits (so pirates matter).
- Stations provide light **defensive fire** against pirates near the station.

### Market + basic services gating
- **Buying/Selling is disabled unless docked at that station.**
- While docked you can use a basic **Repair** service.

## Controls (quick)
- Translate: **W/S** (forward/back), **A/D** (strafe), **R/F** (up/down)
- Rotate: **Arrow keys** (pitch/yaw), **Q/E** (roll)
- Boost: **LShift**   Brake: **X**
- Dampers: **Z** on / **C** off
- Target station: **T** (cycle), clear target: **Y**
- Request clearance: **L**
- Dock/Undock: **G**
- Autopilot approach: **P**
- Fire laser: **Left Mouse**
- UI: **TAB** Galaxy, **F1** Flight, **F2** Market, **F3** Contacts
- Save/Load: **F5 / F9**
- Pause: **Space**

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
- Supercruise + nav assist (“7-second rule”)
- FSD/hyperspace loop (fuel + cooldown + arrival at station)
- Missions (courier/delivery/bounty) + reputation
- NPC traffic that actually moves station inventories/prices
