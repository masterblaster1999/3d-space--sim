#pragma once

#include "stellar/core/Types.h"

namespace stellar::proc {

// Deterministic 2D value noise and fBm (fractal brownian motion) helpers.
// Useful for later: terrain, gas giant bands, nebula density fields, etc.

// Lattice noise at integer coordinates, returns [0,1).
double valueNoise2D(int x, int y, stellar::core::u64 seed);

// Continuous value noise at real coordinates, returns ~[0,1).
double valueNoise2D(double x, double y, stellar::core::u64 seed);

// Fractal brownian motion: sum of multiple octaves of value noise.
double fbm2D(double x, double y, stellar::core::u64 seed,
            int octaves = 5, double lacunarity = 2.0, double gain = 0.5);

} // namespace stellar::proc
