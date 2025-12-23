#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Orbit.h"

#include <cstdint>
#include <string>
#include <vector>

namespace stellar::sim {

enum class StarClass {
  O, B, A, F, G, K, M
};

enum class PlanetType {
  Rocky,
  Ice,
  GasGiant,
  IceGiant,
  AsteroidBelt
};

struct Star {
  StarClass starClass = StarClass::G;
  double massSolar = 1.0;        // solar masses
  double radiusSolar = 1.0;      // solar radii
  double luminositySolar = 1.0;  // solar luminosities
  double temperatureK = 5778.0;  // effective temperature

  std::string name;
};

struct Planet {
  PlanetType type = PlanetType::Rocky;

  double massEarth = 1.0;     // Earth masses
  double radiusKm = 6371.0;   // km
  double equilibriumTempK = 288.0;

  Orbit orbit;

  std::string name;
};

struct StarSystem {
  std::uint64_t id = 0;
  stellar::math::Vec3d positionLy{}; // position in galaxy (light-years)

  Star primary;
  std::vector<Planet> planets;
};

inline std::string toString(StarClass c) {
  switch (c) {
    case StarClass::O: return "O";
    case StarClass::B: return "B";
    case StarClass::A: return "A";
    case StarClass::F: return "F";
    case StarClass::G: return "G";
    case StarClass::K: return "K";
    case StarClass::M: return "M";
  }
  return "?";
}

inline std::string toString(PlanetType t) {
  switch (t) {
    case PlanetType::Rocky:        return "Rocky";
    case PlanetType::Ice:          return "Ice";
    case PlanetType::GasGiant:     return "Gas Giant";
    case PlanetType::IceGiant:     return "Ice Giant";
    case PlanetType::AsteroidBelt: return "Asteroid Belt";
  }
  return "?";
}

} // namespace stellar::sim
