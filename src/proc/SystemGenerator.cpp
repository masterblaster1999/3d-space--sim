#include "stellar/proc/SystemGenerator.h"

#include "stellar/core/Random.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace stellar::proc {
namespace {

struct StarClassSpec {
  stellar::sim::StarClass cls;
  double weight;    // relative
  double massMin;
  double massMax;
};

constexpr StarClassSpec kStarClasses[] = {
  // Very rough main-sequence distribution (massive stars are rare).
  { stellar::sim::StarClass::O, 0.00003, 16.0, 50.0 },
  { stellar::sim::StarClass::B, 0.00130,  2.1, 16.0 },
  { stellar::sim::StarClass::A, 0.00600,  1.4,  2.1 },
  { stellar::sim::StarClass::F, 0.03000,  1.04, 1.4 },
  { stellar::sim::StarClass::G, 0.07600,  0.80, 1.04 },
  { stellar::sim::StarClass::K, 0.12100,  0.45, 0.80 },
  { stellar::sim::StarClass::M, 0.76567,  0.08, 0.45 },
};

stellar::sim::StarClass pickStarClass(stellar::core::SplitMix64& rng) {
  double total = 0.0;
  for (const auto& s : kStarClasses) total += s.weight;

  const double r = rng.uniform(0.0, total);
  double accum = 0.0;
  for (const auto& s : kStarClasses) {
    accum += s.weight;
    if (r <= accum) return s.cls;
  }
  return stellar::sim::StarClass::G;
}

StarClassSpec specFor(stellar::sim::StarClass cls) {
  for (const auto& s : kStarClasses) {
    if (s.cls == cls) return s;
  }
  return { stellar::sim::StarClass::G, 1.0, 0.8, 1.04 };
}

// Rough main-sequence approximations.
// (These are intentionally simple — good enough for a game starting point.)
double approxLuminositySolar(double massSolar) {
  const double m = std::max(0.08, massSolar);
  if (m < 0.43) {
    return 0.23 * std::pow(m, 2.3);
  }
  if (m < 2.0) {
    return std::pow(m, 4.0);
  }
  if (m < 20.0) {
    return 1.5 * std::pow(m, 3.5);
  }
  return 32000.0 * m; // very rough
}

double approxRadiusSolar(double massSolar) {
  const double m = std::max(0.08, massSolar);
  return std::pow(m, 0.8);
}

double approxTemperatureK(double luminositySolar, double radiusSolar) {
  // From L ∝ R^2 T^4 -> T ∝ (L / R^2)^(1/4)
  const double ratio = std::max(1e-9, luminositySolar / (radiusSolar * radiusSolar));
  return 5778.0 * std::pow(ratio, 0.25);
}

std::string romanNumeral(int n) {
  struct Entry { int value; const char* numeral; };
  constexpr Entry table[] = {
    {1000, "M"}, {900, "CM"}, {500, "D"}, {400, "CD"},
    {100, "C"}, {90, "XC"}, {50, "L"}, {40, "XL"},
    {10, "X"}, {9, "IX"}, {5, "V"}, {4, "IV"}, {1, "I"}
  };

  std::string out;
  for (const auto& e : table) {
    while (n >= e.value) {
      out += e.numeral;
      n -= e.value;
    }
  }
  return out;
}

stellar::sim::PlanetType pickPlanetType(double aAU, double snowLineAU, stellar::core::SplitMix64& rng) {
  // Allow occasional asteroid belts in inner system.
  if (aAU < snowLineAU && rng.chance(0.10)) {
    return stellar::sim::PlanetType::AsteroidBelt;
  }

  if (aAU < snowLineAU) {
    return stellar::sim::PlanetType::Rocky;
  }

  // Beyond snow line: mix of ices and giants
  const double far = aAU / std::max(0.1, snowLineAU);

  const double gasChance = std::clamp(0.10 + 0.15 * (far - 1.0), 0.10, 0.65);
  if (rng.chance(gasChance)) {
    // Split between gas and ice giants
    return rng.chance(0.75) ? stellar::sim::PlanetType::GasGiant : stellar::sim::PlanetType::IceGiant;
  }

  return stellar::sim::PlanetType::Ice;
}

double rockyRadiusKmFromMass(double massEarth) {
  // Very rough scaling (for small rocky planets)
  const double m = std::max(0.05, massEarth);
  return 6371.0 * std::pow(m, 0.27);
}

double icyRadiusKmFromMass(double massEarth) {
  const double m = std::max(0.05, massEarth);
  return 6371.0 * 1.10 * std::pow(m, 0.28);
}

} // namespace

SystemGenerator::SystemGenerator(SystemGenConfig cfg)
  : m_cfg(cfg)
  , m_names(stellar::core::deriveSeed(cfg.galaxySeed, "names"))
{}

stellar::sim::StarSystem SystemGenerator::generate(const GalaxySystemStub& stub) const {
  stellar::sim::StarSystem sys;
  sys.id = stub.id;
  sys.positionLy = stub.positionLy;

  // Seed per-system generator from galaxy seed + stub id.
  stellar::core::SplitMix64 rng(stellar::core::deriveSeed(m_cfg.galaxySeed, stub.id));

  // --- Star ---
  sys.primary.starClass = pickStarClass(rng);
  const auto spec = specFor(sys.primary.starClass);

  sys.primary.massSolar = rng.uniform(spec.massMin, spec.massMax);
  sys.primary.radiusSolar = approxRadiusSolar(sys.primary.massSolar);
  sys.primary.luminositySolar = approxLuminositySolar(sys.primary.massSolar);
  sys.primary.temperatureK = approxTemperatureK(sys.primary.luminositySolar, sys.primary.radiusSolar);

  sys.primary.name = m_names.makeName(rng, 2, 4);

  // --- Planets ---
  const int maxPlanets = (sys.primary.massSolar < 0.35) ? 8 : 12;
  const int targetPlanets = rng.range<int>(0, maxPlanets);

  const double snowLineAU = 2.7 * std::sqrt(std::max(0.001, sys.primary.luminositySolar));
  double aAU = rng.uniform(0.15, 0.45) * std::sqrt(std::max(0.05, sys.primary.luminositySolar));

  bool madeBelt = false;

  for (int i = 0; i < targetPlanets; ++i) {
    if (i > 0) {
      aAU *= rng.uniform(1.35, 2.20);
    }

    // Stop if we're too far out for this lightweight model.
    if (aAU > 80.0) break;

    auto type = pickPlanetType(aAU, snowLineAU, rng);

    // Keep asteroid belts rare (at most one per system for now)
    if (type == stellar::sim::PlanetType::AsteroidBelt) {
      if (madeBelt) type = stellar::sim::PlanetType::Rocky;
      madeBelt = true;
    }

    stellar::sim::Planet p;
    p.type = type;

    // Rough mass / radius distributions
    switch (p.type) {
      case stellar::sim::PlanetType::Rocky: {
        // Bias toward small planets
        const double x = rng.nextDouble01();
        p.massEarth = 0.08 + 8.0 * x * x; // square biases small
        p.radiusKm = rockyRadiusKmFromMass(p.massEarth);
        break;
      }
      case stellar::sim::PlanetType::Ice: {
        const double x = rng.nextDouble01();
        p.massEarth = 0.1 + 20.0 * x * x;
        p.radiusKm = icyRadiusKmFromMass(p.massEarth);
        break;
      }
      case stellar::sim::PlanetType::GasGiant: {
        const double x = rng.nextDouble01();
        p.massEarth = 40.0 + 300.0 * x;
        // Radius saturates for big gas giants; keep it in a plausible range.
        p.radiusKm = 6371.0 * rng.uniform(9.0, 15.0);
        break;
      }
      case stellar::sim::PlanetType::IceGiant: {
        const double x = rng.nextDouble01();
        p.massEarth = 10.0 + 60.0 * x;
        p.radiusKm = 6371.0 * rng.uniform(3.0, 6.0);
        break;
      }
      case stellar::sim::PlanetType::AsteroidBelt: {
        p.massEarth = 0.0001 * rng.uniform(0.1, 5.0);
        p.radiusKm = 0.0;
        break;
      }
    }

    // Orbit
    p.orbit.semiMajorAxisAU = aAU;
    p.orbit.eccentricity = rng.uniform(0.0, 0.20);
    p.orbit.inclinationDeg = rng.uniform(0.0, 5.0);
    p.orbit.longitudeAscendingNodeDeg = rng.uniform(0.0, 360.0);
    p.orbit.argumentPeriapsisDeg = rng.uniform(0.0, 360.0);
    p.orbit.meanAnomalyAtEpochDeg = rng.uniform(0.0, 360.0);

    const double periodYears = std::sqrt((aAU * aAU * aAU) / std::max(0.08, sys.primary.massSolar));
    p.orbit.orbitalPeriodDays = periodYears * 365.25;

    // Equilibrium temperature (very rough)
    p.equilibriumTempK = 278.0 * std::pow(std::max(0.001, sys.primary.luminositySolar), 0.25) / std::sqrt(std::max(0.05, aAU));

    // Name: StarName + Roman numeral
    p.name = sys.primary.name + " " + romanNumeral(i + 1);

    sys.planets.push_back(std::move(p));
  }

  return sys;
}

} // namespace stellar::proc
