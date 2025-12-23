#include "stellar/core/Log.h"
#include "stellar/sim/Universe.h"

#include <charconv>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>

namespace {

void printHelp() {
  std::cout <<
    "stellar_sandbox — Procedural space sim starter\n"
    "\n"
    "Usage:\n"
    "  stellar_sandbox [options]\n"
    "\n"
    "Options:\n"
    "  --help                Show this help\n"
    "  --seed <u64>           Galaxy seed (default: 1)\n"
    "  --systems <n>          Number of systems to generate (default: 20)\n"
    "  --list                List system summaries\n"
    "  --pick <index>         Print details for one system (0-based)\n"
    "  --time <days>          Show planet positions at this time (default: 0)\n"
    "  --log <level>          trace|debug|info|warn|error|off (default: info)\n"
    "\n";
}

template <typename T>
std::optional<T> parseNumber(std::string_view s) {
  T value{};
  auto first = s.data();
  auto last = s.data() + s.size();
  auto [ptr, ec] = std::from_chars(first, last, value);
  if (ec != std::errc{} || ptr != last) return std::nullopt;
  return value;
}

std::optional<stellar::core::LogLevel> parseLogLevel(std::string_view s) {
  using stellar::core::LogLevel;
  if (s == "trace") return LogLevel::Trace;
  if (s == "debug") return LogLevel::Debug;
  if (s == "info")  return LogLevel::Info;
  if (s == "warn")  return LogLevel::Warn;
  if (s == "error") return LogLevel::Error;
  if (s == "off")   return LogLevel::Off;
  return std::nullopt;
}

void printSystemSummary(std::size_t index, const stellar::sim::StarSystem& sys) {
  std::cout
    << "[" << index << "] "
    << sys.primary.name << " (class " << stellar::sim::toString(sys.primary.starClass) << ")"
    << " planets=" << sys.planets.size()
    << " posLy=" << sys.positionLy
    << " id=" << sys.id
    << "\n";
}

void printSystemDetails(const stellar::sim::StarSystem& sys, double timeDays) {
  using std::cout;

  cout << "\n=== Star System: " << sys.primary.name << " ===\n";
  cout << "ID: " << sys.id << "\n";
  cout << "Position (ly): " << sys.positionLy << "\n";
  cout << "\n-- Star --\n";
  cout << "Class: " << stellar::sim::toString(sys.primary.starClass) << "\n";
  cout << std::fixed << std::setprecision(3);
  cout << "Mass (solar): " << sys.primary.massSolar << "\n";
  cout << "Radius (solar): " << sys.primary.radiusSolar << "\n";
  cout << "Luminosity (solar): " << sys.primary.luminositySolar << "\n";
  cout << "Temperature (K): " << sys.primary.temperatureK << "\n";

  cout << "\n-- Planets (" << sys.planets.size() << ") --\n";
  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    const auto& p = sys.planets[i];
    cout << "\n[" << i << "] " << p.name << "\n";
    cout << "Type: " << stellar::sim::toString(p.type) << "\n";
    cout << "Mass (Earth): " << p.massEarth << "\n";
    cout << "Radius (km): " << p.radiusKm << "\n";
    cout << "Semi-major axis (AU): " << p.orbit.semiMajorAxisAU << "\n";
    cout << "Eccentricity: " << p.orbit.eccentricity << "\n";
    cout << "Period (days): " << p.orbit.orbitalPeriodDays << "\n";
    cout << "Equilibrium T (K): " << p.equilibriumTempK << "\n";

    const auto pos = p.orbit.positionAU(timeDays);
    cout << "Position at t=" << timeDays << " days (AU): " << pos << "\n";
  }
  cout << "\n";
}

} // namespace

int main(int argc, char** argv) {
  stellar::core::setLogLevel(stellar::core::LogLevel::Info);

  stellar::sim::UniverseConfig cfg;
  cfg.seed = 1;
  cfg.systemCount = 20;

  bool list = false;
  std::optional<std::size_t> pick;
  double timeDays = 0.0;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];

    auto requireValue = [&](std::string_view opt) -> std::optional<std::string_view> {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << opt << "\n";
        return std::nullopt;
      }
      return std::string_view(argv[++i]);
    };

    if (arg == "--help" || arg == "-h") {
      printHelp();
      return 0;
    } else if (arg == "--seed") {
      const auto v = requireValue(arg);
      if (!v) return 2;
      const auto n = parseNumber<stellar::core::u64>(*v);
      if (!n) { std::cerr << "Invalid seed: " << *v << "\n"; return 2; }
      cfg.seed = *n;
    } else if (arg == "--systems") {
      const auto v = requireValue(arg);
      if (!v) return 2;
      const auto n = parseNumber<std::size_t>(*v);
      if (!n) { std::cerr << "Invalid systems count: " << *v << "\n"; return 2; }
      cfg.systemCount = *n;
    } else if (arg == "--list") {
      list = true;
    } else if (arg == "--pick") {
      const auto v = requireValue(arg);
      if (!v) return 2;
      const auto n = parseNumber<std::size_t>(*v);
      if (!n) { std::cerr << "Invalid index: " << *v << "\n"; return 2; }
      pick = *n;
    } else if (arg == "--time") {
      const auto v = requireValue(arg);
      if (!v) return 2;
      const auto n = parseNumber<double>(*v);
      if (!n) { std::cerr << "Invalid time: " << *v << "\n"; return 2; }
      timeDays = *n;
    } else if (arg == "--log") {
      const auto v = requireValue(arg);
      if (!v) return 2;
      const auto lvl = parseLogLevel(*v);
      if (!lvl) { std::cerr << "Invalid log level: " << *v << "\n"; return 2; }
      stellar::core::setLogLevel(*lvl);
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      printHelp();
      return 2;
    }
  }

  stellar::sim::Universe universe(cfg);
  universe.generate();

  if (universe.size() == 0) {
    std::cerr << "No systems generated.\n";
    return 1;
  }

  if (list) {
    for (std::size_t i = 0; i < universe.size(); ++i) {
      printSystemSummary(i, universe.system(i));
    }
  }

  if (pick) {
    if (*pick >= universe.size()) {
      std::cerr << "Pick index out of range: " << *pick << " (size=" << universe.size() << ")\n";
      return 2;
    }
    printSystemDetails(universe.system(*pick), timeDays);
  }

  if (!list && !pick) {
    // Default behavior: print the first system.
    printSystemDetails(universe.system(0), timeDays);
  }

  return 0;
}
