#include "stellar/core/Args.h"
#include "stellar/core/JsonWriter.h"
#include "stellar/core/Log.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/sim/MissionLogic.h"
#include "stellar/sim/NavRoute.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

using namespace stellar;

static const char* starClassName(sim::StarClass c) {
  switch (c) {
    case sim::StarClass::O: return "O";
    case sim::StarClass::B: return "B";
    case sim::StarClass::A: return "A";
    case sim::StarClass::F: return "F";
    case sim::StarClass::G: return "G";
    case sim::StarClass::K: return "K";
    case sim::StarClass::M: return "M";
    default: return "?";
  }
}

static const char* missionTypeName(sim::MissionType t) {
  using sim::MissionType;
  switch (t) {
    case MissionType::Courier: return "Courier";
    case MissionType::Delivery: return "Delivery";
    case MissionType::BountyScan: return "BountyScan";
    case MissionType::BountyKill: return "BountyKill";
    case MissionType::MultiDelivery: return "MultiDelivery";
    case MissionType::Passenger: return "Passenger";
    case MissionType::Smuggle: return "Smuggle";
    default: return "Unknown";
  }
}

static void printMission(const sim::Mission& m) {
  std::cout << "  [" << missionTypeName(m.type) << "] "
            << "toSystem=" << m.toSystem << " toStation=" << m.toStation;
  if (m.type == sim::MissionType::Delivery || m.type == sim::MissionType::MultiDelivery || m.type == sim::MissionType::Smuggle) {
    std::cout << " cargo=" << econ::commodityCode(m.commodity) << " units=" << std::fixed << std::setprecision(0)
              << m.units;
  }
  if (m.type == sim::MissionType::Passenger) {
    std::cout << " party=" << std::fixed << std::setprecision(0) << m.units << " seats";
  }
  if (m.type == sim::MissionType::BountyScan || m.type == sim::MissionType::BountyKill) {
    std::cout << " targetNpcId=" << m.targetNpcId;
  }
  if (m.type == sim::MissionType::MultiDelivery && m.viaStation != 0) {
    std::cout << " via=" << m.viaSystem << "/" << m.viaStation << " leg=" << (int)m.leg;
  }
  std::cout << " reward=" << std::fixed << std::setprecision(0) << m.reward
            << " deadline=" << std::fixed << std::setprecision(1) << m.deadlineDay
            << (m.cargoProvided ? " (cargo-provided)" : " (source-cargo)")
            << "\n";
}

static void writeMissionJson(core::JsonWriter& j, const sim::Mission& m) {
  j.beginObject();
  j.key("id"); j.value((unsigned long long)m.id);
  j.key("type"); j.value(missionTypeName(m.type));
  j.key("factionId"); j.value((unsigned long long)m.factionId);
  j.key("fromSystem"); j.value((unsigned long long)m.fromSystem);
  j.key("fromStation"); j.value((unsigned long long)m.fromStation);
  j.key("toSystem"); j.value((unsigned long long)m.toSystem);
  j.key("toStation"); j.value((unsigned long long)m.toStation);
  j.key("reward"); j.value(m.reward);
  j.key("deadlineDay"); j.value(m.deadlineDay);
  j.key("cargoProvided"); j.value(m.cargoProvided);
  j.key("completed"); j.value(m.completed);
  j.key("failed"); j.value(m.failed);
  j.key("commodity"); j.value(econ::commodityCode(m.commodity));
  j.key("units"); j.value(m.units);
  j.key("viaSystem"); j.value((unsigned long long)m.viaSystem);
  j.key("viaStation"); j.value((unsigned long long)m.viaStation);
  j.key("leg"); j.value((int)m.leg);
  j.key("targetNpcId"); j.value((unsigned long long)m.targetNpcId);
  j.key("scanned"); j.value(m.scanned);
  j.endObject();
}

static void printHelp() {
  std::cout << "stellar_sandbox\n"
            << "  --seed <u64>           Galaxy seed (default: 1337)\n"
            << "  --pos <x y z>          Query position in ly (default: 0 0 0)\n"
            << "  --radius <ly>          Query radius in ly (default: 50)\n"
            << "  --limit <n>            Max systems (default: 32)\n"
            << "  --day <days>           Economy time (days) for trade quotes (default: 0)\n"
            << "  --json                 Emit machine-readable JSON to stdout (also works with --out)\n"
            << "  --out <path>           Write JSON output to a file instead of stdout\n"
            << "\n"
            << "Trade scanning (headless route planner):\n"
            << "  --trade                Print best trade opportunities from a chosen station\n"
            << "  --fromSys <idx>        Origin system index in the printed list (default: 0)\n"
            << "  --fromStation <idx>    Origin station index inside that system (default: 0)\n"
            << "  --cargoKg <kg>         Cargo capacity (kg) used for trip profit (default: 420)\n"
            << "  --tradeLimit <n>       Max trade opportunities to print (default: 12)\n"
            << "\n"
            << "Jump route planning (A*):\n"
            << "  --route                Plot a jump route between two systems in the printed list\n"
            << "  --toSys <idx>          Destination system index in the printed list (default: 0)\n"
            << "  --jr <ly>              Jump range in ly (default: 18)\n"
            << "\n"
            << "Mission board (headless):\n"
            << "  --missions             Print mission offers for the chosen station\n"
            << "  --rep <r>              Player reputation used for offer count/rewards (default: 0)\n"
            << "  --credits <cr>         Starting credits when simulating acceptance (default: 10000)\n"
            << "  --seats <n>            Passenger seats when simulating acceptance (default: 2)\n"
            << "  --acceptOffer <idx>    Accept offer index and print resulting state\n"
            << "  --advanceDays <d>      Advance time by d days before auto-complete (default: 0)\n"
            << "  --autoComplete         After accepting, attempt to complete at destination\n"
            << "\n"
            << "Save/load (tooling):\n"
            << "  --load <path>          Load a save file before running missions\n"
            << "  --save <path>          Save the resulting state to a save file\n";
}

int main(int argc, char** argv) {
  core::setLogLevel(core::LogLevel::Info);

  core::Args args;
  args.setArity("pos", 3);
  args.parse(argc, argv);

  if (args.hasFlag("help") || args.hasFlag("h")) {
    printHelp();
    return 0;
  }

  core::u64 seed = 1337;
  {
    unsigned long long s = (unsigned long long)seed;
    (void)args.getU64("seed", s);
    seed = (core::u64)s;
  }

  math::Vec3d posLy{0,0,0};
  {
    const auto p = args.values("pos");
    if (p.size() >= 3) {
      posLy.x = std::atof(p[p.size() - 3].c_str());
      posLy.y = std::atof(p[p.size() - 2].c_str());
      posLy.z = std::atof(p[p.size() - 1].c_str());
    }
  }

  double radiusLy = 50.0;
  (void)args.getDouble("radius", radiusLy);

  std::size_t limit = 32;
  {
    unsigned long long v = 0;
    if (args.getU64("limit", v)) limit = (std::size_t)v;
  }

  double timeDays = 0.0;
  (void)args.getDouble("day", timeDays);

  const bool json = args.hasFlag("json");
  std::string outPath;
  (void)args.getString("out", outPath);

  const bool doTrade = args.hasFlag("trade");
  std::size_t fromSysIdx = 0;
  std::size_t fromStationIdx = 0;
  {
    unsigned long long v = 0;
    if (args.getU64("fromSys", v)) fromSysIdx = (std::size_t)v;
    if (args.getU64("fromStation", v)) fromStationIdx = (std::size_t)v;
  }
  double cargoCapacityKg = 420.0;
  (void)args.getDouble("cargoKg", cargoCapacityKg);
  std::size_t tradeLimit = 12;
  {
    unsigned long long v = 0;
    if (args.getU64("tradeLimit", v)) tradeLimit = (std::size_t)v;
  }

  const bool doRoute = args.hasFlag("route");
  std::size_t toSysIdx = 0;
  {
    unsigned long long v = 0;
    if (args.getU64("toSys", v)) toSysIdx = (std::size_t)v;
  }
  double jumpRangeLy = 18.0;
  (void)args.getDouble("jr", jumpRangeLy);

  const bool doMissions = args.hasFlag("missions");
  double rep = 0.0;
  (void)args.getDouble("rep", rep);
  double credits = 10000.0;
  (void)args.getDouble("credits", credits);
  int seats = 2;
  (void)args.getInt("seats", seats);
  int acceptOffer = -1;
  (void)args.getInt("acceptOffer", acceptOffer);
  double advanceDays = 0.0;
  (void)args.getDouble("advanceDays", advanceDays);
  const bool autoComplete = args.hasFlag("autoComplete");
  std::string loadPath;
  std::string savePath;
  (void)args.getString("load", loadPath);
  (void)args.getString("save", savePath);

  sim::Universe u(seed);

  const auto systems = u.queryNearby(posLy, radiusLy, limit);

  // If JSON is requested, we build up a machine-readable object.
  // Otherwise, we keep the original human-readable console output.
  std::unique_ptr<std::ofstream> jsonFile;
  std::ostream* jsonStream = &std::cout;
  if (json && !outPath.empty()) {
    jsonFile = std::make_unique<std::ofstream>(outPath, std::ios::out | std::ios::trunc);
    if (!*jsonFile) {
      std::cerr << "Failed to open --out file: " << outPath << "\n";
      return 1;
    }
    jsonStream = jsonFile.get();
  }

  core::JsonWriter j(*jsonStream, /*pretty=*/true);
  if (json) {
    j.beginObject();
    j.key("seed"); j.value((unsigned long long)seed);
    j.key("query");
    j.beginObject();
    j.key("posLy");
    j.beginArray();
    j.value(posLy.x); j.value(posLy.y); j.value(posLy.z);
    j.endArray();
    j.key("radiusLy"); j.value(radiusLy);
    j.key("limit"); j.value((unsigned long long)limit);
    j.key("day"); j.value(timeDays);
    j.endObject();
    j.key("systems");
    j.beginArray();
    for (const auto& s : systems) {
      const double distLy = (s.posLy - posLy).length();
      j.beginObject();
      j.key("id"); j.value((unsigned long long)s.id);
      j.key("name"); j.value(s.name);
      j.key("class"); j.value(starClassName(s.primaryClass));
      j.key("distLy"); j.value(distLy);
      j.key("planetCount"); j.value((unsigned long long)s.planetCount);
      j.key("stationCount"); j.value((unsigned long long)s.stationCount);
      j.key("factionId"); j.value((unsigned long long)s.factionId);
      j.endObject();
    }
    j.endArray();
  } else {
    std::cout << "Seed: " << seed << "\n";
    std::cout << "Query @ (" << posLy.x << "," << posLy.y << "," << posLy.z << ") radius=" << radiusLy << " ly\n";
    std::cout << "Found " << systems.size() << " systems\n\n";

    for (const auto& s : systems) {
      const math::Vec3d d = s.posLy - posLy;
      const double dist = std::sqrt(d.lengthSq());

      std::cout << std::setw(14) << s.id
                << "  " << std::setw(14) << s.name
                << "  class=" << starClassName(s.primaryClass)
                << "  dist=" << std::fixed << std::setprecision(2) << dist << " ly"
                << "  planets=" << s.planetCount
                << "  stations=" << s.stationCount
                << "  faction=" << s.factionId
                << "\n";
    }
  }

  if (!systems.empty() && !json) {
    std::cout << "\n--- Example system detail: " << systems.front().name << " ---\n";
    const auto& sys = u.getSystem(systems.front().id, &systems.front());
    std::cout << "Star mass=" << sys.star.massSol << " Sol, lum=" << sys.star.luminositySol << " Sol\n";
    std::cout << "Planets:\n";
    for (const auto& p : sys.planets) {
      std::cout << "  - " << p.name
                << " a=" << p.orbit.semiMajorAxisAU << " AU"
                << " e=" << p.orbit.eccentricity
                << " period=" << p.orbit.periodDays << " days\n";
    }
    std::cout << "Stations:\n";
    for (const auto& st : sys.stations) {
      std::cout << "  - " << st.name << " (fee=" << st.feeRate << ")\n";
    }
  }

  if (doTrade && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[trade] origin system has no stations\n";
      return 1;
    }
    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& fromSt = fromSys.stations[fromStationIdx];

    auto& fromEcon = u.stationEconomy(fromSt, timeDays);

    struct Idea {
      sim::SystemId sysId{0};
      sim::StationId stationId{0};
      std::string sysName;
      std::string stationName;
      econ::CommodityId commodity{econ::CommodityId::Food};
      double buyAsk{0.0};
      double sellBid{0.0};
      double units{0.0};
      double netProfitPerUnit{0.0};
      double netTripProfit{0.0};
      double distanceLy{0.0};
    };

    std::vector<Idea> ideas;
    ideas.reserve(systems.size());

    for (const auto& stub : systems) {
      const auto& sys = u.getSystem(stub.id, &stub);
      const double distLy = (stub.posLy - fromStub.posLy).length();

      for (const auto& toSt : sys.stations) {
        if (stub.id == fromStub.id && toSt.id == fromSt.id) continue;

        auto& toEcon = u.stationEconomy(toSt, timeDays);
        const auto routes = econ::bestRoutesForCargo(fromEcon,
                                                     fromSt.economyModel,
                                                     toEcon,
                                                     toSt.economyModel,
                                                     cargoCapacityKg,
                                                     fromSt.feeRate,
                                                     toSt.feeRate,
                                                     0.10,
                                                     1);

        if (routes.empty()) continue;
        const auto& r = routes.front();

        Idea it;
        it.sysId = stub.id;
        it.stationId = toSt.id;
        it.sysName = stub.name;
        it.stationName = toSt.name;
        it.commodity = r.commodity;
        it.buyAsk = r.buyPrice;
        it.sellBid = r.sellPrice;
        it.units = r.unitsPossible;
        it.netProfitPerUnit = r.netProfitPerUnit;
        it.netTripProfit = r.netProfitTotal;
        it.distanceLy = distLy;
        ideas.push_back(std::move(it));
      }
    }

    std::sort(ideas.begin(), ideas.end(), [](const Idea& a, const Idea& b) {
      if (a.netTripProfit != b.netTripProfit) return a.netTripProfit > b.netTripProfit;
      return a.distanceLy < b.distanceLy;
    });

    if (ideas.size() > tradeLimit) ideas.resize(tradeLimit);

    if (json) {
      j.key("trade");
      j.beginObject();
      j.key("day"); j.value(timeDays);
      j.key("cargoKg"); j.value(cargoCapacityKg);
      j.key("from");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationId"); j.value((unsigned long long)fromSt.id);
      j.key("stationName"); j.value(fromSt.name);
      j.key("feeRate"); j.value(fromSt.feeRate);
      j.endObject();
      j.key("ideas");
      j.beginArray();
      for (const auto& it : ideas) {
        j.beginObject();
        j.key("commodity"); j.value(econ::commodityCode(it.commodity));
        j.key("toSystemId"); j.value((unsigned long long)it.sysId);
        j.key("toStationId"); j.value((unsigned long long)it.stationId);
        j.key("toSystemName"); j.value(it.sysName);
        j.key("toStationName"); j.value(it.stationName);
        j.key("distLy"); j.value(it.distanceLy);
        j.key("units"); j.value(it.units);
        j.key("netProfitPerUnit"); j.value(it.netProfitPerUnit);
        j.key("netTripProfit"); j.value(it.netTripProfit);
        j.key("buyAsk"); j.value(it.buyAsk);
        j.key("sellBid"); j.value(it.sellBid);
        j.endObject();
      }
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Trade scan (day=" << timeDays << ", cargoKg=" << cargoCapacityKg << ") ---\n";
      std::cout << "From: " << fromStub.name << " / " << fromSt.name << "  (fee=" << fromSt.feeRate * 100.0 << "%)\n";
      if (ideas.empty()) {
        std::cout << "No profitable routes found inside the query set.\n";
      } else {
        for (const auto& it : ideas) {
          const auto code = econ::commodityCode(it.commodity);
          std::cout << std::setw(6) << code
                    << "  to=" << it.sysName << "/" << it.stationName
                    << "  dist=" << std::fixed << std::setprecision(1) << it.distanceLy << " ly"
                    << "  units=" << std::fixed << std::setprecision(0) << it.units
                    << "  net/unit=" << std::fixed << std::setprecision(2) << it.netProfitPerUnit
                    << "  net/trip=" << std::fixed << std::setprecision(0) << it.netTripProfit
                    << "\n";
        }
        std::cout << "Tip: increase --radius (or move --pos) to scan a larger area.\n";
      }
    }
  }

  if (doRoute && !systems.empty()) {
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);
    toSysIdx = std::min(toSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& toStub = systems[toSysIdx];

    sim::RoutePlanStats stats{};
    const auto route = sim::plotRouteAStarHops(systems, fromStub.id, toStub.id, jumpRangeLy, &stats);

    if (json) {
      j.key("route");
      j.beginObject();
      j.key("fromSysIdx"); j.value((unsigned long long)fromSysIdx);
      j.key("toSysIdx"); j.value((unsigned long long)toSysIdx);
      j.key("fromSystemId"); j.value((unsigned long long)fromStub.id);
      j.key("toSystemId"); j.value((unsigned long long)toStub.id);
      j.key("jumpRangeLy"); j.value(jumpRangeLy);
      j.key("found"); j.value(!route.empty());
      j.key("hops"); j.value((unsigned long long)(route.size() > 1 ? route.size() - 1 : 0));
      j.key("distanceLy"); j.value(stats.distanceLy);
      j.key("visited"); j.value((unsigned long long)stats.visited);
      j.key("expansions"); j.value((unsigned long long)stats.expansions);
      j.key("path");
      j.beginArray();
      for (const auto id : route) j.value((unsigned long long)id);
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Route plan (A* hops, jr=" << std::fixed << std::setprecision(1) << jumpRangeLy << " ly) ---\n";
      std::cout << "From: [" << fromSysIdx << "] " << fromStub.name << " (" << fromStub.id << ")\n";
      std::cout << "To:   [" << toSysIdx << "] " << toStub.name << " (" << toStub.id << ")\n";

      if (route.empty()) {
        std::cout << "No route found inside the queried node set.\n";
        std::cout << "Tip: increase --radius (and/or --limit) so intermediate systems are available.\n";
      } else {
        std::unordered_map<sim::SystemId, const sim::SystemStub*> stubById;
        stubById.reserve(systems.size());
        for (const auto& s : systems) stubById[s.id] = &s;

        std::cout << "Hops: " << (route.size() - 1)
                  << "  Dist: " << std::fixed << std::setprecision(2) << stats.distanceLy
                  << " ly  (visited=" << stats.visited << ")\n";

        for (std::size_t i = 0; i < route.size(); ++i) {
          const auto id = route[i];
          const auto it = stubById.find(id);
          const auto* s = (it != stubById.end()) ? it->second : nullptr;

          std::cout << "  " << std::setw(2) << i << ": ";
          if (s) {
            std::cout << s->name << " (" << s->id << ")";
          } else {
            std::cout << "(unknown " << (unsigned long long)id << ")";
          }

          if (i > 0 && s) {
            const auto itPrev = stubById.find(route[i - 1]);
            const auto* p = (itPrev != stubById.end()) ? itPrev->second : nullptr;
            if (p) {
              const double d = (s->posLy - p->posLy).length();
              std::cout << "  +" << std::fixed << std::setprecision(2) << d << " ly";
            }
          }
          std::cout << "\n";
        }
      }
    }
  }

  if (doMissions && !systems.empty()) {
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);
    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[missions] chosen system has no stations\n";
      if (json) {
        j.key("missions");
        j.beginObject();
        j.key("error"); j.value("chosen system has no stations");
        j.endObject();
      }
    } else {
      fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
      const auto& dockedStation = fromSys.stations[fromStationIdx];

      sim::SaveGame save{};
      if (!loadPath.empty()) {
        if (!sim::loadFromFile(loadPath, save)) {
          std::cerr << "Failed to load save from: " << loadPath << "\n";
          return 1;
        }
      }

      // Tool defaults / overrides.
      save.seed = seed;
      save.timeDays = timeDays;
      save.currentSystem = fromStub.id;
      save.dockedStation = dockedStation.id;
      save.credits = credits;
      save.cargoCapacityKg = cargoCapacityKg;
      save.passengerSeats = std::max(0, seats);

      // Ensure offers exist + are deterministic.
      sim::refreshMissionOffers(u, fromSys, dockedStation, timeDays, rep, save);

      if (json) {
        j.key("missions");
        j.beginObject();
        j.key("station");
        j.beginObject();
        j.key("systemId"); j.value((unsigned long long)fromStub.id);
        j.key("systemName"); j.value(fromStub.name);
        j.key("stationId"); j.value((unsigned long long)dockedStation.id);
        j.key("stationName"); j.value(dockedStation.name);
        j.key("rep"); j.value(rep);
        j.endObject();
        j.key("offers");
        j.beginArray();
        for (const auto& m : save.missionOffers) writeMissionJson(j, m);
        j.endArray();
      } else {
        std::cout << "\n--- Mission board (day=" << timeDays << ", rep=" << rep << ") ---\n";
        std::cout << "At: " << fromStub.name << " / " << dockedStation.name << "\n";
        if (save.missionOffers.empty()) {
          std::cout << "No mission offers.\n";
        } else {
          for (std::size_t i = 0; i < save.missionOffers.size(); ++i) {
            std::cout << "Offer #" << i << ":\n";
            printMission(save.missionOffers[i]);
          }
        }
      }

      std::string acceptError;
      bool accepted = false;
      sim::Mission acceptedMission{};

      if (acceptOffer >= 0) {
        const std::size_t idx = (std::size_t)acceptOffer;
        if (idx >= save.missionOffers.size()) {
          acceptError = "acceptOffer index out of range";
        } else {
          acceptedMission = save.missionOffers[idx];
          accepted = sim::acceptMissionOffer(u, dockedStation, timeDays, rep, save, idx, &acceptError);
        }
      }

      if (json) {
        if (acceptOffer >= 0) {
          j.key("accept");
          j.beginObject();
          j.key("ok"); j.value(accepted);
          if (!acceptError.empty()) {
            j.key("error"); j.value(acceptError);
          }
          j.key("credits"); j.value(save.credits);
          j.key("cargoCapacityKg"); j.value(save.cargoCapacityKg);
          j.key("passengerSeats"); j.value(save.passengerSeats);
          j.key("activeMissions");
          j.beginArray();
          for (const auto& m : save.missions) writeMissionJson(j, m);
          j.endArray();
          j.endObject();
        }

        if (accepted && autoComplete) {
          // Best-effort completion simulation.
          sim::MissionDockResult r{};
          save.timeDays += advanceDays;

          const auto completeAt = [&](sim::SystemId sysId, sim::StationId stId) {
            const auto& sys = u.getSystem(sysId, nullptr);
            const sim::Station* st = nullptr;
            for (const auto& s : sys.stations) {
              if (s.id == stId) { st = &s; break; }
            }
            if (!st) return;
            save.currentSystem = sysId;
            save.dockedStation = stId;
            r = sim::tryCompleteMissionsAtDock(u, sys, *st, save.timeDays, save);
          };

          if (acceptedMission.type == sim::MissionType::MultiDelivery && acceptedMission.viaStation != 0) {
            completeAt(acceptedMission.viaSystem, acceptedMission.viaStation);
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          } else {
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          }

          j.key("autoComplete");
          j.beginObject();
          j.key("advancedDays"); j.value(advanceDays);
          j.key("completed"); j.value(r.completed);
          j.key("progressedMultiLeg"); j.value(r.progressedMultiLeg);
          j.key("credits"); j.value(save.credits);
          j.endObject();
        }

        j.endObject(); // missions
      } else {
        if (acceptOffer >= 0) {
          if (!accepted) {
            std::cout << "\n[accept] failed: " << acceptError << "\n";
          } else {
            std::cout << "\n[accept] ok. credits=" << std::fixed << std::setprecision(0) << save.credits
                      << " activeMissions=" << save.missions.size() << "\n";
          }
        }

        if (accepted && autoComplete) {
          save.timeDays += advanceDays;
          sim::MissionDockResult r{};
          const auto completeAt = [&](sim::SystemId sysId, sim::StationId stId) {
            const auto& sys = u.getSystem(sysId, nullptr);
            const sim::Station* st = nullptr;
            for (const auto& s : sys.stations) {
              if (s.id == stId) { st = &s; break; }
            }
            if (!st) {
              std::cout << "[autoComplete] station not found: " << stId << "\n";
              return;
            }
            save.currentSystem = sysId;
            save.dockedStation = stId;
            r = sim::tryCompleteMissionsAtDock(u, sys, *st, save.timeDays, save);
          };

          if (acceptedMission.type == sim::MissionType::MultiDelivery && acceptedMission.viaStation != 0) {
            completeAt(acceptedMission.viaSystem, acceptedMission.viaStation);
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          } else {
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          }

          std::cout << "[autoComplete] advancedDays=" << advanceDays
                    << " completed=" << r.completed
                    << " progressedMultiLeg=" << r.progressedMultiLeg
                    << " credits=" << std::fixed << std::setprecision(0) << save.credits
                    << "\n";
        }
      }

      if (!savePath.empty()) {
        if (!sim::saveToFile(save, savePath)) {
          std::cerr << "Failed to write save to: " << savePath << "\n";
          return 1;
        }
        if (!json) {
          std::cout << "[save] wrote " << savePath << "\n";
        }
      }
    }
  }

  if (json) {
    j.endObject();
    *jsonStream << "\n";
  }

  return 0;
}
