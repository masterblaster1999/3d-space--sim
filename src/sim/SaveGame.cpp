#include "stellar/sim/SaveGame.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <fstream>
#include <sstream>

namespace stellar::sim {

bool saveToFile(const SaveGame& s, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Error, "SaveGame: failed to open file for writing: " + path);
    return false;
  }

  f.setf(std::ios::fixed);
  f.precision(8);

  f << "StellarForgeSave " << s.version << "\n";
  f << "seed " << s.seed << "\n";
  f << "timeDays " << s.timeDays << "\n";
  f << "currentSystem " << s.currentSystem << "\n";
  f << "dockedStation " << s.dockedStation << "\n";

  f << "shipPosKm " << s.shipPosKm.x << " " << s.shipPosKm.y << " " << s.shipPosKm.z << "\n";
  f << "shipVelKmS " << s.shipVelKmS.x << " " << s.shipVelKmS.y << " " << s.shipVelKmS.z << "\n";
  f << "shipOrient " << s.shipOrient.w << " " << s.shipOrient.x << " " << s.shipOrient.y << " " << s.shipOrient.z << "\n";
  f << "shipAngVel " << s.shipAngVelRadS.x << " " << s.shipAngVelRadS.y << " " << s.shipAngVelRadS.z << "\n";

  f << "credits " << s.credits << "\n";

  // Ship meta / progression
  f << "fuel " << s.fuel << "\n";
  f << "fuelMax " << s.fuelMax << "\n";
  f << "fsdRangeLy " << s.fsdRangeLy << "\n";
  f << "hull " << s.hull << "\n";
  f << "shield " << s.shield << "\n";
  f << "heat " << s.heat << "\n";
  f << "cargoCapacityKg " << s.cargoCapacityKg << "\n";
  f << "fsdReadyDay " << s.fsdReadyDay << "\n";

  // Loadout / progression
  f << "shipHull " << (int)s.shipHull << "\n";
  f << "thrusterMk " << (int)s.thrusterMk << "\n";
  f << "shieldMk " << (int)s.shieldMk << "\n";
  f << "distributorMk " << (int)s.distributorMk << "\n";
  f << "weaponPrimary " << (int)s.weaponPrimary << "\n";
  f << "weaponSecondary " << (int)s.weaponSecondary << "\n";

  f << "cargo";
  for (double u : s.cargo) f << " " << u;
  f << "\n";

  // Exploration
  f << "explorationDataCr " << s.explorationDataCr << "\n";
  f << "scannedKeys " << s.scannedKeys.size() << "\n";
  for (core::u64 k : s.scannedKeys) {
    f << "scan " << k << "\n";
  }

  // Missions
  f << "nextMissionId " << s.nextMissionId << "\n";
  f << "missions " << s.missions.size() << "\n";
  for (const auto& m : s.missions) {
    f << "mission "
      << m.id << " "
      << static_cast<int>(m.type) << " "
      << m.fromSystem << " " << m.fromStation << " "
      << m.toSystem << " " << m.toStation << " "
      << static_cast<int>(m.commodity) << " "
      << m.units << " "
      << m.targetNpcId << " "
      << m.reward << " "
      << m.deadlineDay << " "
      << (m.completed ? 1 : 0) << " "
      << (m.failed ? 1 : 0) << " "
      << (m.cargoProvided ? 1 : 0) << " "
      // Optional / newer fields (kept at end for backward compatibility)
      << m.factionId << " "
      << m.viaSystem << " " << m.viaStation << " "
      << static_cast<int>(m.leg) << " "
      << (m.scanned ? 1 : 0)
      << "\n";
  }

  // Mission board (cached offers)
  f << "mission_offers_station " << s.missionOffersStationId << "\n";
  f << "mission_offers_day " << s.missionOffersDayStamp << "\n";
  f << "mission_offers " << s.missionOffers.size() << "\n";
  for (const auto& m : s.missionOffers) {
    f << "offer "
      << m.id << " "
      << static_cast<int>(m.type) << " "
      << m.fromSystem << " " << m.fromStation << " "
      << m.toSystem << " " << m.toStation << " "
      << static_cast<int>(m.commodity) << " "
      << m.units << " "
      << m.targetNpcId << " "
      << m.reward << " "
      << m.deadlineDay << " "
      << (m.completed ? 1 : 0) << " "
      << (m.failed ? 1 : 0) << " "
      << (m.cargoProvided ? 1 : 0) << " "
      // Optional / newer fields (kept at end for backward compatibility)
      << m.factionId << " "
      << m.viaSystem << " " << m.viaStation << " "
      << static_cast<int>(m.leg) << " "
      << (m.scanned ? 1 : 0)
      << "\n";
  }

  // Reputation
  f << "reputation " << s.reputation.size() << "\n";
  for (const auto& r : s.reputation) {
    f << "rep " << r.factionId << " " << r.rep << "\n";
  }

  // Bounties
  f << "bounties " << s.bounties.size() << "\n";
  for (const auto& b : s.bounties) {
    f << "bounty " << b.factionId << " " << b.bountyCr << "\n";
  }

  // Bounty vouchers (earned from destroying criminals; redeemed later)
  f << "bounty_vouchers " << s.bountyVouchers.size() << "\n";
  for (const auto& v : s.bountyVouchers) {
    f << "voucher " << v.factionId << " " << v.bountyCr << "\n";
  }

  f << "station_overrides " << s.stationOverrides.size() << "\n";
  for (const auto& ov : s.stationOverrides) {
    f << "station " << ov.stationId << "\n";
    f << "lastUpdateDay " << ov.state.lastUpdateDay << "\n";
    f << "lastSampleDay " << ov.state.lastSampleDay << "\n";

    f << "inventory";
    for (double v : ov.state.inventory) f << " " << v;
    f << "\n";

    // history per commodity
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      const auto& hist = ov.state.history[i];
      f << "history " << i << " " << hist.size();
      for (const auto& p : hist) f << " " << p.day << " " << p.price;
      f << "\n";
    }

    f << "endstation\n";
  }

  return true;
}

static bool expectToken(std::istream& in, const char* tok) {
  std::string s;
  if (!(in >> s)) return false;
  return s == tok;
}

bool loadFromFile(const std::string& path, SaveGame& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "SaveGame: file not found: " + path);
    return false;
  }

  std::string header;
  if (!(f >> header)) return false;
  if (header != "StellarForgeSave") {
    stellar::core::log(stellar::core::LogLevel::Error, "SaveGame: bad header");
    return false;
  }

  int version = 0;
  if (!(f >> version)) return false;
  out = SaveGame{};
  out.version = version;

  std::string key;
  while (f >> key) {
    if (key == "seed") {
      f >> out.seed;
    } else if (key == "timeDays") {
      f >> out.timeDays;
    } else if (key == "currentSystem") {
      f >> out.currentSystem;
    } else if (key == "dockedStation") {
      f >> out.dockedStation;
    } else if (key == "shipPosKm") {
      f >> out.shipPosKm.x >> out.shipPosKm.y >> out.shipPosKm.z;
    } else if (key == "shipVelKmS") {
      f >> out.shipVelKmS.x >> out.shipVelKmS.y >> out.shipVelKmS.z;
    } else if (key == "shipOrient") {
      f >> out.shipOrient.w >> out.shipOrient.x >> out.shipOrient.y >> out.shipOrient.z;
    } else if (key == "shipAngVel") {
      f >> out.shipAngVelRadS.x >> out.shipAngVelRadS.y >> out.shipAngVelRadS.z;
    } else if (key == "credits") {
      f >> out.credits;
    } else if (key == "fuel") {
      f >> out.fuel;
    } else if (key == "fuelMax") {
      f >> out.fuelMax;
    } else if (key == "fsdRangeLy") {
      f >> out.fsdRangeLy;
    } else if (key == "hull") {
      f >> out.hull;
    } else if (key == "shield") {
      f >> out.shield;
    } else if (key == "heat") {
      f >> out.heat;
    } else if (key == "cargoCapacityKg") {
      f >> out.cargoCapacityKg;
    } else if (key == "fsdReadyDay") {
      f >> out.fsdReadyDay;
    } else if (key == "shipHull") {
      int v = 0;
      f >> v;
      out.shipHull = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "thrusterMk") {
      int v = 0;
      f >> v;
      out.thrusterMk = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "shieldMk") {
      int v = 0;
      f >> v;
      out.shieldMk = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "distributorMk") {
      int v = 0;
      f >> v;
      out.distributorMk = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "weaponPrimary") {
      int v = 0;
      f >> v;
      out.weaponPrimary = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "weaponSecondary") {
      int v = 0;
      f >> v;
      out.weaponSecondary = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "cargo") {
      for (std::size_t i = 0; i < econ::kCommodityCount; ++i) f >> out.cargo[i];
    } else if (key == "explorationDataCr") {
      f >> out.explorationDataCr;
    } else if (key == "scannedKeys") {
      std::size_t n = 0;
      f >> n;
      out.scannedKeys.clear();
      out.scannedKeys.reserve(n);

      // Read scan entries robustly.
      // If the count is stale/corrupt, fall back to stopping when we no longer
      // see a "scan" tag, *without* consuming the next top-level key.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "scan") {
          f.clear();
          f.seekg(pos);
          break;
        }
        core::u64 k = 0;
        if (!(f >> k)) break;
        out.scannedKeys.push_back(k);
      }
    } else if (key == "nextMissionId") {
      f >> out.nextMissionId;
    } else if (key == "missions") {
      std::size_t n = 0;
      f >> n;
      out.missions.clear();
      out.missions.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        std::string tok;
        if (!(f >> tok) || tok != "mission") return false;

        // Parse the rest of the mission line so we can be backward compatible
        // with older saves that did not include newer fields.
        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        Mission m{};
        int type = 0;
        int commodity = 0;
        int completed = 0;
        int failed = 0;
        int cargoProvided = 0;

        iss >> m.id
            >> type
            >> m.fromSystem >> m.fromStation
            >> m.toSystem >> m.toStation
            >> commodity
            >> m.units
            >> m.targetNpcId
            >> m.reward
            >> m.deadlineDay
            >> completed
            >> failed
            >> cargoProvided;

        m.type = static_cast<MissionType>(type);
        m.commodity = static_cast<econ::CommodityId>(commodity);
        m.completed = (completed != 0);
        m.failed = (failed != 0);
        m.cargoProvided = (cargoProvided != 0);

        // Optional fields (save version >= 3)
        int leg = 0;
        int scanned = 0;
        if (iss >> m.factionId >> m.viaSystem >> m.viaStation >> leg >> scanned) {
          m.leg = static_cast<core::u8>(std::clamp(leg, 0, 255));
          m.scanned = (scanned != 0);
        }

        out.missions.push_back(std::move(m));
      }
    } else if (key == "mission_offers_station") {
      f >> out.missionOffersStationId;
    } else if (key == "mission_offers_day") {
      f >> out.missionOffersDayStamp;
    } else if (key == "mission_offers") {
      std::size_t n = 0;
      f >> n;
      out.missionOffers.clear();
      out.missionOffers.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        std::string tok;
        if (!(f >> tok) || tok != "offer") return false;

        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        Mission m{};
        int type = 0;
        int commodity = 0;
        int completed = 0;
        int failed = 0;
        int cargoProvided = 0;

        iss >> m.id
            >> type
            >> m.fromSystem >> m.fromStation
            >> m.toSystem >> m.toStation
            >> commodity
            >> m.units
            >> m.targetNpcId
            >> m.reward
            >> m.deadlineDay
            >> completed
            >> failed
            >> cargoProvided;

        m.type = static_cast<MissionType>(type);
        m.commodity = static_cast<econ::CommodityId>(commodity);
        m.completed = (completed != 0);
        m.failed = (failed != 0);
        m.cargoProvided = (cargoProvided != 0);

        int leg = 0;
        int scanned = 0;
        if (iss >> m.factionId >> m.viaSystem >> m.viaStation >> leg >> scanned) {
          m.leg = static_cast<core::u8>(std::clamp(leg, 0, 255));
          m.scanned = (scanned != 0);
        }

        out.missionOffers.push_back(std::move(m));
      }
} else if (key == "reputation") {
      std::size_t n = 0;
      f >> n;
      out.reputation.clear();
      out.reputation.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        std::string tok;
        if (!(f >> tok) || tok != "rep") return false;
        FactionReputation r{};
        f >> r.factionId >> r.rep;
        out.reputation.push_back(std::move(r));
      }
    } else if (key == "bounties") {
      std::size_t n = 0;
      f >> n;
      out.bounties.clear();
      out.bounties.reserve(n);

      // Same robustness strategy as scannedKeys.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "bounty") {
          f.clear();
          f.seekg(pos);
          break;
        }
        FactionBounty b{};
        if (!(f >> b.factionId >> b.bountyCr)) break;
        out.bounties.push_back(std::move(b));
      }
    } else if (key == "bounty_vouchers") {
      std::size_t n = 0;
      f >> n;
      out.bountyVouchers.clear();
      out.bountyVouchers.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "voucher") {
          f.clear();
          f.seekg(pos);
          break;
        }
        FactionBounty v{};
        if (!(f >> v.factionId >> v.bountyCr)) break;
        out.bountyVouchers.push_back(std::move(v));
      }
    } else if (key == "station_overrides") {
      std::size_t n = 0;
      f >> n;
      out.stationOverrides.clear();
      out.stationOverrides.reserve(n);

      for (std::size_t si = 0; si < n; ++si) {
        StationEconomyOverride ov{};
        if (!expectToken(f, "station")) return false;
        f >> ov.stationId;

        // Parse station block until endstation
        std::string sk;
        while (f >> sk) {
          if (sk == "endstation") break;
          if (sk == "lastUpdateDay") {
            f >> ov.state.lastUpdateDay;
          } else if (sk == "lastSampleDay") {
            f >> ov.state.lastSampleDay;
          } else if (sk == "inventory") {
            for (std::size_t i = 0; i < econ::kCommodityCount; ++i) f >> ov.state.inventory[i];
          } else if (sk == "history") {
            std::size_t cid = 0;
            std::size_t count = 0;
            f >> cid >> count;
            if (cid >= econ::kCommodityCount) return false;
            auto& hist = ov.state.history[cid];
            hist.clear();
            hist.reserve(count);
            for (std::size_t j = 0; j < count; ++j) {
              econ::PricePoint p{};
              f >> p.day >> p.price;
              hist.push_back(p);
            }
          } else {
            // Unknown token: attempt to skip line
            std::string line;
            std::getline(f, line);
          }
        }

        out.stationOverrides.push_back(std::move(ov));
      }
    } else {
      // Unknown key: skip rest of line
      std::string line;
      std::getline(f, line);
    }
  }

  return true;
}

} // namespace stellar::sim
