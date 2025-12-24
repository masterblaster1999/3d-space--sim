#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/econ/Market.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/render/Camera.h"
#include "stellar/render/Gl.h"
#include "stellar/render/LineRenderer.h"
#include "stellar/render/Mesh.h"
#include "stellar/render/MeshRenderer.h"
#include "stellar/render/Texture.h"
#include "stellar/sim/Orbit.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/Universe.h"

#include <SDL.h>
#include <SDL_opengl.h>

#include <imgui.h>
#include <backends/imgui_impl_opengl3.h>
#include <backends/imgui_impl_sdl2.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace stellar;

static constexpr double kAU_KM = 149597870.7;
static constexpr double kSOLAR_RADIUS_KM = 695700.0;
static constexpr double kEARTH_RADIUS_KM = 6371.0;

// Rendering scale:
// The sim uses kilometers. For rendering, we scale down by this factor.
static constexpr double kRENDER_UNIT_KM = 1.0e6; // 1 unit = 1 million km

static void matToFloat(const math::Mat4d& m, float out[16]) {
  for (int i = 0; i < 16; ++i) out[i] = static_cast<float>(m.m[i]);
}

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

static const char* stationTypeName(econ::StationType t) {
  switch (t) {
    case econ::StationType::Outpost: return "Outpost";
    case econ::StationType::Agricultural: return "Agricultural";
    case econ::StationType::Mining: return "Mining";
    case econ::StationType::Refinery: return "Refinery";
    case econ::StationType::Industrial: return "Industrial";
    case econ::StationType::Research: return "Research";
    case econ::StationType::TradeHub: return "Trade Hub";
    case econ::StationType::Shipyard: return "Shipyard";
    default: return "?";
  }
}

static double cargoMassKg(const std::array<double, econ::kCommodityCount>& cargo) {
  double kg = 0.0;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const auto cid = (econ::CommodityId)i;
    const double units = cargo[i];
    if (units <= 0.0) continue;
    kg += units * econ::commodityDef(cid).massKg;
  }
  return kg;
}

static double clampRep(double r) {
  return std::clamp(r, -100.0, 100.0);
}

static double repNorm(double r) {
  return clampRep(r) / 100.0;
}

static double applyRepToFee(double baseFeeRate, double rep) {
  // Positive rep reduces fees, negative rep increases fees.
  const double kMaxEffect = 0.35; // +/- 35%
  const double eff = baseFeeRate * (1.0 - kMaxEffect * repNorm(rep));
  return std::clamp(eff, 0.0, 0.25);
}

static math::Vec3d toRenderU(const math::Vec3d& km) { return km * (1.0 / kRENDER_UNIT_KM); }

static math::Quatd quatFromTo(const math::Vec3d& from, const math::Vec3d& to) {
  math::Vec3d f = from.normalized();
  math::Vec3d t = to.normalized();
  const double c = std::clamp(math::dot(f, t), -1.0, 1.0);

  if (c > 0.999999) return math::Quatd::identity();

  if (c < -0.999999) {
    math::Vec3d axis = math::cross({1,0,0}, f);
    if (axis.lengthSq() < 1e-12) axis = math::cross({0,1,0}, f);
    return math::Quatd::fromAxisAngle(axis, math::kPi);
  }

  const math::Vec3d axis = math::cross(f, t);
  const double s = std::sqrt((1.0 + c) * 2.0);
  const double invs = 1.0 / s;
  return math::Quatd{ s * 0.5, axis.x * invs, axis.y * invs, axis.z * invs }.normalized();
}

static math::Vec3d stationPosKm(const sim::Station& st, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(st.orbit, timeDays);
  return posAU * kAU_KM;
}

static math::Vec3d stationVelKmS(const sim::Station& st, double timeDays) {
  // Numeric derivative over a short window.
  const double epsDays = 0.001; // 86.4s
  const math::Vec3d a = stationPosKm(st, timeDays - epsDays);
  const math::Vec3d b = stationPosKm(st, timeDays + epsDays);
  const double dt = epsDays * 2.0 * 86400.0;
  return (b - a) * (1.0 / dt);
}

static math::Quatd stationOrient(const sim::Station& st, const math::Vec3d& posKm, double timeDays) {
  // Station local +Z points outward from the slot.
  // We point it away from the star (radial outward) for intuitive docking.
  math::Vec3d outward = posKm.normalized();
  if (outward.lengthSq() < 1e-12) outward = {0,0,1};
  math::Quatd base = quatFromTo({0,0,1}, outward);

  // Add a small spin to make orientation/geometry feel alive.
  // Spin around the station's local +Z axis.
  const double spinRadPerDay = math::degToRad(35.0);
  math::Quatd spin = math::Quatd::fromAxisAngle({0,0,1}, std::fmod(timeDays * spinRadPerDay, 2.0 * math::kPi));
  return (base * spin).normalized();
}

static render::InstanceData makeInst(const math::Vec3d& posU,
                                    const math::Vec3d& scaleU,
                                    const math::Quatd& qWorld,
                                    float cr, float cg, float cb) {
  // shader expects quat as x,y,z,w
  const math::Quatd q = qWorld.normalized();
  return render::InstanceData{
    (float)posU.x, (float)posU.y, (float)posU.z,
    (float)scaleU.x, (float)scaleU.y, (float)scaleU.z,
    (float)q.x, (float)q.y, (float)q.z, (float)q.w,
    cr, cg, cb
  };
}

static render::InstanceData makeInstUniform(const math::Vec3d& posU,
                                           double s,
                                           float cr, float cg, float cb) {
  return makeInst(posU, {s,s,s}, math::Quatd::identity(), cr,cg,cb);
}

static bool projectToScreen(const math::Vec3d& worldU,
                            const math::Mat4d& view,
                            const math::Mat4d& proj,
                            int w, int h,
                            ImVec2& outPx) {
  // clip = proj * view * vec4(world, 1)
  const math::Mat4d vp = proj * view;
  const double x = worldU.x, y = worldU.y, z = worldU.z;

  const double cx = vp.m[0]*x + vp.m[4]*y + vp.m[8]*z + vp.m[12];
  const double cy = vp.m[1]*x + vp.m[5]*y + vp.m[9]*z + vp.m[13];
  const double cz = vp.m[2]*x + vp.m[6]*y + vp.m[10]*z + vp.m[14];
  const double cw = vp.m[3]*x + vp.m[7]*y + vp.m[11]*z + vp.m[15];

  if (cw <= 1e-6) return false;

  const double ndcX = cx / cw;
  const double ndcY = cy / cw;
  // const double ndcZ = cz / cw;

  // off-screen (still allow a bit of slack)
  if (ndcX < -1.2 || ndcX > 1.2 || ndcY < -1.2 || ndcY > 1.2) return false;

  outPx.x = (float)((ndcX * 0.5 + 0.5) * (double)w);
  outPx.y = (float)((-ndcY * 0.5 + 0.5) * (double)h);
  return true;
}

static bool beginStationSelectorHUD(const sim::StarSystem& sys, int& stationIndex, bool docked, sim::StationId dockedId) {
  bool changed = false;

  ImGui::Begin("Dock / Station");

  ImGui::Text("System: %s  (Star %s, planets %d, stations %d)",
              sys.stub.name.c_str(),
              starClassName(sys.stub.primaryClass),
              sys.stub.planetCount,
              sys.stub.stationCount);

  if (!sys.stations.empty()) {
    std::vector<const char*> names;
    names.reserve(sys.stations.size());
    for (const auto& st : sys.stations) names.push_back(st.name.c_str());

    int old = stationIndex;
    ImGui::Combo("Station", &stationIndex, names.data(), (int)names.size());
    changed = (old != stationIndex);

    const auto& st = sys.stations[(std::size_t)stationIndex];
    ImGui::SameLine();
    ImGui::TextDisabled("(%s, fee %.1f%%)", stationTypeName(st.type), st.feeRate * 100.0);

    if (docked) {
      ImGui::Text("Docked at: %s", (st.id == dockedId) ? st.name.c_str() : "(other station)");
    } else {
      ImGui::TextDisabled("Not docked");
    }
  } else {
    ImGui::Text("No stations in system.");
  }

  ImGui::End();
  return changed;
}

struct ToastMsg {
  std::string text;
  double ttl{3.0}; // seconds remaining
};

static void toast(std::vector<ToastMsg>& toasts, std::string msg, double ttlSec=3.0) {
  toasts.push_back({std::move(msg), ttlSec});
}

struct ClearanceState {
  bool granted{false};
  double expiresDays{0.0};
  double cooldownUntilDays{0.0};
};

enum class TargetKind : int { None=0, Station=1, Planet=2, Contact=3 };

struct Target {
  TargetKind kind{TargetKind::None};
  std::size_t index{0}; // station/planet/contact index
};

struct Contact {
  core::u64 id{0};
  std::string name;
  sim::Ship ship{};
  bool pirate{false};

  bool missionTarget{false};

  double shield{60.0};
  double hull{70.0};

  double fireCooldown{0.0}; // seconds
  bool alive{true};
};

static void applyDamage(double dmg, double& shield, double& hull) {
  if (shield > 0.0) {
    const double s = std::min(shield, dmg);
    shield -= s;
    dmg -= s;
  }
  if (dmg > 0.0) {
    hull -= dmg;
  }
}

static void emitStationGeometry(const sim::Station& st,
                                const math::Vec3d& stPosKm,
                                const math::Quatd& stQ,
                                std::vector<render::InstanceData>& outCubeInstances) {
  // Render in render-units
  const math::Vec3d posU = toRenderU(stPosKm);

  // Station scale factor: exaggerate a bit so it's readable at prototype camera distances.
  const double s = std::max(0.8, (st.radiusKm / kRENDER_UNIT_KM) * 1800.0);
  const math::Vec3d baseScaleU = {s, s, s};

  auto addPart = [&](const math::Vec3d& localPosU, const math::Vec3d& localScaleU, float r, float g, float b) {
    const math::Vec3d worldPosU = posU + stQ.rotate(localPosU);
    const math::Quatd worldQ = stQ;
    outCubeInstances.push_back(makeInst(worldPosU, localScaleU, worldQ, r,g,b));
  };

  // Central body (elongated)
  addPart({0,0,0}, {baseScaleU.x*0.9, baseScaleU.y*0.9, baseScaleU.z*1.3}, 0.65f,0.68f,0.72f);

  // Ring (8 segments)
  const double ringR = baseScaleU.x * 1.25;
  const double segLen = baseScaleU.z * 0.55;
  for (int i = 0; i < 8; ++i) {
    const double ang = (double)i / 8.0 * 2.0 * math::kPi;
    const double cx = std::cos(ang), sx = std::sin(ang);
    const math::Vec3d lp{ ringR * cx, ringR * sx, 0.0 };
    addPart(lp, {baseScaleU.x*0.22, baseScaleU.y*0.22, segLen}, 0.45f,0.55f,0.70f);
  }

  // Docking "mail slot" frame on +Z face
  const double frameZ = baseScaleU.z * 1.3;
  const double frameW = baseScaleU.x * 1.15;
  const double frameH = baseScaleU.y * 0.60;
  const double thickness = baseScaleU.x * 0.12;

  // Left/right pylons
  addPart({-frameW*0.5, 0, frameZ}, {thickness, frameH, thickness}, 0.85f,0.80f,0.55f);
  addPart({+frameW*0.5, 0, frameZ}, {thickness, frameH, thickness}, 0.85f,0.80f,0.55f);

  // Top/bottom bars
  addPart({0, +frameH*0.5, frameZ}, {frameW, thickness, thickness}, 0.85f,0.80f,0.55f);
  addPart({0, -frameH*0.5, frameZ}, {frameW, thickness, thickness}, 0.85f,0.80f,0.55f);

  // A small "light strip" above slot
  addPart({0, +frameH*0.65, frameZ}, {frameW*0.8, thickness*0.6, thickness*0.6}, 0.25f,0.85f,0.35f);
}

static bool insideStationHullExceptSlot(const sim::Station& st,
                                       const math::Vec3d& relLocalKm) {
  // Local station hull approximation: box (wx,wy,wz), with a rectangular slot tunnel cut out.
  const double wx = st.radiusKm * 0.70;
  const double wy = st.radiusKm * 0.70;
  const double wz = st.radiusKm * 1.10;

  const bool insideBox = (std::abs(relLocalKm.x) < wx) && (std::abs(relLocalKm.y) < wy) && (std::abs(relLocalKm.z) < wz);
  if (!insideBox) return false;

  // Slot tunnel cutout near +Z face (entrance at +wz)
  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  const double zEntrance = wz;
  const double zMin = zEntrance - st.slotDepthKm;

  const bool insideTunnel =
    (std::abs(relLocalKm.x) < slotHalfW) &&
    (std::abs(relLocalKm.y) < slotHalfH) &&
    (relLocalKm.z <= zEntrance) &&
    (relLocalKm.z >= zMin);

  return !insideTunnel;
}

static bool dockingSlotConditions(const sim::Station& st,
                                 const math::Vec3d& relLocalKm,
                                 const math::Vec3d& shipVelRelLocalKmS,
                                 const math::Vec3d& shipForwardLocal,
                                 const math::Vec3d& shipUpLocal,
                                 bool clearanceGranted) {
  if (!clearanceGranted) return false;

  const double wx = st.radiusKm * 0.70;
  const double wy = st.radiusKm * 0.70;
  const double wz = st.radiusKm * 1.10;
  const double zEntrance = wz;
  const double zMin = zEntrance - st.slotDepthKm;

  // Must be inside tunnel volume (i.e. have flown into the mail-slot)
  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  const bool insideTunnel =
    (std::abs(relLocalKm.x) < slotHalfW) &&
    (std::abs(relLocalKm.y) < slotHalfH) &&
    (relLocalKm.z <= zEntrance - 0.05 * st.radiusKm) &&
    (relLocalKm.z >= zMin + 0.10 * st.radiusKm);

  if (!insideTunnel) return false;

  const double relSpeed = shipVelRelLocalKmS.length();
  if (relSpeed > st.maxApproachSpeedKmS) return false;

  // Orientation: ship forward should point into the station (-Z local),
  // and ship up should be roughly +Y local (roll alignment).
  const double fwdAlign = math::dot(shipForwardLocal.normalized(), math::Vec3d{0,0,-1});
  const double upAlign = math::dot(shipUpLocal.normalized(), math::Vec3d{0,1,0});
  return (fwdAlign > 0.92 && upAlign > 0.70);
}

static double systemDistanceLy(const sim::SystemStub& a, const sim::SystemStub& b) {
  return (a.posLy - b.posLy).length();
}

// Very small/fast route plotter for early gameplay.
// Uses A* on hop count (minimizes number of jumps) with an admissible heuristic
// (ceil(remainingDistance / maxJump)).
static std::vector<sim::SystemId> plotRouteAStarHops(const std::vector<sim::SystemStub>& nodes,
                                                     sim::SystemId startId,
                                                     sim::SystemId goalId,
                                                     double maxJumpLy) {
  if (startId == 0 || goalId == 0) return {};
  if (maxJumpLy <= 0.0) return {};
  if (nodes.empty()) return {};

  std::unordered_map<sim::SystemId, std::size_t> idx;
  idx.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) idx[nodes[i].id] = i;

  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) return {};

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  std::vector<int> cameFrom(N, -1);
  std::vector<int> gScore(N, std::numeric_limits<int>::max());

  auto heuristic = [&](std::size_t i) -> int {
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    return (int)std::ceil(d / maxJumpLy);
  };

  struct QN {
    int f;
    int g;
    std::size_t i;
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const { return a.f > b.f; }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0;
  open.push({heuristic(start), 0, start});

  std::vector<char> closed(N, 0);

  while (!open.empty()) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    if (cur.i == goal) {
      std::vector<sim::SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());
      return path;
    }

    // Neighbors: any node within maxJumpLy.
    for (std::size_t j = 0; j < N; ++j) {
      if (j == cur.i) continue;
      if (closed[j]) continue;
      const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
      if (d > maxJumpLy + 1e-9) continue;

      const int tentative = gScore[cur.i] + 1;
      if (tentative < gScore[j]) {
        gScore[j] = tentative;
        cameFrom[j] = (int)cur.i;
        const int f = tentative + heuristic(j);
        open.push({f, tentative, j});
      }
    }
  }

  return {};
}

int main(int argc, char** argv) {
  (void)argc; (void)argv;

  core::setLogLevel(core::LogLevel::Info);

  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0) {
    core::log(core::LogLevel::Error, std::string("SDL_Init failed: ") + SDL_GetError());
    return 1;
  }

  // GL 3.3 core
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  SDL_Window* window = SDL_CreateWindow(
      "Stellar Forge (prototype)",
      SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
      1280, 720,
      SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

  if (!window) {
    core::log(core::LogLevel::Error, std::string("SDL_CreateWindow failed: ") + SDL_GetError());
    return 1;
  }

  SDL_GLContext glContext = SDL_GL_CreateContext(window);
  SDL_GL_MakeCurrent(window, glContext);
  SDL_GL_SetSwapInterval(1);

  if (!render::gl::load()) {
    core::log(core::LogLevel::Error, "Failed to load OpenGL functions.");
    return 1;
  }

  core::log(core::LogLevel::Info, std::string("OpenGL: ") + render::gl::glVersionString());

  glEnable(GL_DEPTH_TEST);

  // --- ImGui setup ---
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  ImGui::StyleColorsDark();

  ImGui_ImplSDL2_InitForOpenGL(window, glContext);
  ImGui_ImplOpenGL3_Init("#version 330 core");

  // --- Render assets ---
  render::Mesh sphere = render::Mesh::makeUvSphere(48, 24);
  render::Mesh cube   = render::Mesh::makeCube();

  render::Texture2D checker;
  checker.createChecker(256, 256, 16);

  render::MeshRenderer meshRenderer;
  std::string err;
  if (!meshRenderer.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }
  meshRenderer.setTexture(&checker);

  render::LineRenderer lineRenderer;
  if (!lineRenderer.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }

  // --- Universe / sim state ---
  core::u64 seed = 1337;
  sim::Universe universe(seed);

  sim::SystemStub currentStub{};
  if (auto s = universe.findClosestSystem({0,0,0}, 80.0)) {
    currentStub = *s;
  } else {
    // Should be rare; fall back
    currentStub = sim::SystemStub{};
    currentStub.id = 0;
    currentStub.seed = seed;
    currentStub.name = "Origin";
    currentStub.posLy = {0,0,0};
    currentStub.planetCount = 6;
    currentStub.stationCount = 1;
  }

  const sim::StarSystem* currentSystem = &universe.getSystem(currentStub.id, &currentStub);

  sim::Ship ship;
  const double kPlayerBaseLinAccelKmS2 = 0.08;
  const double kPlayerBaseAngAccelRadS2 = 1.2;
  ship.setMaxLinearAccelKmS2(kPlayerBaseLinAccelKmS2);
  ship.setMaxAngularAccelRadS2(kPlayerBaseAngAccelRadS2);

  // Player combat state
  double playerShield = 100.0;
  double playerHull = 100.0;
  double playerLaserCooldown = 0.0; // seconds

  // Time
  double timeDays = 0.0;
  double timeScale = 60.0; // simulated seconds per real second
  bool paused = false;

  // Economy
  double credits = 2500.0;
  std::array<double, econ::kCommodityCount> cargo{};
  double cargoCapacityKg = 420.0;
  int selectedStationIndex = 0;

  // Ship meta / progression
  double fuel = 45.0;
  double fuelMax = 45.0;
  double fsdRangeLy = 18.0;
  double fsdReadyDay = 0.0;

  // Missions + reputation
  core::u64 nextMissionId = 1;
  std::vector<sim::Mission> missions;
  std::unordered_map<core::u32, double> repByFaction;

  // Cached mission board offers (regenerated when docking / day changes)
  std::vector<sim::Mission> missionOffers;
  sim::StationId missionOffersStationId = 0;
  int missionOffersDayStamp = -1;

  // Docking state
  bool docked = false;
  sim::StationId dockedStationId = 0;
  std::unordered_map<sim::StationId, ClearanceState> clearances;

  // Contacts (pirates etc.)
  std::vector<Contact> contacts;
  core::SplitMix64 rng(seed ^ 0xC0FFEEu);
  double nextPirateSpawnDays = 0.01; // soon after start

  // Beams (for laser visuals)
  struct Beam { math::Vec3d aU, bU; float r,g,b; double ttl; };
  std::vector<Beam> beams;

  // Save/load
  const std::string savePath = "savegame.txt";

  // UI state
  bool showGalaxy = true;
  bool showShip = true;
  bool showEconomy = true;
  bool showMissions = true;
  bool showContacts = true;

  // Flight assistance
  bool autopilot = false;

  // Mouse flight (optional): mouse-driven pitch/yaw while the cursor is captured.
  // Toggle with Right Mouse or M. Press Esc to release capture.
  bool mouseFlight = false;
  bool invertMouseY = false;
  float mouseSensitivity = 0.010f; // normalized torque input per pixel (tune in UI)

  // Supercruise-lite
  bool supercruise = false;
  bool supercruiseAssist = true;

  // FSD / hyperspace (system-to-system)
  enum class FsdState { Idle, Charging, Travelling };
  FsdState fsdState = FsdState::Idle;
  sim::SystemId fsdTargetSystem{0};
  double fsdChargeRemainingSec = 0.0;
  double fsdTravelRemainingSec = 0.0;
  double fsdFuelCost = 0.0;
  double fsdJumpDistanceLy = 0.0;

  // Galaxy navigation / route plotting
  sim::SystemId galaxySelectedSystemId = 0;
  std::vector<sim::SystemId> navRoute;
  std::size_t navRouteHop = 0;
  bool navAutoRun = false;

  // Bounty scan interaction
  bool scanning = false;
  double scanProgressSec = 0.0;

  // Target
  Target target{};

  std::vector<ToastMsg> toasts;

  auto setMouseFlight = [&](bool enabled) {
    if (mouseFlight == enabled) return;
    mouseFlight = enabled;

    SDL_SetRelativeMouseMode(enabled ? SDL_TRUE : SDL_FALSE);
    SDL_ShowCursor(enabled ? SDL_DISABLE : SDL_ENABLE);

    // Flush any pending relative motion so we don't get a "jump" on toggle.
    int dx = 0, dy = 0;
    (void)SDL_GetRelativeMouseState(&dx, &dy);

    if (enabled) {
      toast(toasts, "Mouse flight: ON (RMB/M to toggle, Esc to release)", 2.5);
    } else {
      toast(toasts, "Mouse flight: OFF", 1.5);
    }
  };

  auto respawnNearStation = [&](const sim::StarSystem& sys, std::size_t stationIdx) {
    if (sys.stations.empty()) {
      ship.setPositionKm({0,0,-8000.0});
      ship.setVelocityKmS({0,0,0});
      ship.setAngularVelocityRadS({0,0,0});
      ship.setOrientation(math::Quatd::identity());
      return;
    }
    stationIdx = std::min(stationIdx, sys.stations.size() - 1);
    const auto& st = sys.stations[stationIdx];

    const math::Vec3d stPos = stationPosKm(st, timeDays);
    const math::Quatd stQ = stationOrient(st, stPos, timeDays);
    const math::Vec3d axis = stQ.rotate({0,0,1}); // outward
    const double startDistKm = st.radiusKm * 14.0;
    ship.setPositionKm(stPos + axis * startDistKm);
    ship.setVelocityKmS(stationVelKmS(st, timeDays));
    ship.setAngularVelocityRadS({0,0,0});
    // face toward the slot (into station)
    ship.setOrientation(quatFromTo({0,0,1}, -axis));
  };

  // Spawn near first station for immediate gameplay.
  respawnNearStation(*currentSystem, 0);

  galaxySelectedSystemId = currentSystem->stub.id;

  auto findFaction = [&](core::u32 factionId) -> const sim::Faction* {
    for (const auto& f : universe.factions()) {
      if (f.id == factionId) return &f;
    }
    return nullptr;
  };

  auto factionName = [&](core::u32 factionId) -> std::string {
    if (factionId == 0) return "Independent";
    if (auto f = findFaction(factionId)) return f->name;
    return "Faction " + std::to_string(factionId);
  };

  auto getRep = [&](core::u32 factionId) -> double {
    auto it = repByFaction.find(factionId);
    return it == repByFaction.end() ? 0.0 : it->second;
  };

  auto addRep = [&](core::u32 factionId, double delta) {
    if (factionId == 0) return;
    repByFaction[factionId] = clampRep(getRep(factionId) + delta);
  };

  auto effectiveFeeRate = [&](const sim::Station& st) -> double {
    return applyRepToFee(st.feeRate, getRep(st.factionId));
  };

  // FSD / jump parameters
  const double kFsdFuelBase = 2.0;
  const double kFsdFuelPerLy = 0.5;
  const double kFsdChargeSec = 4.0;

  auto fsdBaseRangeLyNow = [&]() -> double {
    const double cap = std::max(1.0, cargoCapacityKg);
    const double load = std::clamp(cargoMassKg(cargo) / cap, 0.0, 1.0);
    // Cargo load reduces effective range a bit (keeps hauling interesting).
    return std::max(0.0, fsdRangeLy * (1.0 - 0.25 * load));
  };

  auto fsdFuelLimitedRangeLy = [&]() -> double {
    if (fuel <= kFsdFuelBase) return 0.0;
    return std::max(0.0, (fuel - kFsdFuelBase) / kFsdFuelPerLy);
  };

  auto fsdCurrentRangeLy = [&]() -> double {
    return std::min(fsdBaseRangeLyNow(), fsdFuelLimitedRangeLy());
  };

  auto fsdFuelCostFor = [&](double distanceLy) -> double {
    return kFsdFuelBase + distanceLy * kFsdFuelPerLy;
  };

  auto isMassLocked = [&]() -> bool {
    if (!currentSystem) return false;
    const math::Vec3d shipPos = ship.positionKm();

    // Stations (treat as heavy bodies / traffic control).
    for (const auto& st : currentSystem->stations) {
      const math::Vec3d stPos = stationPosKm(st, timeDays);
      const double distKm = (shipPos - stPos).length();
      const double lockKm = st.radiusKm * 12.0;
      if (distKm < lockKm) return true;
    }

    // Planets (simple gravity-well proxy).
    for (const auto& p : currentSystem->planets) {
      const math::Vec3d pPos = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
      const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const double distKm = (shipPos - pPos).length();
      const double lockKm = std::max(5000.0, rKm * 20.0);
      if (distKm < lockKm) return true;
    }

    return false;
  };

  auto startFsdJumpTo = [&](sim::SystemId destId) {
    if (destId == 0 || destId == currentSystem->stub.id) {
      toast(toasts, "No destination selected.", 1.8);
      return;
    }
    if (docked) {
      toast(toasts, "Can't jump while docked.", 2.0);
      return;
    }
    if (supercruise) {
      toast(toasts, "Disengage supercruise before jumping.", 2.0);
      return;
    }
    if (fsdState != FsdState::Idle) {
      toast(toasts, "FSD already busy.", 2.0);
      return;
    }
    if (timeDays < fsdReadyDay) {
      toast(toasts, "FSD cooling down...", 1.8);
      return;
    }
    if (isMassLocked()) {
      toast(toasts, "Mass-locked: move away from bodies/stations.", 2.2);
      return;
    }

    const sim::StarSystem& destSys = universe.getSystem(destId);
    const double distLy = (destSys.stub.posLy - currentSystem->stub.posLy).length();

    const double rangeLy = fsdBaseRangeLyNow();
    if (distLy > rangeLy + 1e-9) {
      toast(toasts, "Out of jump range (plot a multi-jump route).", 2.5);
      return;
    }

    const double fuelCost = fsdFuelCostFor(distLy);
    if (fuel < fuelCost) {
      toast(toasts, "Not enough fuel for this jump.", 2.5);
      return;
    }

    // Consume fuel on charge complete (so you can still cancel cleanly).
    fsdTargetSystem = destId;
    fsdJumpDistanceLy = distLy;
    fsdFuelCost = fuelCost;
    fsdChargeRemainingSec = kFsdChargeSec;
    fsdTravelRemainingSec = 0.0;
    fsdState = FsdState::Charging;

    autopilot = false;
    scanning = false;
    scanProgressSec = 0.0;
    beams.clear();

    toast(toasts, "FSD charging...", 2.0);
  };

  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

  SDL_SetRelativeMouseMode(SDL_FALSE);
  SDL_ShowCursor(SDL_ENABLE);

  while (running) {
    // Timing
    auto now = std::chrono::high_resolution_clock::now();
    const double dtReal = std::chrono::duration<double>(now - last).count();
    last = now;

    // Per-frame mouse deltas (used for mouse flight).
    double mouseDx = 0.0;
    double mouseDy = 0.0;

    // Events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      ImGui_ImplSDL2_ProcessEvent(&event);

      if (event.type == SDL_QUIT) running = false;
      if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE) running = false;

      if (event.type == SDL_KEYDOWN && !event.key.repeat) {
        if (event.key.keysym.sym == SDLK_ESCAPE) {
          // Escape either releases mouse flight (if active) or quits.
          if (mouseFlight) {
            setMouseFlight(false);
          } else {
            running = false;
          }
        }

        if (event.key.keysym.sym == SDLK_m) {
          // Toggle mouse flight.
          if (!io.WantCaptureKeyboard) {
            setMouseFlight(!mouseFlight);
          }
        }

        if (event.key.keysym.sym == SDLK_F5) {
          sim::SaveGame s{};
          s.seed = universe.seed();
          s.timeDays = timeDays;
          s.currentSystem = currentSystem->stub.id;
          s.dockedStation = docked ? dockedStationId : 0;

          s.shipPosKm = ship.positionKm();
          s.shipVelKmS = ship.velocityKmS();
          s.shipOrient = ship.orientation();
          s.shipAngVelRadS = ship.angularVelocityRadS();

          s.credits = credits;
          s.cargo = cargo;
          s.cargoCapacityKg = cargoCapacityKg;
          s.fuel = fuel;
          s.fuelMax = fuelMax;
          s.fsdRangeLy = fsdRangeLy;
          s.hull = std::clamp(playerHull / 100.0, 0.0, 1.0);
          s.fsdReadyDay = fsdReadyDay;

          s.nextMissionId = nextMissionId;
          s.missions = missions;

          s.reputation.clear();
          s.reputation.reserve(repByFaction.size());
          for (const auto& kv : repByFaction) {
            sim::FactionReputation r{};
            r.factionId = kv.first;
            r.rep = kv.second;
            s.reputation.push_back(r);
          }
          s.stationOverrides = universe.exportStationOverrides();

          if (sim::saveToFile(s, savePath)) {
            toast(toasts, "Saved to " + savePath, 2.5);
          }
        }

        if (event.key.keysym.sym == SDLK_F9) {
          sim::SaveGame s{};
          if (sim::loadFromFile(savePath, s)) {
            universe = sim::Universe(s.seed);
            universe.importStationOverrides(s.stationOverrides);

            timeDays = s.timeDays;

            const sim::StarSystem& sys = universe.getSystem(s.currentSystem);
            currentStub = sys.stub;
            currentSystem = &sys;

            ship.setPositionKm(s.shipPosKm);
            ship.setVelocityKmS(s.shipVelKmS);
            ship.setOrientation(s.shipOrient);
            ship.setAngularVelocityRadS(s.shipAngVelRadS);

            credits = s.credits;
            cargo = s.cargo;

            cargoCapacityKg = s.cargoCapacityKg;
            fuel = s.fuel;
            fuelMax = s.fuelMax;
            fsdRangeLy = s.fsdRangeLy;
            fsdReadyDay = s.fsdReadyDay;
            playerHull = std::clamp(s.hull, 0.0, 1.0) * 100.0;

            nextMissionId = s.nextMissionId;
            missions = s.missions;

            repByFaction.clear();
            for (const auto& r : s.reputation) repByFaction[r.factionId] = r.rep;

            docked = (s.dockedStation != 0);
            dockedStationId = s.dockedStation;
            selectedStationIndex = 0;
            if (docked) {
              for (std::size_t i = 0; i < sys.stations.size(); ++i) {
                if (sys.stations[i].id == s.dockedStation) selectedStationIndex = (int)i;
              }
            }

            // clear transient runtime things
            contacts.clear();
            beams.clear();
            autopilot = false;
            supercruise = false;
            fsdState = FsdState::Idle;
            fsdTargetSystem = 0;
            navRoute.clear();
            navRouteHop = 0;
            navAutoRun = false;
            scanning = false;
            scanProgressSec = 0.0;
            clearances.clear();
            target = Target{};

            galaxySelectedSystemId = currentSystem->stub.id;

            toast(toasts, "Loaded " + savePath, 2.5);
          }
        }

        if (event.key.keysym.sym == SDLK_TAB) showGalaxy = !showGalaxy;
        if (event.key.keysym.sym == SDLK_F1) showShip = !showShip;
        if (event.key.keysym.sym == SDLK_F2) showEconomy = !showEconomy;
        if (event.key.keysym.sym == SDLK_F4) showMissions = !showMissions;
        if (event.key.keysym.sym == SDLK_F3) showContacts = !showContacts;

        if (event.key.keysym.sym == SDLK_SPACE) paused = !paused;

        if (event.key.keysym.sym == SDLK_p) autopilot = !autopilot;

        if (event.key.keysym.sym == SDLK_h) {
          // Supercruise-lite toggle (in-system travel assist).
          if (docked) {
            toast(toasts, "Can't engage supercruise while docked.", 2.0);
          } else if (supercruise) {
            supercruise = false;
            toast(toasts, "Supercruise disengaged.", 1.5);
          } else {
            if (target.kind == TargetKind::None) {
              toast(toasts, "Set a target (T/B/N) to use supercruise.", 2.5);
            } else {
              supercruise = true;
              autopilot = false;
              beams.clear();
              toast(toasts, "Supercruise engaged.", 1.5);
            }
          }
        }

        if (event.key.keysym.sym == SDLK_k) {
          // Bounty scan action
          if (scanning) {
            scanning = false;
            scanProgressSec = 0.0;
            toast(toasts, "Scan cancelled.", 1.5);
          } else {
            scanning = true;
            scanProgressSec = 0.0;
            toast(toasts, "Scanning...", 1.5);
          }
        }

        if (event.key.keysym.sym == SDLK_j) {
          // Jump to next hop in the plotted route, otherwise jump to selected system.
          if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
            startFsdJumpTo(navRoute[navRouteHop + 1]);
          } else {
            startFsdJumpTo(galaxySelectedSystemId);
          }
        }
        if (event.key.keysym.sym == SDLK_t) {
          // cycle station targets
          if (!currentSystem->stations.empty()) {
            if (target.kind != TargetKind::Station) {
              target.kind = TargetKind::Station;
              target.index = 0;
            } else {
              target.index = (target.index + 1) % currentSystem->stations.size();
            }
          }
        }

        if (event.key.keysym.sym == SDLK_b) {
          // Cycle planet targets (bodies) for supercruise.
          if (!currentSystem->planets.empty()) {
            if (target.kind != TargetKind::Planet) {
              target.kind = TargetKind::Planet;
              target.index = 0;
            } else {
              target.index = (target.index + 1) % currentSystem->planets.size();
            }
          }
        }

        if (event.key.keysym.sym == SDLK_n) {
          // Cycle contacts (prefers living pirates).
          if (!contacts.empty()) {
            std::size_t start = (target.kind == TargetKind::Contact) ? (target.index + 1) : 0;
            bool found = false;

            // First pass: living pirates.
            for (std::size_t i = 0; i < contacts.size(); ++i) {
              std::size_t idx = (start + i) % contacts.size();
              if (contacts[idx].alive && contacts[idx].pirate) {
                target.kind = TargetKind::Contact;
                target.index = idx;
                found = true;
                break;
              }
            }

            // Second pass: any living contact.
            if (!found) {
              for (std::size_t i = 0; i < contacts.size(); ++i) {
                std::size_t idx = (start + i) % contacts.size();
                if (contacts[idx].alive) {
                  target.kind = TargetKind::Contact;
                  target.index = idx;
                  found = true;
                  break;
                }
              }
            }

            if (!found) toast(toasts, "No living contacts.", 2.0);
          }
        }
        if (event.key.keysym.sym == SDLK_y) target = Target{};

        if (event.key.keysym.sym == SDLK_l) {
          // Request docking clearance from targeted station
          if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
            const auto& st = currentSystem->stations[target.index];
            const math::Vec3d stPos = stationPosKm(st, timeDays);
            const double dist = (ship.positionKm() - stPos).length();

            auto& cs = clearances[st.id];

            if (dist > st.commsRangeKm) {
              toast(toasts, "Out of comms range for clearance.", 2.5);
            } else if (timeDays < cs.cooldownUntilDays) {
              toast(toasts, "Clearance channel busy. Try again soon.", 2.5);
            } else {
              // Simple logic: usually granted; sometimes denied (simulate traffic / capacity).
              const double pGrant = 0.82;
              cs.granted = rng.chance(pGrant);
              if (cs.granted) {
                cs.expiresDays = timeDays + (12.0 * 60.0) / 86400.0; // 12 minutes
                toast(toasts, "Docking clearance GRANTED.", 3.0);
              } else {
                cs.cooldownUntilDays = timeDays + (90.0) / 86400.0; // 90 seconds
                toast(toasts, "Docking clearance DENIED. (traffic)", 3.0);
              }
            }
          } else {
            toast(toasts, "No station targeted for clearance.", 2.0);
          }
        }

        if (event.key.keysym.sym == SDLK_g) {
          if (docked) {
            // Undock: place ship just outside the slot
            const auto it = std::find_if(currentSystem->stations.begin(), currentSystem->stations.end(),
                                         [&](const sim::Station& s){ return s.id == dockedStationId; });
            if (it != currentSystem->stations.end()) {
              const auto& st = *it;
              const math::Vec3d stPos = stationPosKm(st, timeDays);
              const math::Quatd stQ = stationOrient(st, stPos, timeDays);
              const math::Vec3d axis = stQ.rotate({0,0,1});
              ship.setPositionKm(stPos + axis * (st.radiusKm * 1.8));
              ship.setVelocityKmS(stationVelKmS(st, timeDays));
              ship.setAngularVelocityRadS({0,0,0});
              ship.setOrientation(quatFromTo({0,0,1}, -axis));

              docked = false;
              dockedStationId = 0;
              toast(toasts, "Undocked.", 2.0);
            } else {
              docked = false;
              dockedStationId = 0;
            }
          } else {
            // Dock attempt: must be inside slot tunnel, aligned, and have clearance.
            if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
              const auto& st = currentSystem->stations[target.index];
              auto& cs = clearances[st.id];
              const bool clearanceValid = cs.granted && (timeDays <= cs.expiresDays);

              const math::Vec3d stPos = stationPosKm(st, timeDays);
              const math::Quatd stQ = stationOrient(st, stPos, timeDays);
              const math::Vec3d stV = stationVelKmS(st, timeDays);

              const math::Vec3d relWorldKm = ship.positionKm() - stPos;
              const math::Vec3d relLocalKm = stQ.conjugate().rotate(relWorldKm);
              const math::Vec3d relVelWorld = ship.velocityKmS() - stV;
              const math::Vec3d relVelLocal = stQ.conjugate().rotate(relVelWorld);

              const math::Vec3d fwdLocal = stQ.conjugate().rotate(ship.forward());
              const math::Vec3d upLocal = stQ.conjugate().rotate(ship.up());

              if (!clearanceValid) {
                toast(toasts, "No valid clearance. Press L to request.", 2.5);
              } else if (!dockingSlotConditions(st, relLocalKm, relVelLocal, fwdLocal, upLocal, clearanceValid)) {
                toast(toasts, "Docking failed: align and enter the slot under speed limit.", 2.8);
              } else {
                // Docked: lock ship to a point inside hangar.
                docked = true;
                dockedStationId = st.id;
                selectedStationIndex = (int)target.index;

                ship.setVelocityKmS(stV);
                ship.setAngularVelocityRadS({0,0,0});

                toast(toasts, "Docked at " + st.name, 2.5);
              }
            } else {
              toast(toasts, "Target a station (T) before docking.", 2.5);
            }
          }
        }
      }

      if (event.type == SDL_MOUSEMOTION) {
        if (mouseFlight && !io.WantCaptureMouse) {
          mouseDx += (double)event.motion.xrel;
          mouseDy += (double)event.motion.yrel;
        }
      }

      if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_RIGHT) {
        if (!io.WantCaptureMouse) {
          setMouseFlight(!mouseFlight);
        }
      }

      if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT) {
        if (!io.WantCaptureMouse) {
          // Fire laser
          if (playerLaserCooldown <= 0.0 && !docked && fsdState == FsdState::Idle && !supercruise) {
            playerLaserCooldown = 0.18; // rate of fire
            const double rangeKm = 120000.0;
            const double dmg = 18.0;

            const math::Vec3d aKm = ship.positionKm();
            const math::Vec3d dir = ship.forward().normalized();

            // Find best hit (simple cone + nearest along ray).
            int bestIdx = -1;
            double bestT = rangeKm;

            for (int i = 0; i < (int)contacts.size(); ++i) {
              auto& c = contacts[(std::size_t)i];
              if (!c.alive) continue;

              const math::Vec3d to = c.ship.positionKm() - aKm;
              const double dist = to.length();
              if (dist > rangeKm) continue;

              const math::Vec3d toN = to / std::max(1e-9, dist);
              const double aim = math::dot(dir, toN);
              if (aim < 0.995) continue; // narrow cone

              // approximate hit at closest along ray
              const double t = dist * aim;
              if (t < bestT) { bestT = t; bestIdx = i; }
            }

            const math::Vec3d bKm = aKm + dir * bestT;
            beams.push_back({toRenderU(aKm), toRenderU(bKm), 1.0f, 0.25f, 0.25f, 0.10});

            if (bestIdx >= 0) {
              auto& hit = contacts[(std::size_t)bestIdx];
              applyDamage(dmg, hit.shield, hit.hull);
              if (hit.hull <= 0.0) {
                const core::u64 deadId = hit.id;
                hit.alive = false;
                toast(toasts, "Target destroyed. +500 cr", 2.5);
                credits += 500.0;

                // Bounty kill missions
                for (auto& m : missions) {
                  if (m.completed || m.failed) continue;
                  if (m.type == sim::MissionType::BountyKill && m.targetNpcId == deadId && m.toSystem == currentSystem->stub.id) {
                    m.completed = true;
                    credits += m.reward;
                    addRep(m.factionId, +2.0);
                    toast(toasts, "Mission complete: bounty target eliminated. +" + std::to_string((int)m.reward) + " cr", 3.0);
                  }
                }
              }
            }
          }
        }
      }
    }

    // Input (6DOF)
    sim::ShipInput input{};
    const Uint8* keys = SDL_GetKeyboardState(nullptr);

    const bool captureKeys = io.WantCaptureKeyboard;
    if (!captureKeys && !docked) {
      input.thrustLocal.z += (keys[SDL_SCANCODE_W] ? 1.0 : 0.0);
      input.thrustLocal.z -= (keys[SDL_SCANCODE_S] ? 1.0 : 0.0);

      input.thrustLocal.x += (keys[SDL_SCANCODE_D] ? 1.0 : 0.0);
      input.thrustLocal.x -= (keys[SDL_SCANCODE_A] ? 1.0 : 0.0);

      input.thrustLocal.y += (keys[SDL_SCANCODE_R] ? 1.0 : 0.0);
      input.thrustLocal.y -= (keys[SDL_SCANCODE_F] ? 1.0 : 0.0);

      input.torqueLocal.x += (keys[SDL_SCANCODE_UP] ? 1.0 : 0.0);
      input.torqueLocal.x -= (keys[SDL_SCANCODE_DOWN] ? 1.0 : 0.0);

      input.torqueLocal.y += (keys[SDL_SCANCODE_RIGHT] ? 1.0 : 0.0);
      input.torqueLocal.y -= (keys[SDL_SCANCODE_LEFT] ? 1.0 : 0.0);

      input.torqueLocal.z += (keys[SDL_SCANCODE_E] ? 1.0 : 0.0);
      input.torqueLocal.z -= (keys[SDL_SCANCODE_Q] ? 1.0 : 0.0);

      // Mouse flight (captured mouse) feeds pitch/yaw.
      if (mouseFlight) {
        const double invY = invertMouseY ? -1.0 : 1.0;
        input.torqueLocal.y += (double)mouseDx * (double)mouseSensitivity;
        input.torqueLocal.x += (-(double)mouseDy) * (double)mouseSensitivity * invY;
      }

      input.boost = keys[SDL_SCANCODE_LSHIFT] != 0;
      input.brake = keys[SDL_SCANCODE_X] != 0;

      static bool dampers = true;
      if (keys[SDL_SCANCODE_Z]) dampers = true;
      if (keys[SDL_SCANCODE_C]) dampers = false;
      input.dampers = dampers;
    }

    // Autopilot: station approach assist (keeps you aligned to slot and under approach speed).
    if (autopilot && !docked && target.kind == TargetKind::Station && target.index < currentSystem->stations.size() && !captureKeys) {
      const auto& st = currentSystem->stations[target.index];

      const math::Vec3d stPos = stationPosKm(st, timeDays);
      const math::Quatd stQ = stationOrient(st, stPos, timeDays);
      const math::Vec3d stV = stationVelKmS(st, timeDays);

      const math::Vec3d axisOut = stQ.rotate({0,0,1});
      const math::Vec3d desiredPoint = stPos + axisOut * (st.radiusKm * 2.2); // just outside slot

      const math::Vec3d rel = desiredPoint - ship.positionKm();
      const double dist = rel.length();

      // Desired velocity: approach slowly as we get closer.
      const double maxV = st.maxApproachSpeedKmS * 0.9;
      const double vMag = std::min(maxV, 0.004 * dist); // km/s
      const math::Vec3d desiredVel = stV + (dist > 1e-6 ? rel.normalized() * vMag : math::Vec3d{0,0,0});

      const math::Vec3d dv = desiredVel - ship.velocityKmS();

      // Convert desired acceleration direction to local thrust.
      const math::Vec3d accelWorldDir = dv.normalized();
      const math::Vec3d accelLocal = ship.orientation().conjugate().rotate(accelWorldDir);

      input.thrustLocal = accelLocal;
      input.brake = false;
      input.boost = false;
      input.dampers = true;

      // Orientation control: face into the slot (-axisOut)
      const math::Vec3d desiredFwdWorld = (-axisOut).normalized();
      const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(desiredFwdWorld);

      // Yaw/pitch to bring +Z toward desired direction.
      const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
      const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

      input.torqueLocal.x = std::clamp(pitchErr * 2.0, -1.0, 1.0);
      input.torqueLocal.y = std::clamp(yawErr * 2.0, -1.0, 1.0);
      input.torqueLocal.z = 0.0;

      // Auto-request clearance when close enough.
      auto& cs = clearances[st.id];
      const double shipDist = (ship.positionKm() - stPos).length();
      if (shipDist < st.commsRangeKm && !(cs.granted && timeDays <= cs.expiresDays) && timeDays >= cs.cooldownUntilDays) {
        cs.granted = true;
        cs.expiresDays = timeDays + (12.0 * 60.0) / 86400.0;
        toast(toasts, "Autopilot: requested clearance (auto-granted).", 2.0);
      }
    }

    // Cargo affects normal handling: a fully-loaded ship accelerates and turns a bit worse.
    const double cargoLoad = std::clamp(cargoMassKg(cargo) / std::max(1.0, cargoCapacityKg), 0.0, 1.0);
    const double linHandlingMult = (1.0 - 0.25 * cargoLoad);
    const double angHandlingMult = (1.0 - 0.20 * cargoLoad);

    // Supercruise-lite: fast in-system travel assist to current target.
    // Inspired by the "supercruise" style travel in games like Elite Dangerous.
    if (supercruise && !docked && !captureKeys && fsdState == FsdState::Idle) {
      // Simple interdiction: don't allow supercruise if a pirate is nearby.
      bool interdicted = false;
      for (const auto& c : contacts) {
        if (!c.alive || !c.pirate) continue;
        const double d = (c.ship.positionKm() - ship.positionKm()).length();
        if (d < 120000.0) { interdicted = true; break; }
      }

      if (interdicted) {
        supercruise = false;
        // Restore normal handling immediately (we may have been in supercruise last frame).
        ship.setMaxLinearAccelKmS2(kPlayerBaseLinAccelKmS2 * linHandlingMult);
        ship.setMaxAngularAccelRadS2(kPlayerBaseAngAccelRadS2 * angHandlingMult);
        toast(toasts, "Supercruise blocked: hostile contact nearby.", 2.0);
      } else {
        std::optional<math::Vec3d> destPosKm;
        math::Vec3d destVelKmS{0,0,0};
        double dropKm = 60000.0;

        if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
          const auto& st = currentSystem->stations[target.index];
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          const math::Quatd stQ = stationOrient(st, stPos, timeDays);
          const math::Vec3d axisOut = stQ.rotate({0,0,1});
          destPosKm = stPos + axisOut * (st.radiusKm * 2.2);
          destVelKmS = stationVelKmS(st, timeDays);
          dropKm = std::max(45000.0, st.commsRangeKm * 0.35);
        } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
          const auto& p = currentSystem->planets[target.index];
          const math::Vec3d pPos = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
          destPosKm = pPos;
          const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
          dropKm = std::max(80000.0, rKm * 30.0);
        }

        if (!destPosKm) {
          supercruise = false;
          ship.setMaxLinearAccelKmS2(kPlayerBaseLinAccelKmS2 * linHandlingMult);
          ship.setMaxAngularAccelRadS2(kPlayerBaseAngAccelRadS2 * angHandlingMult);
          toast(toasts, "Supercruise requires a station/planet target.", 2.0);
        } else {
          const math::Vec3d rel = *destPosKm - ship.positionKm();
          const double dist = rel.length();
          const math::Vec3d dir = (dist > 1e-6) ? (rel / dist) : math::Vec3d{0,0,1};

          if (dist < dropKm) {
            supercruise = false;
            // Drop to a safe approach speed relative to destination.
            ship.setVelocityKmS(destVelKmS + dir * 0.12);
            ship.setMaxLinearAccelKmS2(kPlayerBaseLinAccelKmS2 * linHandlingMult);
            ship.setMaxAngularAccelRadS2(kPlayerBaseAngAccelRadS2 * angHandlingMult);
            toast(toasts, "Dropped from supercruise.", 1.6);
          } else {
            // Crank the ship up for supercruise, but keep controls simple.
            ship.setMaxLinearAccelKmS2(6.0);
            ship.setMaxAngularAccelRadS2(1.0);

            const double desiredSpeed = std::clamp(dist * 0.0008, 40.0, supercruiseMaxSpeedKmS);
            const math::Vec3d desiredVel = destVelKmS + dir * desiredSpeed;
            const math::Vec3d dv = desiredVel - ship.velocityKmS();
            const math::Vec3d accelWorldDir = (dv.lengthSq() > 1e-9) ? dv.normalized() : math::Vec3d{0,0,0};
            input.thrustLocal = ship.orientation().conjugate().rotate(accelWorldDir);
            input.dampers = true;
            input.brake = false;
            input.boost = true;

            // Face the destination direction.
            const math::Vec3d desiredFwdWorld = dir;
            const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(desiredFwdWorld);
            const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
            const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

            input.torqueLocal.x = std::clamp(pitchErr * 2.0, -1.0, 1.0);
            input.torqueLocal.y = std::clamp(yawErr * 2.0, -1.0, 1.0);
            input.torqueLocal.z = 0.0;
          }
        }
      }
    } else {
      // Restore normal handling when not in supercruise.
      ship.setMaxLinearAccelKmS2(kPlayerBaseLinAccelKmS2 * linHandlingMult);
      ship.setMaxAngularAccelRadS2(kPlayerBaseAngAccelRadS2 * angHandlingMult);
    }

    // Sim step
    const double dtSim = dtReal * timeScale;
    if (!paused) {
      // --- FSD (system-to-system travel) ---
      bool fsdJustArrived = false;
      if (fsdState == FsdState::Charging) {
        fsdChargeRemainingSec = std::max(0.0, fsdChargeRemainingSec - dtReal);
        if (fsdChargeRemainingSec <= 0.0) {
          // Consume fuel at the moment we enter hyperspace.
          fuel = std::clamp(fuel - fsdFuelCost, 0.0, fuelMax);
          fsdState = FsdState::Travelling;
          fsdTravelRemainingSec = 1.25 + fsdJumpDistanceLy * 0.08;
          toast(toasts, "FSD: entering hyperspace...", 1.8);
        }
      } else if (fsdState == FsdState::Travelling) {
        fsdTravelRemainingSec = std::max(0.0, fsdTravelRemainingSec - dtReal);
        if (fsdTravelRemainingSec <= 0.0) {
          const auto& nextSystem = universe.getSystem(fsdTargetSystem);
          currentStub = nextSystem.stub;
          currentSystem = &nextSystem;

          // Clear transient state.
          docked = false;
          dockedStationId = 0;
          selectedStationIndex = 0;
          clearances.clear();
          contacts.clear();
          beams.clear();
          scanning = false;
          scanProgressSec = 0.0;
          target = Target{};

          // Spawn near the first station.
          respawnNearStation(*currentSystem, 0);
          galaxySelectedSystemId = currentSystem->stub.id;

          // Cooldown (sim time)
          fsdReadyDay = timeDays + (25.0 / 86400.0);

          // Advance route cursor if this jump matches the plotted route.
          if (!navRoute.empty() && navRouteHop + 1 < navRoute.size() && navRoute[navRouteHop + 1] == currentSystem->stub.id) {
            navRouteHop++;
            if (navRouteHop + 1 >= navRoute.size()) {
              navAutoRun = false;
              toast(toasts, "Route complete.", 2.0);
            }
          }

          toast(toasts, std::string("Arrived in ") + currentSystem->stub.name + ".", 2.0);

          fsdState = FsdState::Idle;
          fsdTargetSystem = 0;
          fsdJumpDistanceLy = 0.0;
          fsdFuelCost = 0.0;
          fsdJustArrived = true;
        }
      }

      // Auto-run route (hands-free multi-jump). We only auto-trigger when safe.
      if (navAutoRun && fsdState == FsdState::Idle && timeDays >= fsdReadyDay && !docked && !supercruise) {
        if (!navRoute.empty() && navRouteHop + 1 < navRoute.size() && !isMassLocked()) {
          startFsdJumpTo(navRoute[navRouteHop + 1]);
        }
      }

      // --- Ship physics ---
      if (!docked) {
        if (fsdState != FsdState::Idle || fsdJustArrived) {
          sim::ShipInput hold{};
          hold.dampers = true;
          hold.brake = true;
          ship.step(dtSim, hold);
        } else {
          ship.step(dtSim, input);
        }
      } else {
        // While docked, keep ship attached to station.
        auto it = std::find_if(currentSystem->stations.begin(), currentSystem->stations.end(),
                               [&](const sim::Station& s){ return s.id == dockedStationId; });
        if (it != currentSystem->stations.end()) {
          const auto& st = *it;
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          const math::Quatd stQ = stationOrient(st, stPos, timeDays);
          const math::Vec3d stV = stationVelKmS(st, timeDays);

          const double wz = st.radiusKm * 1.10;
          const double dockZ = wz - st.slotDepthKm - st.radiusKm * 0.25;
          const math::Vec3d dockLocal{0,0, dockZ};
          ship.setPositionKm(stPos + stQ.rotate(dockLocal));
          ship.setVelocityKmS(stV);
          ship.setOrientation(stQ * math::Quatd::fromAxisAngle({0,1,0}, math::kPi)); // face outward-ish
        }
      }

      // --- Fuel burn (sim-time) ---
      if (!docked && fsdState == FsdState::Idle) {
        if (input.boost) {
          fuel -= dtSim * 0.0025;
        }
        if (supercruise) {
          const double v = ship.velocityKmS().length();
          fuel -= dtSim * (0.0020 + v * 0.00000015);
        }

        if (fuel < 0.0) fuel = 0.0;
        if (fuel > fuelMax) fuel = fuelMax;

        if (fuel <= 0.0 && supercruise) {
          supercruise = false;
          toast(toasts, "Fuel depleted: supercruise disengaged.", 2.0);
        }
      }

      const bool combatSimEnabled = (fsdState == FsdState::Idle && !supercruise);
      if (combatSimEnabled) {
      // Spawn pirates occasionally (simple "threat" loop).
      if (timeDays >= nextPirateSpawnDays && contacts.size() < 8) {
        nextPirateSpawnDays = timeDays + (rng.range(120.0, 220.0) / 86400.0); // every ~2-4 minutes

        Contact p{};
        p.id = std::max<core::u64>(1, rng.nextU64());
        p.pirate = true;
        p.name = "Pirate " + std::to_string((int)contacts.size() + 1);
        p.ship.setMaxLinearAccelKmS2(0.06);
        p.ship.setMaxAngularAccelRadS2(0.9);

        // Spawn somewhere near the player, but not too close.
        const math::Vec3d randDir = math::Vec3d{rng.range(-1.0,1.0), rng.range(-0.3,0.3), rng.range(-1.0,1.0)}.normalized();
        const double distKm = rng.range(50000.0, 120000.0);
        p.ship.setPositionKm(ship.positionKm() + randDir * distKm);
        p.ship.setVelocityKmS(ship.velocityKmS());
        p.ship.setOrientation(quatFromTo({0,0,1}, (-randDir).normalized()));

        contacts.push_back(std::move(p));
        toast(toasts, "Contact: pirate detected!", 3.0);
      }

      // Ensure mission bounty targets exist in their target system.
      if (!docked && !missions.empty()) {
        for (const auto& m : missions) {
          if (m.completed || m.failed) continue;
          if (!((m.type == sim::MissionType::BountyScan) || (m.type == sim::MissionType::BountyKill))) {
            continue;
          }
          if (m.toSystem != currentSystem->stub.id) continue;
          if (m.targetNpcId == 0) continue;

          bool present = false;
          for (auto& c : contacts) {
            if (c.alive && c.id == m.targetNpcId) {
              c.missionTarget = true;
              present = true;
              break;
            }
          }
          if (present) continue;

          // Spawn a distinct target pirate.
          Contact tgt{};
          tgt.id = m.targetNpcId;
          tgt.pirate = true;
          tgt.missionTarget = true;
          tgt.name = "Bounty Target";
          tgt.shield = 80.0;
          tgt.hull = 90.0;
          tgt.ship.setMaxLinearAccelKmS2(0.07);
          tgt.ship.setMaxAngularAccelRadS2(1.0);

          const math::Vec3d randDir = math::Vec3d{rng.range(-1.0,1.0), rng.range(-0.2,0.2), rng.range(-1.0,1.0)}.normalized();
          const double distKm = rng.range(60000.0, 140000.0);
          tgt.ship.setPositionKm(ship.positionKm() + randDir * distKm);
          tgt.ship.setVelocityKmS(ship.velocityKmS());
          tgt.ship.setOrientation(quatFromTo({0,0,1}, (-randDir).normalized()));

          contacts.push_back(std::move(tgt));
          break; // spawn at most one target per frame
        }
      }

      // Contacts AI + combat
      for (auto& c : contacts) {
        if (!c.alive) continue;

        // cooldowns
        c.fireCooldown = std::max(0.0, c.fireCooldown - dtSim);
        if (c.pirate) {
          // chase player
          sim::ShipInput ai{};
          ai.dampers = true;

          const math::Vec3d to = ship.positionKm() - c.ship.positionKm();
          const double dist = to.length();
          const math::Vec3d toN = (dist > 1e-6) ? to / dist : math::Vec3d{0,0,1};

          // Face player
          const math::Vec3d desiredFwdWorld = toN;
          const math::Vec3d desiredFwdLocal = c.ship.orientation().conjugate().rotate(desiredFwdWorld);
          const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
          const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);
          ai.torqueLocal.x = std::clamp(pitchErr * 1.8, -1.0, 1.0);
          ai.torqueLocal.y = std::clamp(yawErr * 1.8, -1.0, 1.0);

          // Thrust: try to keep a standoff distance.
          const double desiredDist = 35000.0;
          double vAim = 0.0;
          if (dist > desiredDist) vAim = std::min(0.22, 0.000004 * (dist - desiredDist));
          if (dist < desiredDist*0.6) vAim = -0.12;

          const math::Vec3d desiredVel = ship.velocityKmS() + toN * vAim;
          const math::Vec3d dv = desiredVel - c.ship.velocityKmS();
          const math::Vec3d accelWorldDir = (dv.lengthSq() > 1e-9) ? dv.normalized() : math::Vec3d{0,0,0};
          ai.thrustLocal = c.ship.orientation().conjugate().rotate(accelWorldDir);

          c.ship.step(dtSim, ai);

          // Fire if aligned
          if (c.fireCooldown <= 0.0 && dist < 90000.0) {
            const double aim = math::dot(c.ship.forward().normalized(), toN);
            if (aim > 0.992) {
              c.fireCooldown = 0.35;
              const double dmg = 11.0;
              applyDamage(dmg, playerShield, playerHull);

              const math::Vec3d aKm = c.ship.positionKm();
              const math::Vec3d bKm = ship.positionKm();
              beams.push_back({toRenderU(aKm), toRenderU(bKm), 0.95f, 0.45f, 0.10f, 0.08});
            }
          }
        }
      }

      // Station turrets: in station vicinity, help against pirates.
      for (const auto& st : currentSystem->stations) {
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const double zoneKm = st.radiusKm * 25.0;
        const double distShip = (ship.positionKm() - stPos).length();
        if (distShip > zoneKm) continue; // only if player is nearby

        // Shoot a pirate every so often
        for (auto& c : contacts) {
          if (!c.alive || !c.pirate) continue;
          const double d = (c.ship.positionKm() - stPos).length();
          if (d < zoneKm) {
            // apply light damage (feel of station defenses)
            applyDamage(4.0 * dtSim, c.shield, c.hull);
            if (c.hull <= 0.0) {
              c.alive = false;
              credits += 250.0;
              toast(toasts, "Station defenses destroyed a pirate (+250).", 2.5);
            }
          }
        }
      }

      }

      // Collisions (player with station hull)
      if (!docked) {
        for (const auto& st : currentSystem->stations) {
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          const math::Quatd stQ = stationOrient(st, stPos, timeDays);
          const math::Vec3d relLocal = stQ.conjugate().rotate(ship.positionKm() - stPos);

          if (insideStationHullExceptSlot(st, relLocal)) {
            // Damage based on relative speed and push outward slightly.
            const math::Vec3d stV = stationVelKmS(st, timeDays);
            const double relSpeed = (ship.velocityKmS() - stV).length();
            applyDamage(relSpeed * 18.0, playerShield, playerHull);

            // Push out along local axis with largest penetration.
            const double wx = st.radiusKm * 0.70;
            const double wy = st.radiusKm * 0.70;
            const double wz = st.radiusKm * 1.10;

            double dx = wx - std::abs(relLocal.x);
            double dy = wy - std::abs(relLocal.y);
            double dz = wz - std::abs(relLocal.z);

            math::Vec3d pushLocal{0,0,0};
            if (dx <= dy && dx <= dz) pushLocal.x = (relLocal.x >= 0 ? 1 : -1) * (dx + 200.0);
            else if (dy <= dz) pushLocal.y = (relLocal.y >= 0 ? 1 : -1) * (dy + 200.0);
            else pushLocal.z = (relLocal.z >= 0 ? 1 : -1) * (dz + 200.0);

            ship.setPositionKm(ship.positionKm() + stQ.rotate(pushLocal));
            ship.setVelocityKmS(stV); // kill relative motion on impact
            toast(toasts, "Collision!", 1.2);
            break;
          }
        }
      }

      // --- Bounty scan progress ---
      if (scanning && !docked && fsdState == FsdState::Idle && !supercruise) {
        bool valid = false;
        if (target.kind == TargetKind::Contact && target.index < (int)contacts.size()) {
          auto& c = contacts[(std::size_t)target.index];
          if (c.alive) {
            const double dist = (c.ship.positionKm() - ship.positionKm()).length();
            if (dist <= scanRangeKm) {
              // Find an active bounty-scan mission for this target.
              for (auto& m : missions) {
                if (m.completed || m.failed) continue;
                if (m.type != sim::MissionType::BountyScan) continue;
                if (m.toSystem != currentSystem->stub.id) continue;
                if (m.targetNpcId != c.id) continue;

                valid = true;
                scanProgressSec += dtReal;
                if (scanProgressSec >= scanRequiredSec && !m.scanned) {
                  m.scanned = true;
                  m.completed = true;
                  credits += m.reward;
                  addRep(m.factionId, +2.0);
                  toast(toasts, "Mission complete: bounty scan uploaded! +" + std::to_string((int)m.reward) + " cr", 3.0);
                  scanning = false;
                  scanProgressSec = 0.0;
                }
                break;
              }
            }
          }
        }

        if (!valid) {
          // Cancel silently if the target / mission is no longer valid.
          scanning = false;
          scanProgressSec = 0.0;
        }
      } else if (!scanning) {
        scanProgressSec = 0.0;
      }

      timeDays += dtSim / 86400.0;

      // --- Mission deadlines / docked completion ---
      if (!missions.empty()) {
        for (auto& m : missions) {
          if (m.completed || m.failed) continue;
          if (m.deadlineDay > 0.0 && timeDays > m.deadlineDay) {
            m.failed = true;
            addRep(m.factionId, -4.0);
            toast(toasts, "Mission failed: deadline missed.", 3.0);
          }
        }

        if (docked && dockedStationId != 0) {
          const sim::StationId here = dockedStationId;
          const sim::SystemId sysId = currentSystem->stub.id;
          for (auto& m : missions) {
            if (m.completed || m.failed) continue;

            const bool atFinal = (sysId == m.toSystem && here == m.toStation);
            const bool atVia = (m.viaSystem != 0 && sysId == m.viaSystem && here == m.viaStation);

            if (m.type == sim::MissionType::Courier) {
              if (atFinal) {
                m.completed = true;
                credits += m.reward;
                addRep(m.factionId, +2.0);
                toast(toasts, "Mission complete: courier delivery. +" + std::to_string((int)m.reward) + " cr", 3.0);
              }
            } else if (m.type == sim::MissionType::Delivery || m.type == sim::MissionType::MultiDelivery) {
              if (m.type == sim::MissionType::MultiDelivery && m.viaSystem != 0 && m.leg == 0 && atVia) {
                m.leg = 1;
                toast(toasts, "Multi-hop delivery: leg 1/2 complete.", 2.5);
              } else if (atFinal && (m.viaSystem == 0 || m.leg >= 1)) {
                const auto cid = (econ::CommodityId)m.commodity;
                const double have = cargo[(std::size_t)cid];
                if (have + 1e-6 >= m.units) {
                  cargo[(std::size_t)cid] -= m.units;
                  if (cargo[(std::size_t)cid] < 1e-6) cargo[(std::size_t)cid] = 0.0;

                  m.completed = true;
                  credits += m.reward;
                  addRep(m.factionId, +2.0);
                  toast(toasts, "Mission complete: delivery. +" + std::to_string((int)m.reward) + " cr", 3.0);
                }
              }
            }
          }
        }
      }

      playerLaserCooldown = std::max(0.0, playerLaserCooldown - dtSim);

      // Shield regen (slow)
      if (!paused && playerHull > 0.0 && playerShield < 100.0) {
        playerShield = std::min(100.0, playerShield + 2.5 * (dtSim / 60.0));
      }
    }

    // Death / respawn
    if (playerHull <= 0.0) {
      toast(toasts, "Ship destroyed! Respawning (lost cargo, -10% credits).", 4.0);
      playerHull = 100.0;
      playerShield = 60.0;
      credits *= 0.90;
      cargo.fill(0.0);
      fuel = fuelMax;
      supercruise = false;
      scanning = false;
      scanProgressSec = 0.0;
      fsdState = FsdState::Idle;
      fsdTargetSystem = 0;
      fsdChargeRemainingSec = 0.0;
      fsdTravelRemainingSec = 0.0;
      navAutoRun = false;
      beams.clear();
      contacts.clear();
      docked = true;
      if (!currentSystem->stations.empty()) {
        dockedStationId = currentSystem->stations.front().id;
        selectedStationIndex = 0;
      } else {
        dockedStationId = 0;
      }
      // Will snap to dock position next tick.
    }

    // Update beam TTL
    for (auto& b : beams) b.ttl -= dtReal;
    beams.erase(std::remove_if(beams.begin(), beams.end(), [](const Beam& b){ return b.ttl <= 0.0; }), beams.end());

    // Toast TTL
    for (auto& t : toasts) t.ttl -= dtReal;
    toasts.erase(std::remove_if(toasts.begin(), toasts.end(), [](const ToastMsg& t){ return t.ttl <= 0.0; }), toasts.end());

    // ---- Camera follow (third-person) ----
    render::Camera cam;
    int w = 1280, h = 720;
    SDL_GetWindowSize(window, &w, &h);
    const double aspect = (h > 0) ? (double)w / (double)h : 16.0/9.0;

    cam.setPerspective(math::degToRad(60.0), aspect, 0.01, 20000.0);

    const math::Vec3d shipPosU = toRenderU(ship.positionKm());
    const math::Vec3d back = ship.forward() * (-6.0);
    const math::Vec3d up = ship.up() * (2.0);

    cam.setPosition(shipPosU + back + up);
    cam.setOrientation(ship.orientation());

    const math::Mat4d view = cam.viewMatrix();
    const math::Mat4d proj = cam.projectionMatrix();

    float viewF[16], projF[16];
    matToFloat(view, viewF);
    matToFloat(proj, projF);

    meshRenderer.setViewProj(viewF, projF);
    lineRenderer.setViewProj(viewF, projF);

    // ---- Build instances (star + planets) ----
    std::vector<render::InstanceData> spheres;
    spheres.reserve(1 + currentSystem->planets.size());

    // Star at origin
    {
      const double starRadiusKm = currentSystem->star.radiusSol * kSOLAR_RADIUS_KM;
      const float starScale = (float)std::max(0.8, (starRadiusKm / kRENDER_UNIT_KM) * 3.0);
      spheres.push_back(makeInstUniform({0,0,0}, starScale, 1.0f, 0.95f, 0.75f));
    }

    // Planets
    for (const auto& p : currentSystem->planets) {
      const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, timeDays);
      const math::Vec3d posKm = posAU * kAU_KM;
      const math::Vec3d posU = toRenderU(posKm);

      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const float scale = (float)std::max(0.25, (radiusKm / kRENDER_UNIT_KM) * 200.0);

      // Simple color palette by type
      float cr=0.6f, cg=0.6f, cb=0.6f;
      switch (p.type) {
        case sim::PlanetType::Rocky: cr=0.6f; cg=0.55f; cb=0.5f; break;
        case sim::PlanetType::Desert: cr=0.8f; cg=0.7f; cb=0.35f; break;
        case sim::PlanetType::Ocean: cr=0.25f; cg=0.45f; cb=0.85f; break;
        case sim::PlanetType::Ice: cr=0.7f; cg=0.85f; cb=0.95f; break;
        case sim::PlanetType::GasGiant: cr=0.7f; cg=0.55f; cb=0.35f; break;
        default: break;
      }

      spheres.push_back(makeInstUniform(posU, scale, cr,cg,cb));
    }

    // Orbit lines (planets)
    std::vector<render::LineVertex> lines;
    lines.reserve(currentSystem->planets.size() * 128 + currentSystem->stations.size() * 64 + beams.size() * 2);

    for (const auto& p : currentSystem->planets) {
      const int seg = 96;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * p.orbit.periodDays;
        const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, t);
        const math::Vec3d posU = toRenderU(posAU * kAU_KM);

        if (s > 0) {
          lines.push_back({(float)prev.x,(float)prev.y,(float)prev.z, 0.22f,0.22f,0.25f});
          lines.push_back({(float)posU.x,(float)posU.y,(float)posU.z, 0.22f,0.22f,0.25f});
        }
        prev = posU;
      }
    }

    // Station orbit lines + slot axis lines
    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      const auto& st = currentSystem->stations[i];

      // orbit (fewer segs, stations are secondary)
      const int seg = 48;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * st.orbit.periodDays;
        const math::Vec3d posU = toRenderU(sim::orbitPosition3DAU(st.orbit, t) * kAU_KM);
        if (s > 0) {
          lines.push_back({(float)prev.x,(float)prev.y,(float)prev.z, 0.16f,0.18f,0.22f});
          lines.push_back({(float)posU.x,(float)posU.y,(float)posU.z, 0.16f,0.18f,0.22f});
        }
        prev = posU;
      }

      // approach corridor axis for targeted station
      if (target.kind == TargetKind::Station && target.index == i) {
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const math::Quatd stQ = stationOrient(st, stPos, timeDays);
        const math::Vec3d axis = stQ.rotate({0,0,1});
        const math::Vec3d aU = toRenderU(stPos + axis * (st.radiusKm * 1.1));
        const math::Vec3d bU = toRenderU(stPos + axis * (st.radiusKm * 1.1 + st.approachLengthKm));
        lines.push_back({(float)aU.x,(float)aU.y,(float)aU.z, 0.35f,0.85f,0.45f});
        lines.push_back({(float)bU.x,(float)bU.y,(float)bU.z, 0.35f,0.85f,0.45f});
      }
    }

    // Laser beams
    for (const auto& b : beams) {
      lines.push_back({(float)b.aU.x,(float)b.aU.y,(float)b.aU.z, b.r,b.g,b.b});
      lines.push_back({(float)b.bU.x,(float)b.bU.y,(float)b.bU.z, b.r,b.g,b.b});
    }

    // Station geometry (cubes)
    std::vector<render::InstanceData> cubes;
    cubes.reserve(1 + currentSystem->stations.size() * 18 + contacts.size());

    for (const auto& st : currentSystem->stations) {
      const math::Vec3d stPos = stationPosKm(st, timeDays);
      const math::Quatd stQ = stationOrient(st, stPos, timeDays);
      emitStationGeometry(st, stPos, stQ, cubes);
    }

    // Ship instance (cube, rotated)
    cubes.push_back(makeInst(toRenderU(ship.positionKm()),
                             {0.35, 0.20, 0.60},
                             ship.orientation(),
                             0.90f, 0.90f, 1.00f));

    // Contacts (pirates)
    for (const auto& c : contacts) {
      if (!c.alive) continue;
      const float r = c.pirate ? 1.0f : 0.6f;
      const float g = c.pirate ? 0.25f : 0.7f;
      const float b = c.pirate ? 0.25f : 0.8f;
      cubes.push_back(makeInst(toRenderU(c.ship.positionKm()),
                               {0.25, 0.18, 0.45},
                               c.ship.orientation(),
                               r,g,b));
    }

    // ---- Render ---
    glViewport(0, 0, w, h);
    glClearColor(0.01f, 0.01f, 0.02f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Lines (orbits, corridor, beams)
    lineRenderer.drawLines(lines);

    // Star + planets
    meshRenderer.setMesh(&sphere);
    meshRenderer.drawInstances(spheres);

    // Cubes (stations, ship, contacts)
    meshRenderer.setMesh(&cube);
    meshRenderer.drawInstances(cubes);

    // ---- UI ----
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame(window);
    ImGui::NewFrame();

    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    // HUD overlay: target marker + crosshair + toasts
    {
      ImDrawList* draw = ImGui::GetForegroundDrawList();
      const ImVec2 center((float)w * 0.5f, (float)h * 0.5f);
      draw->AddLine({center.x - 8, center.y}, {center.x + 8, center.y}, IM_COL32(160,160,170,140), 1.0f);
      draw->AddLine({center.x, center.y - 8}, {center.x, center.y + 8}, IM_COL32(160,160,170,140), 1.0f);

      // Target marker
      std::optional<math::Vec3d> tgtKm{};
      std::string tgtLabel;

      if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
        const auto& st = currentSystem->stations[target.index];
        tgtKm = stationPosKm(st, timeDays);
        tgtLabel = st.name;
      } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
        const auto& p = currentSystem->planets[target.index];
        tgtKm = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
        tgtLabel = p.name;
      } else if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
        const auto& c = contacts[target.index];
        if (c.alive) {
          tgtKm = c.ship.positionKm();
          tgtLabel = c.name;
        }
      }

      if (tgtKm) {
        ImVec2 px{};
        if (projectToScreen(toRenderU(*tgtKm), view, proj, w, h, px)) {
          const double distKm = (*tgtKm - ship.positionKm()).length();
          draw->AddCircle({px.x, px.y}, 14.0f, IM_COL32(255,170,80,190), 1.5f);
          std::string s = tgtLabel + "  " + std::to_string((int)distKm) + " km";
          draw->AddText({px.x + 18, px.y - 8}, IM_COL32(255,210,170,210), s.c_str());
        }
      }

      // Toasts
      float y = 18.0f;
      for (const auto& t : toasts) {
        draw->AddText({18.0f, y}, IM_COL32(240,240,240,220), t.text.c_str());
        y += 18.0f;
      }
    }

    if (showShip) {
      ImGui::Begin("Ship / Flight");

      const auto pos = ship.positionKm();
      const auto vel = ship.velocityKmS();
      const auto wv  = ship.angularVelocityRadS();

      ImGui::Text("Time: %.2f days  (x%.1f) %s", timeDays, timeScale, paused ? "[PAUSED]" : "");
      ImGui::SliderFloat("Time scale (sim sec / real sec)", (float*)&timeScale, 0.0f, 2000.0f);

      ImGui::Separator();
      ImGui::Text("Hull: %.0f / 100   Shield: %.0f / 100", playerHull, playerShield);
      ImGui::Text("Laser cooldown: %.2fs", playerLaserCooldown);
      ImGui::Text("Fuel: %.1f / %.1f", fuel, fuelMax);

      {
        const double cap = std::max(1.0, cargoCapacityKg);
        const double load = std::clamp(cargoMassKg(cargo) / cap, 0.0, 1.0);
        const double linMult = (1.0 - 0.25 * load);
        const double angMult = (1.0 - 0.20 * load);
        ImGui::Text("Cargo load: %.0f%%   Handling: x%.2f accel, x%.2f turn", load * 100.0, linMult, angMult);
      }

      const double jrMax = fsdBaseRangeLy();
      const double jrNow = fsdCurrentRangeLy();
      ImGui::Text("FSD range: %.1f ly (current) / %.1f ly (max)", jrNow, jrMax);
      if (fsdState == FsdState::Charging) {
        ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.25f, 1.0f), "FSD CHARGING (%.1fs)", fsdChargeRemainingSec);
      } else if (fsdState == FsdState::Travelling) {
        ImGui::TextColored(ImVec4(0.7f, 0.9f, 1.0f, 1.0f), "IN HYPERSPACE (%.1fs)", fsdTravelRemainingSec);
      } else {
        const bool ready = (timeDays >= fsdReadyDay);
        ImGui::Text("FSD: %s  (J to jump)", ready ? "READY" : "COOLDOWN");
      }

      if (docked) {
        ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.5f, 1.0f), "DOCKED");
      } else {
        ImGui::TextDisabled("Not docked");
      }

      ImGui::Checkbox("Autopilot (P)", &autopilot);

      // Mouse flight settings
      bool mf = mouseFlight;
      if (ImGui::Checkbox("Mouse flight (M/RMB)", &mf)) {
        setMouseFlight(mf);
      }
      ImGui::SameLine();
      ImGui::Checkbox("Invert mouse Y", &invertMouseY);
      ImGui::SliderFloat("Mouse sensitivity", &mouseSensitivity, 0.002f, 0.030f, "%.4f");

      if (supercruise) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(0.6f, 0.9f, 1.0f, 1.0f), "SUPERCRUISE");
      }

      if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
        ImGui::Text("Route: hop %d/%d  (J: next hop)", (int)navRouteHop + 1, (int)navRoute.size() - 1);
      }

      if (scanning) {
        ImGui::TextColored(ImVec4(0.95f, 0.9f, 0.6f, 1.0f), "Scanning... %.0f%% (K to cancel)", (scanProgressSec / scanDurationSec) * 100.0);
      }

      if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
        const auto& st = currentSystem->stations[target.index];
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const double dist = (pos - stPos).length();
        auto it = clearances.find(st.id);
        bool clearance = (it != clearances.end() && it->second.granted && timeDays <= it->second.expiresDays);
        ImGui::Separator();
        ImGui::Text("Target: %s (%.0f km)", st.name.c_str(), dist);
        ImGui::Text("Clearance: %s", clearance ? "GRANTED" : "NONE");
        ImGui::TextDisabled("Docking: request clearance (L), fly through slot, then press G.");
      } else {
        ImGui::TextDisabled("Target: (none)  [T stations, B planets, N contacts]");
      }

      ImGui::Separator();
      ImGui::Text("Pos (km):   [%.1f %.1f %.1f]", pos.x, pos.y, pos.z);
      ImGui::Text("Vel (km/s): [%.3f %.3f %.3f] |v|=%.3f", vel.x, vel.y, vel.z, vel.length());
      ImGui::Text("AngVel (rad/s): [%.3f %.3f %.3f]", wv.x, wv.y, wv.z);

      ImGui::TextDisabled("Controls:");
      ImGui::BulletText("Translate: WASD + R/F (up/down)");
      ImGui::BulletText("Rotate: Arrow keys + Q/E roll");
      ImGui::BulletText("Mouse flight: Right Mouse / M (Esc to release)");
      ImGui::BulletText("Boost: LShift   Brake: X");
      ImGui::BulletText("Dampers: Z (on) / C (off)");
      ImGui::BulletText("Target: T (stations), B (planets), N (contacts), Y (clear)");
      ImGui::BulletText("Docking: L (request clearance), G (dock/undock)");
      ImGui::BulletText("Supercruise-lite: H (to target)");
      ImGui::BulletText("FSD jump: J (route hop or selected system)");
      ImGui::BulletText("Scan (bounty missions): K");
      ImGui::BulletText("Fire: Left Mouse");
      ImGui::BulletText("Pause: Space   Save: F5   Load: F9");
      ImGui::BulletText("Missions: F4");

      ImGui::End();
    }

    if (showEconomy) {
      beginStationSelectorHUD(*currentSystem, selectedStationIndex, docked, dockedStationId);

      if (!currentSystem->stations.empty()) {
        const auto& station = currentSystem->stations[(std::size_t)selectedStationIndex];
        auto& stEcon = universe.stationEconomy(station, timeDays);
        const double rep = getRep(station.factionId);
        const double feeEff = effectiveFeeRate(station);
        double cargoKgNow = cargoMassKg(cargo);

        ImGui::Begin("Market Details");
        ImGui::Text("Credits: %.2f", credits);

        ImGui::Text("Station: %s (%s)", station.name.c_str(), stationTypeName(station.type));
        ImGui::Text("Faction: %s   Rep: %.0f", factionName(station.factionId).c_str(), rep);
        ImGui::Text("Fees: base %.2f%%  effective %.2f%%", station.feeRate * 100.0, feeEff * 100.0);
        ImGui::Text("Cargo: %.0f / %.0f kg", cargoKgNow, cargoCapacityKg);

        const bool canTrade = docked && (station.id == dockedStationId);
        if (!canTrade) {
          ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.35f, 1.0f), "Trade disabled: dock at this station to buy/sell.");
        }

        // Simple docked services
        if (canTrade) {
          const double hullMissing = std::max(0.0, 100.0 - playerHull);
          const double repairBase = hullMissing * 12.0;
          const double repairCost = repairBase * (1.0 + feeEff);

          if (ImGui::Button("Repair hull")) {
            if (hullMissing <= 0.01) {
              toast(toasts, "Hull already at 100%.", 2.0);
            } else if (credits >= repairCost) {
              credits -= repairCost;
              playerHull = 100.0;
              toast(toasts, "Ship repaired.", 2.0);
            } else {
              toast(toasts, "Not enough credits for repair.", 2.0);
            }
          }
          ImGui::SameLine();
          ImGui::TextDisabled("(%.0f cr)", repairCost);

          // Refuel service: buys Fuel commodity into the ship tank.
          const auto fuelQuote = econ::quote(stEcon, station.economyModel, econ::CommodityId::Fuel, 0.10);
          const double fuelNeed = std::max(0.0, fuelMax - fuel);
          const double fuelAvail = std::max(0.0, fuelQuote.inventory);
          const double fuelBuy = std::min(fuelNeed, fuelAvail);
          const double fuelCost = fuelBuy * fuelQuote.ask * (1.0 + feeEff);

          if (ImGui::Button("Refuel")) {
            if (fuelNeed <= 0.01) {
              toast(toasts, "Fuel tank already full.", 2.0);
            } else if (fuelBuy <= 0.01) {
              toast(toasts, "Station is out of fuel.", 2.0);
            } else if (credits >= fuelCost) {
              credits -= fuelCost;
              fuel += fuelBuy;
              // Reduce station inventory (best-effort).
              stEcon.inventory[(std::size_t)econ::CommodityId::Fuel] = std::max(0.0, fuelAvail - fuelBuy);
              toast(toasts, "Refueled.", 2.0);
            } else {
              toast(toasts, "Not enough credits to refuel.", 2.0);
            }
          }
          ImGui::SameLine();
          ImGui::TextDisabled("(%.1f units, %.0f cr)", fuelBuy, fuelCost);

          // Shipyard upgrades
          if (station.type == econ::StationType::Shipyard) {
            ImGui::Separator();
            ImGui::Text("Shipyard Upgrades");

            const double cargoUpgradeKg = 200.0;
            const double cargoUpgradeCost = 2000.0 * (1.0 + feeEff);
            if (ImGui::Button("Cargo racks +200kg")) {
              if (credits >= cargoUpgradeCost) {
                credits -= cargoUpgradeCost;
                cargoCapacityKg += cargoUpgradeKg;
                toast(toasts, "Cargo capacity upgraded.", 2.0);
              } else {
                toast(toasts, "Not enough credits.", 2.0);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", cargoUpgradeCost);

            const double fuelUpgrade = 10.0;
            const double fuelUpgradeCost = 2500.0 * (1.0 + feeEff);
            if (ImGui::Button("Fuel tank +10")) {
              if (credits >= fuelUpgradeCost) {
                credits -= fuelUpgradeCost;
                fuelMax += fuelUpgrade;
                toast(toasts, "Fuel tank upgraded.", 2.0);
              } else {
                toast(toasts, "Not enough credits.", 2.0);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", fuelUpgradeCost);

            const double rangeUpgrade = 2.0;
            const double rangeUpgradeCost = 9000.0 * (1.0 + feeEff);
            if (ImGui::Button("FSD tuning +2ly")) {
              if (credits >= rangeUpgradeCost) {
                credits -= rangeUpgradeCost;
                fsdRangeLy += rangeUpgrade;
                toast(toasts, "FSD range improved.", 2.0);
              } else {
                toast(toasts, "Not enough credits.", 2.0);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", rangeUpgradeCost);
          }
        }

        static int selectedCommodity = 0;
        ImGui::SliderInt("Plot commodity", &selectedCommodity, 0, (int)econ::kCommodityCount - 1);

        // Table
        if (ImGui::BeginTable("market", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
          ImGui::TableSetupColumn("Commodity");
          ImGui::TableSetupColumn("Inv");
          ImGui::TableSetupColumn("Ask");
          ImGui::TableSetupColumn("Bid");
          ImGui::TableSetupColumn("Cargo");
          ImGui::TableSetupColumn("Trade");
          ImGui::TableHeadersRow();

          for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
            const auto cid = (econ::CommodityId)i;
            const auto q = econ::quote(stEcon, station.economyModel, cid, 0.10);

            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", std::string(econ::commodityName(cid)).c_str());

            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%.0f", q.inventory);

            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.2f", q.ask);

            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.2f", q.bid);

            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.0f", cargo[i]);

            ImGui::TableSetColumnIndex(5);
            ImGui::PushID((int)i);

            static float qty[ (int)econ::kCommodityCount ] = {};
            if (qty[i] <= 0.0f) qty[i] = 10.0f;
            ImGui::SetNextItemWidth(70);
            ImGui::InputFloat("##qty", &qty[i], 1.0f, 10.0f, "%.0f");

            ImGui::SameLine();
            ImGui::BeginDisabled(!canTrade);
            if (ImGui::SmallButton("Buy")) {
              const double buyUnits = std::max(0.0, (double)qty[i]);
              const double addKg = buyUnits * econ::commodityDef(cid).massKg;
              if (cargoKgNow + addKg > cargoCapacityKg + 1e-6) {
                toast(toasts, "Cargo hold full (mass limit).", 2.0);
              } else {
                auto tr = econ::buy(stEcon, station.economyModel, cid, buyUnits, credits, 0.10, feeEff);
                if (tr.ok) {
                  cargo[i] += buyUnits;
                  cargoKgNow += addKg;
                }
              }
            }

            ImGui::SameLine();
            if (ImGui::SmallButton("Sell")) {
              const double sellUnits = std::min<double>(qty[i], cargo[i]);
              if (sellUnits > 0.0) {
                auto tr = econ::sell(stEcon, station.economyModel, cid, sellUnits, credits, 0.10, feeEff);
                if (tr.ok) {
                  cargo[i] -= sellUnits;
                  cargoKgNow = std::max(0.0, cargoKgNow - sellUnits * econ::commodityDef(cid).massKg);
                }
              }
            }
            ImGui::EndDisabled();

            ImGui::PopID();
          }

          ImGui::EndTable();
        }

        // Price history plot for selected commodity
        const std::size_t cidx = (std::size_t)selectedCommodity;
        const auto& hist = stEcon.history[cidx];
        if (!hist.empty()) {
          std::vector<float> vals;
          vals.reserve(hist.size());
          for (const auto& p : hist) vals.push_back((float)p.price);

          ImGui::PlotLines("Price history", vals.data(), (int)vals.size(), 0, nullptr, 0.0f, 0.0f, ImVec2(0, 120));
        } else {
          ImGui::TextDisabled("No history yet (time needs to advance).");
        }

        ImGui::End();
      }
    }

    if (showMissions) {
      ImGui::Begin("Missions");

      auto stationNameById = [&](sim::SystemId sysId, sim::StationId stId) -> std::string {
        if (stId == 0) return "";
        const auto& sys = universe.getSystem(sysId);
        for (const auto& st : sys.stations) {
          if (st.id == stId) return st.name;
        }
        return "Station #" + std::to_string((std::uint64_t)stId);
      };

      auto systemNameById = [&](sim::SystemId sysId) -> std::string {
        if (sysId == 0) return "";
        return universe.getSystem(sysId).stub.name;
      };

      auto describeMission = [&](const sim::Mission& m) -> std::string {
        std::string out;
        switch (m.type) {
          case sim::MissionType::Courier: {
            out = "Courier: Deliver data to " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::Delivery: {
            const auto cid = (econ::CommodityId)m.commodity;
            out = "Delivery: Deliver " + std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name)
                + " to " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::MultiDelivery: {
            const auto cid = (econ::CommodityId)m.commodity;
            out = "Multi-hop: Deliver " + std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name);
            if (m.viaSystem != 0 && m.viaStation != 0) {
              out += " via " + stationNameById(m.viaSystem, m.viaStation) + " (" + systemNameById(m.viaSystem) + ")";
            }
            out += " -> " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::BountyScan: {
            out = "Bounty: Scan wanted pirate in " + systemNameById(m.toSystem);
          } break;
          case sim::MissionType::BountyKill: {
            out = "Bounty: Eliminate wanted pirate in " + systemNameById(m.toSystem);
          } break;
          default:
            out = "Mission";
            break;
        }
        return out;
      };

      if (ImGui::BeginTabBar("missions_tabs")) {
        if (ImGui::BeginTabItem("Active")) {
          if (missions.empty()) {
            ImGui::TextDisabled("No missions accepted.");
          }

          for (auto& m : missions) {
            ImGui::PushID((int)m.id);
            const bool active = !(m.completed || m.failed);
            const char* status = m.completed ? "COMPLETED" : (m.failed ? "FAILED" : "ACTIVE");
            ImGui::Text("[%s] %s", status, describeMission(m).c_str());
            ImGui::TextDisabled("Reward %.0f cr | Faction: %s (rep %.1f)", m.reward, factionName(m.factionId).c_str(), getRep(m.factionId));

            if (m.deadlineDay > 0.0) {
              const double hrsLeft = (m.deadlineDay - timeDays) * 24.0;
              ImGui::TextDisabled("Deadline: day %.2f (%.1f h left)", m.deadlineDay, hrsLeft);
            }

            if (active && m.type == sim::MissionType::MultiDelivery && m.viaSystem != 0 && m.viaStation != 0) {
              ImGui::TextDisabled("Progress: leg %d / 2", (int)m.leg + 1);
            }

            if (active) {
              if (ImGui::SmallButton("Set destination")) {
                sim::SystemId destSys = m.toSystem;
                if (m.type == sim::MissionType::MultiDelivery && m.viaSystem != 0 && m.leg == 0) {
                  destSys = m.viaSystem;
                }
                galaxySelectedSystemId = destSys;
                showGalaxy = true;
                toast(toasts, "Galaxy: mission destination selected.", 2.0);
              }
              ImGui::SameLine();
              if (ImGui::SmallButton("Abandon")) {
                m.failed = true;
                addRep(m.factionId, -2.0);
                toast(toasts, "Mission abandoned.", 2.0);
              }
            }

            ImGui::Separator();
            ImGui::PopID();
          }

          ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Mission Board")) {
          if (!docked || dockedStationId == 0) {
            ImGui::TextDisabled("Dock at a station to browse missions.");
            ImGui::EndTabItem();
          } else {
            // Resolve docked station
            int dockedIdx = -1;
            for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
              if (currentSystem->stations[i].id == dockedStationId) {
                dockedIdx = (int)i;
                break;
              }
            }

            if (dockedIdx < 0) {
              ImGui::TextDisabled("Docked station not found.");
              ImGui::EndTabItem();
            } else {
              const auto& st = currentSystem->stations[(std::size_t)dockedIdx];
              const double rep = getRep(st.factionId);
              ImGui::Text("Station: %s", st.name.c_str());
              ImGui::Text("Faction: %s (rep %.1f)", factionName(st.factionId).c_str(), rep);

              const int dayStamp = (int)std::floor(timeDays);
              if (missionOffersStationId != st.id || missionOffersDayStamp != dayStamp) {
                missionOffersStationId = st.id;
                missionOffersDayStamp = dayStamp;
                missionOffers.clear();

                // Deterministic board per-station per-day.
                core::SplitMix64 mrng((core::u64)st.id * 1469598103934665603ull ^ (core::u64)dayStamp * 1099511628211ull ^ (core::u64)st.factionId);

                auto candidates = universe.queryNearby(currentSystem->stub.posLy, 160.0, 128);
                // Filter out current system & systems without stations.
                std::vector<sim::SystemStub> dests;
                dests.reserve(candidates.size());
                for (const auto& s : candidates) {
                  if (s.id == currentSystem->stub.id) continue;
                  if (s.stationCount == 0) continue;
                  dests.push_back(s);
                }

                const int baseCount = 6 + (rep >= 25.0 ? 1 : 0) + (rep >= 50.0 ? 1 : 0);
                const int offerCount = std::clamp(baseCount, 4, 9);
                for (int i = 0; i < offerCount && !dests.empty(); ++i) {
                  const int pick = (int)(mrng.nextU32() % (core::u32)dests.size());
                  const auto destStub = dests[(std::size_t)pick];
                  const auto& destSys = universe.getSystem(destStub.id, &destStub);
                  const auto& destSt = destSys.stations[(std::size_t)(mrng.nextU32() % (core::u32)destSys.stations.size())];

                  const double distLy = (destStub.posLy - currentSystem->stub.posLy).length();

                  sim::Mission m{};
                  m.id = 0;
                  m.factionId = st.factionId;
                  m.fromSystem = currentSystem->stub.id;
                  m.fromStation = st.id;
                  m.toSystem = destStub.id;
                  m.toStation = destSt.id;
                  m.deadlineDay = timeDays + 1.0 + distLy / 20.0;

                  // Weighted mission types.
                  const double r = mrng.nextUnit();
                  if (r < 0.35) {
                    m.type = sim::MissionType::Courier;
                    m.reward = 350.0 + distLy * 110.0;
                  } else if (r < 0.70) {
                    m.type = sim::MissionType::Delivery;
                    const auto cid = (econ::CommodityId)(mrng.nextU32() % econ::kCommodityCount);
                    m.commodity = (core::u8)cid;
                    m.units = 25 + (mrng.nextU32() % 120);
                    m.reward = 250.0 + distLy * 120.0 + (double)m.units * 6.0;
                    m.cargoProvided = (mrng.nextUnit() < 0.25);
                  } else if (r < 0.85 && dests.size() >= 2) {
                    m.type = sim::MissionType::MultiDelivery;
                    const auto cid = (econ::CommodityId)(mrng.nextU32() % econ::kCommodityCount);
                    m.commodity = (core::u8)cid;
                    m.units = 20 + (mrng.nextU32() % 100);
                    m.reward = 400.0 + distLy * 150.0 + (double)m.units * 8.0;
                    m.cargoProvided = (mrng.nextUnit() < 0.35);

                    // Pick a via system/station from remaining candidates.
                    const int vpick = (int)(mrng.nextU32() % (core::u32)dests.size());
                    const auto viaStub = dests[(std::size_t)vpick];
                    const auto& viaSys = universe.getSystem(viaStub.id, &viaStub);
                    const auto& viaSt = viaSys.stations[(std::size_t)(mrng.nextU32() % (core::u32)viaSys.stations.size())];
                    m.viaSystem = viaStub.id;
                    m.viaStation = viaSt.id;
                    m.deadlineDay = timeDays + 2.0 + distLy / 18.0;
                  } else if (r < 0.93) {
                    m.type = sim::MissionType::BountyScan;
                    m.targetNpcId = std::max<core::u64>(1, mrng.nextU64());
                    m.reward = 900.0 + distLy * 80.0;
                    m.deadlineDay = timeDays + 1.5 + distLy / 22.0;
                  } else {
                    m.type = sim::MissionType::BountyKill;
                    m.targetNpcId = std::max<core::u64>(1, mrng.nextU64());
                    m.reward = 1400.0 + distLy * 90.0;
                    m.deadlineDay = timeDays + 1.8 + distLy / 20.0;
                  }

                  // Small reward scaling with reputation (positive only).
                  const double repScale = 1.0 + std::clamp(rep, 0.0, 100.0) / 100.0 * 0.10;
                  m.reward *= repScale;

                  missionOffers.push_back(m);
                }
              }

              if (missionOffers.empty()) {
                ImGui::TextDisabled("No offers available right now.");
              }

              for (std::size_t i = 0; i < missionOffers.size(); ++i) {
                auto& offer = missionOffers[i];
                ImGui::PushID((int)i);

                ImGui::TextWrapped("%s", describeMission(offer).c_str());
                if (offer.deadlineDay > 0.0) {
                  const double hrsLeft = (offer.deadlineDay - timeDays) * 24.0;
                  ImGui::TextDisabled("Reward %.0f cr | Deadline in %.1f h", offer.reward, hrsLeft);
                } else {
                  ImGui::TextDisabled("Reward %.0f cr", offer.reward);
                }

                bool canAccept = (missions.size() < 16);
                if (offer.cargoProvided) {
                  const auto cid = (econ::CommodityId)offer.commodity;
                  const double addKg = (double)offer.units * econ::commodityDef(cid).massKg;
                  canAccept = canAccept && (cargoMassKg(cargo) + addKg <= cargoCapacityKg + 1e-6);
                }

                if (!canAccept) ImGui::BeginDisabled();
                if (ImGui::SmallButton("Accept")) {
                  sim::Mission m = offer;
                  m.id = nextMissionId++;
                  m.completed = false;
                  m.failed = false;
                  m.leg = 0;
                  m.scanned = false;

                  if (m.cargoProvided) {
                    const auto cid = (econ::CommodityId)m.commodity;
                    const double addKg = (double)m.units * econ::commodityDef(cid).massKg;
                    if (cargoMassKg(cargo) + addKg > cargoCapacityKg + 1e-6) {
                      toast(toasts, "Cargo too full to accept this mission.", 2.0);
                    } else {
                      cargo[(std::size_t)cid] += (double)m.units;
                      toast(toasts, "Mission accepted (cargo provided).", 2.0);
                      missions.push_back(m);
                      missionOffers.erase(missionOffers.begin() + (std::ptrdiff_t)i);
                      --i;
                    }
                  } else {
                    toast(toasts, "Mission accepted.", 2.0);
                    missions.push_back(m);
                    missionOffers.erase(missionOffers.begin() + (std::ptrdiff_t)i);
                    --i;
                  }
                }
                if (!canAccept) ImGui::EndDisabled();

                ImGui::Separator();
                ImGui::PopID();
              }

              ImGui::EndTabItem();
            }
          }
        }

        ImGui::EndTabBar();
      }

      ImGui::End();
    }

    if (showContacts) {
      ImGui::Begin("Contacts / Combat");

      int alivePirates = 0;
      for (const auto& c : contacts) if (c.alive && c.pirate) ++alivePirates;
      ImGui::Text("Alive pirates: %d", alivePirates);

      if (ImGui::Button("Panic: clear pirates")) {
        for (auto& c : contacts) c.alive = false;
        toast(toasts, "Contacts cleared.", 2.0);
      }

      ImGui::Separator();

      for (std::size_t i = 0; i < contacts.size(); ++i) {
        const auto& c = contacts[i];
        if (!c.alive) continue;

        ImGui::PushID((int)i);
        ImGui::Text("%s %s%s", c.name.c_str(), c.pirate ? "[PIRATE]" : "", c.missionTarget ? " [BOUNTY]" : "");
        ImGui::SameLine();
        if (ImGui::SmallButton("Target")) {
          target.kind = TargetKind::Contact;
          target.index = i;
        }
        ImGui::TextDisabled("Hull %.0f  Shield %.0f", c.hull, c.shield);
        ImGui::PopID();
      }

      ImGui::End();
    }

    if (showGalaxy) {
      ImGui::Begin("Galaxy / Streaming");

      const auto center = currentSystem->stub.posLy;
      static float radius = 200.0f;
      ImGui::SliderFloat("Query radius (ly)", &radius, 20.0f, 1200.0f);

      auto nearby = universe.queryNearby(center, radius, 128);

      ImGui::Text("Nearby systems: %d", (int)nearby.size());

      // Mini-map canvas
      const ImVec2 canvasSize = ImVec2(420, 420);
      ImGui::BeginChild("map", canvasSize, true, ImGuiWindowFlags_NoScrollbar);

      ImDrawList* draw = ImGui::GetWindowDrawList();
      const ImVec2 p0 = ImGui::GetCursorScreenPos();
      const ImVec2 p1 = ImVec2(p0.x + canvasSize.x, p0.y + canvasSize.y);
      const ImVec2 centerPx = ImVec2((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);

      // Background
      draw->AddRectFilled(p0, p1, IM_COL32(10, 10, 14, 255));
      draw->AddRect(p0, p1, IM_COL32(80, 80, 95, 255));

      auto toPx = [&](const math::Vec3d& posLy) -> ImVec2 {
        const math::Vec3d d = posLy - center;
        const float sx = (float)(d.x / (double)radius) * (canvasSize.x * 0.5f);
        const float sy = (float)(d.y / (double)radius) * (canvasSize.y * 0.5f);
        return ImVec2(centerPx.x + sx, centerPx.y + sy);
      };

      // Build a quick lookup of stub positions for route drawing.
      std::unordered_map<sim::SystemId, math::Vec3d> stubPosById;
      stubPosById.reserve(nearby.size());
      for (const auto& s : nearby) {
        stubPosById[s.id] = s.posLy;
      }

      // Jump range feedback (max range circle).
      const double jrMaxLy = fsdBaseRangeLy();
      const float jrPx = (float)(jrMaxLy / (double)radius) * (canvasSize.x * 0.5f);
      if (jrPx > 2.0f && jrPx < 5000.0f) {
        draw->AddCircle(centerPx, jrPx, IM_COL32(90, 150, 240, 120), 96, 1.5f);
      }

      // Plotted route overlay.
      if (navRoute.size() >= 2) {
        for (std::size_t i = 0; i + 1 < navRoute.size(); ++i) {
          auto itA = stubPosById.find(navRoute[i]);
          auto itB = stubPosById.find(navRoute[i + 1]);
          if (itA == stubPosById.end() || itB == stubPosById.end()) continue;
          const ImVec2 a = toPx(itA->second);
          const ImVec2 b = toPx(itB->second);
          draw->AddLine(a, b, IM_COL32(255, 140, 80, 200), 2.0f);
        }
      }

      // Star lanes: connect each system to 3 nearest neighbors in XY
      const int k = 3;
      for (std::size_t i = 0; i < nearby.size(); ++i) {
        struct N { std::size_t j; double d2; };
        std::vector<N> ns;
        ns.reserve(nearby.size());
        for (std::size_t j = 0; j < nearby.size(); ++j) if (j != i) {
          const auto di = nearby[j].posLy - nearby[i].posLy;
          const double d2 = di.x*di.x + di.y*di.y + di.z*di.z;
          ns.push_back({j, d2});
        }
        std::sort(ns.begin(), ns.end(), [](const N& a, const N& b){ return a.d2 < b.d2; });
        const int count = std::min<int>(k, (int)ns.size());

        const ImVec2 a = toPx(nearby[i].posLy);
        for (int n = 0; n < count; ++n) {
          const ImVec2 b = toPx(nearby[ns[n].j].posLy);
          draw->AddLine(a, b, IM_COL32(50, 80, 120, 100), 1.0f);
        }
      }

      // Systems
      sim::SystemId selected = galaxySelectedSystemId;
      for (const auto& s : nearby) {
        const ImVec2 p = toPx(s.posLy);
        const bool isCurrent = (s.id == currentSystem->stub.id);
        const bool isSel = (s.id == selected);

        const double distLy = (s.posLy - currentSystem->stub.posLy).length();
        const bool inJumpRange = (distLy <= jrMaxLy + 1e-9);

        ImU32 col = isCurrent ? IM_COL32(255, 240, 160, 255) : IM_COL32(170, 170, 190, 255);
        if (s.factionId != 0) col = IM_COL32(160, 220, 170, 255);
        if (isSel) col = IM_COL32(255, 120, 120, 255);

        draw->AddCircleFilled(p, isCurrent ? 5.5f : 4.0f, col);

        // Jump range ring
        if (!isCurrent && inJumpRange) {
          draw->AddCircle(p, 6.5f, IM_COL32(80, 140, 220, 140), 24, 1.0f);
        }

        // Click detection
        const float rClick = 6.0f;
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
          const ImVec2 mp = ImGui::GetIO().MousePos;
          const float dx = mp.x - p.x;
          const float dy = mp.y - p.y;
          if (dx*dx + dy*dy <= rClick*rClick) selected = s.id;
        }
      }

      galaxySelectedSystemId = selected;

      ImGui::EndChild();

      ImGui::Separator();
      if (galaxySelectedSystemId == 0) {
        ImGui::Text("Select a system on the map.");
      } else {
        const auto& selSys = universe.getSystem(galaxySelectedSystemId);
        const double distLy = (selSys.stub.posLy - currentSystem->stub.posLy).length();
        const double jrMaxLy = fsdBaseRangeLy();
        const double jrNowLy = fsdCurrentRangeLy();

        ImGui::Text("Selected: %s", selSys.stub.name.c_str());
        ImGui::Text("Distance: %.1f ly", distLy);
        ImGui::Text("Jump range: %.1f ly max, %.1f ly current-fuel", jrMaxLy, jrNowLy);

        if (galaxySelectedSystemId == currentSystem->stub.id) {
          ImGui::TextDisabled("This is your current system.");
        } else {
          const double fuelCost = fsdFuelCostFor(distLy);
          ImGui::Text("Direct jump: %s (fuel cost %.1f)", (distLy <= jrMaxLy + 1e-6) ? "IN RANGE" : "OUT OF RANGE", fuelCost);

          if (ImGui::Button("Plot route")) {
            navRoute = plotRouteAStarHops(nearby, currentSystem->stub.id, galaxySelectedSystemId, jrMaxLy);
            navRouteHop = 0;
            navAutoRun = false;
            if (navRoute.empty()) {
              toast(toasts, "No route found. Try increasing the scan radius.", 3.0);
            } else {
              toast(toasts, "Route plotted: " + std::to_string((int)navRoute.size() - 1) + " jumps.", 2.5);
            }
          }
          ImGui::SameLine();
          if (ImGui::Button("Clear route")) {
            navRoute.clear();
            navRouteHop = 0;
            navAutoRun = false;
          }

          ImGui::Checkbox("Auto-run route", &navAutoRun);

          if (!navRoute.empty()) {
            const int totalJumps = (int)navRoute.size() - 1;
            const int remaining = std::max(0, totalJumps - navRouteHop);
            double totalDist = 0.0;
            double totalFuel = 0.0;
            for (std::size_t i = 0; i + 1 < navRoute.size(); ++i) {
              const auto& a = universe.getSystem(navRoute[i]).stub;
              const auto& b = universe.getSystem(navRoute[i + 1]).stub;
              const double d = (b.posLy - a.posLy).length();
              totalDist += d;
              totalFuel += fsdFuelCostFor(d);
            }
            ImGui::Text("Route: %d jumps (%d remaining)", totalJumps, remaining);
            ImGui::Text("Total: %.1f ly, est fuel %.1f", totalDist, totalFuel);
            if (navRouteHop + 1 < navRoute.size()) {
              const auto& nextSys = universe.getSystem(navRoute[navRouteHop + 1]);
              const double hopDist = (nextSys.stub.posLy - currentSystem->stub.posLy).length();
              ImGui::Text("Next hop: %s (%.1f ly)", nextSys.stub.name.c_str(), hopDist);
            }
          }

          if (ImGui::Button("Engage FSD (J)")) {
            if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
              startFsdJumpTo(navRoute[navRouteHop + 1]);
            } else {
              startFsdJumpTo(galaxySelectedSystemId);
            }
          }
        }
      }

      ImGui::TextDisabled("TAB toggles this window, F1 Flight, F2 Market, F3 Contacts");

      ImGui::End();
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    SDL_GL_SwapWindow(window);
  }

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  SDL_GL_DeleteContext(glContext);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
