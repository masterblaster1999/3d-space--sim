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
#include <cmath>
#include <cstring>
#include <optional>
#include <string>
#include <unordered_map>
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
  std::string name;
  sim::Ship ship{};
  bool pirate{false};

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
  ship.setMaxLinearAccelKmS2(0.08);
  ship.setMaxAngularAccelRadS2(1.2);

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
  int selectedStationIndex = 0;

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
  bool showContacts = true;

  // Flight assistance
  bool autopilot = false;

  // Target
  Target target{};

  std::vector<ToastMsg> toasts;

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

  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

  SDL_SetRelativeMouseMode(SDL_FALSE);

  while (running) {
    // Timing
    auto now = std::chrono::high_resolution_clock::now();
    const double dtReal = std::chrono::duration<double>(now - last).count();
    last = now;

    // Events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      ImGui_ImplSDL2_ProcessEvent(&event);

      if (event.type == SDL_QUIT) running = false;
      if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE) running = false;

      if (event.type == SDL_KEYDOWN && !event.key.repeat) {
        if (event.key.keysym.sym == SDLK_ESCAPE) running = false;

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

            toast(toasts, "Loaded " + savePath, 2.5);
          }
        }

        if (event.key.keysym.sym == SDLK_TAB) showGalaxy = !showGalaxy;
        if (event.key.keysym.sym == SDLK_F1) showShip = !showShip;
        if (event.key.keysym.sym == SDLK_F2) showEconomy = !showEconomy;
        if (event.key.keysym.sym == SDLK_F3) showContacts = !showContacts;

        if (event.key.keysym.sym == SDLK_SPACE) paused = !paused;

        if (event.key.keysym.sym == SDLK_p) autopilot = !autopilot;
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

      if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT) {
        if (!io.WantCaptureMouse) {
          // Fire laser
          if (playerLaserCooldown <= 0.0 && !docked) {
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
                hit.alive = false;
                toast(toasts, "Target destroyed. +500 cr", 2.5);
                credits += 500.0;
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

    // Sim step
    const double dtSim = dtReal * timeScale;
    if (!paused) {
      if (!docked) {
        ship.step(dtSim, input);
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

      // Spawn pirates occasionally (simple "threat" loop).
      if (timeDays >= nextPirateSpawnDays && contacts.size() < 8) {
        nextPirateSpawnDays = timeDays + (rng.range(120.0, 220.0) / 86400.0); // every ~2-4 minutes

        Contact p{};
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

      timeDays += dtSim / 86400.0;

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

      if (docked) {
        ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.5f, 1.0f), "DOCKED");
      } else {
        ImGui::TextDisabled("Not docked");
      }

      ImGui::Checkbox("Autopilot (P)", &autopilot);

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
        ImGui::TextDisabled("Target: (none)  [T cycles stations]");
      }

      ImGui::Separator();
      ImGui::Text("Pos (km):   [%.1f %.1f %.1f]", pos.x, pos.y, pos.z);
      ImGui::Text("Vel (km/s): [%.3f %.3f %.3f] |v|=%.3f", vel.x, vel.y, vel.z, vel.length());
      ImGui::Text("AngVel (rad/s): [%.3f %.3f %.3f]", wv.x, wv.y, wv.z);

      ImGui::TextDisabled("Controls:");
      ImGui::BulletText("Translate: WASD + R/F (up/down)");
      ImGui::BulletText("Rotate: Arrow keys + Q/E roll");
      ImGui::BulletText("Boost: LShift   Brake: X");
      ImGui::BulletText("Dampers: Z (on) / C (off)");
      ImGui::BulletText("Target: T (cycle stations), Y (clear)");
      ImGui::BulletText("Docking: L (request clearance), G (dock/undock)");
      ImGui::BulletText("Fire: Left Mouse");
      ImGui::BulletText("Pause: Space   Save: F5   Load: F9");

      ImGui::End();
    }

    if (showEconomy) {
      beginStationSelectorHUD(*currentSystem, selectedStationIndex, docked, dockedStationId);

      if (!currentSystem->stations.empty()) {
        const auto& station = currentSystem->stations[(std::size_t)selectedStationIndex];
        auto& stEcon = universe.stationEconomy(station, timeDays);

        ImGui::Begin("Market Details");
        ImGui::Text("Credits: %.2f", credits);

        const bool canTrade = docked && (station.id == dockedStationId);
        if (!canTrade) {
          ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.35f, 1.0f), "Trade disabled: dock at this station to buy/sell.");
        }

        // Simple docked services
        if (canTrade) {
          if (ImGui::Button("Repair hull (500 cr)")) {
            if (credits >= 500.0) {
              credits -= 500.0;
              playerHull = 100.0;
              toast(toasts, "Ship repaired.", 2.0);
            } else {
              toast(toasts, "Not enough credits.", 2.0);
            }
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
              auto tr = econ::buy(stEcon, station.economyModel, cid, qty[i], credits, 0.10, station.feeRate);
              if (tr.ok) cargo[i] += qty[i];
            }

            ImGui::SameLine();
            if (ImGui::SmallButton("Sell")) {
              const double sellUnits = std::min<double>(qty[i], cargo[i]);
              if (sellUnits > 0.0) {
                auto tr = econ::sell(stEcon, station.economyModel, cid, sellUnits, credits, 0.10, station.feeRate);
                if (tr.ok) cargo[i] -= sellUnits;
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
        ImGui::Text("%s %s", c.name.c_str(), c.pirate ? "[PIRATE]" : "");
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
      static sim::SystemId selected = 0;
      for (const auto& s : nearby) {
        const ImVec2 p = toPx(s.posLy);
        const bool isCurrent = (s.id == currentSystem->stub.id);
        const bool isSel = (s.id == selected);

        ImU32 col = isCurrent ? IM_COL32(255, 240, 160, 255) : IM_COL32(170, 170, 190, 255);
        if (s.factionId != 0) col = IM_COL32(160, 220, 170, 255);
        if (isSel) col = IM_COL32(255, 120, 120, 255);

        draw->AddCircleFilled(p, isCurrent ? 5.5f : 4.0f, col);

        // Click detection
        const float rClick = 6.0f;
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
          const ImVec2 mp = ImGui::GetIO().MousePos;
          const float dx = mp.x - p.x;
          const float dy = mp.y - p.y;
          if (dx*dx + dy*dy <= rClick*rClick) selected = s.id;
        }
      }

      ImGui::EndChild();

      if (selected != 0 && selected != currentSystem->stub.id) {
        if (ImGui::Button("Jump to selected system")) {
          auto it = std::find_if(nearby.begin(), nearby.end(), [&](const sim::SystemStub& s){ return s.id == selected; });
          if (it != nearby.end()) {
            currentStub = *it;
            const auto& sys = universe.getSystem(currentStub.id, &currentStub);
            currentSystem = &sys;

            // Clear transient state on jump
            contacts.clear();
            clearances.clear();
            docked = false;
            dockedStationId = 0;
            autopilot = false;
            target = Target{};
            nextPirateSpawnDays = timeDays + 0.02;

            // spawn near first station
            respawnNearStation(*currentSystem, 0);

            toast(toasts, "Arrived in " + currentSystem->stub.name, 2.5);
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
