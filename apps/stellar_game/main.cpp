#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/core/Hash.h"
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
#include <cstdio>
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
static const char* planetTypeName(sim::PlanetType t) {
  switch (t) {
    case sim::PlanetType::Rocky:    return "Rocky";
    case sim::PlanetType::Desert:   return "Desert";
    case sim::PlanetType::Ocean:    return "Ocean";
    case sim::PlanetType::Ice:      return "Ice";
    case sim::PlanetType::GasGiant: return "Gas Giant";
    default: return "Unknown";
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

static math::Vec3d planetPosKm(const sim::Planet& p, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, timeDays);
  return posAU * kAU_KM;
}

static math::Vec3d planetVelKmS(const sim::Planet& p, double timeDays) {
  // Numeric derivative over a short window.
  const double epsDays = 0.001; // 86.4s
  const math::Vec3d a = planetPosKm(p, timeDays - epsDays);
  const math::Vec3d b = planetPosKm(p, timeDays + epsDays);
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

static bool segmentHitsSphere(const math::Vec3d& aKm,
                             const math::Vec3d& bKm,
                             const math::Vec3d& centerKm,
                             double radiusKm) {
  const math::Vec3d ab = bKm - aKm;
  const double abLenSq = ab.lengthSq();
  if (abLenSq < 1e-12) {
    return (aKm - centerKm).lengthSq() <= radiusKm * radiusKm;
  }

  const double t = std::clamp(math::dot(centerKm - aKm, ab) / abLenSq, 0.0, 1.0);
  const math::Vec3d closest = aKm + ab * t;
  return (closest - centerKm).lengthSq() <= radiusKm * radiusKm;
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

enum class TargetKind : int { None=0, Station=1, Planet=2, Contact=3, Star=4 };

struct Target {
  TargetKind kind{TargetKind::None};
  std::size_t index{0}; // station/planet/contact index
};

enum class ContactRole : int { Pirate=0, Trader=1, Police=2 };

static const char* contactRoleName(ContactRole r) {
  switch (r) {
    case ContactRole::Pirate: return "Pirate";
    case ContactRole::Trader: return "Trader";
    case ContactRole::Police: return "Police";
    default: return "?";
  }
}

struct Contact {
  core::u64 id{0};
  std::string name;
  sim::Ship ship{};

  ContactRole role{ContactRole::Pirate};
  core::u32 factionId{0}; // police/trader affiliation; 0 for pirates / independent
  std::size_t homeStationIndex{0}; // best-effort "patrol" anchor (index into system stations)

  bool missionTarget{false};

  // Behavior flags
  bool hostileToPlayer{false}; // police become hostile when you're wanted / shoot them
  double fleeUntilDays{0.0};   // traders flee after being attacked

  // Combat stats
  double shield{60.0};
  double hull{70.0};

  // Traders can have an approximate loot value (paid on destruction for now).
  double cargoValueCr{0.0};

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
  // Note: Docking is only available in Dear ImGui's "docking" branch.
  // The project pins ImGui's mainline tag by default, so keep docking optional.
  // (If you later switch to the docking branch, feel free to re-enable this.)
  // io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
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
  double playerCannonCooldown = 0.0; // seconds

  // Time
  double timeDays = 0.0;
  double timeScale = 60.0; // simulated seconds per real second
  bool paused = false;

  // Economy
  double credits = 2500.0;
  double explorationDataCr = 0.0; // sellable scan data
  std::array<double, econ::kCommodityCount> cargo{};
  double cargoCapacityKg = 420.0;
  int selectedStationIndex = 0;

  // Ship meta / progression
  double fuel = 45.0;
  double fuelMax = 45.0;
  double fsdRangeLy = 18.0;
  double fsdReadyDay = 0.0;

  // Heat (0..100). Kept intentionally simple for early gameplay feedback.
  double heat = 0.0;

  // Missions + reputation
  core::u64 nextMissionId = 1;
  std::vector<sim::Mission> missions;
  std::unordered_map<core::u32, double> repByFaction;

  // Law / crime (per-faction bounties)
  std::unordered_map<core::u32, double> bountyByFaction;
  double policeAlertUntilDays = 0.0;

  // Exploration / discovery
  std::unordered_set<core::u64> scannedKeys; // scanned bodies/stations in the universe (player-local)

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
  double nextTraderSpawnDays = 0.008;
  double nextPoliceSpawnDays = 0.006;

  // Beams (for laser visuals)
  struct Beam { math::Vec3d aU, bU; float r,g,b; double ttl; };
  std::vector<Beam> beams;

  // Projectiles (kinetic cannons / slugs)
  struct Projectile {
    math::Vec3d prevKm;
    math::Vec3d posKm;
    math::Vec3d velKmS;
    float r{1}, g{1}, b{1};
    double ttlSim{0.0};      // simulated seconds remaining
    double radiusKm{450.0};  // collision radius
    double dmg{0.0};
    bool fromPlayer{false};
    core::u64 shooterId{0};
  };
  std::vector<Projectile> projectiles;

  // Save/load
  const std::string savePath = "savegame.txt";

  // UI state
  bool showGalaxy = true;
  bool showShip = true;
  bool showEconomy = true;
  bool showMissions = true;
  bool showContacts = true;
  bool showScanner = true;

  // Optional mouse steering (relative mouse mode). Toggle with M.
  bool mouseSteer = false;
  float mouseSensitivity = 0.0025f; // torque intent per pixel
  bool mouseInvertY = false;

  // Flight assistance
  bool autopilot = false;
  int autopilotPhase = 0; // 0=staging, 1=corridor

  // Local reference frame ("space is local" feel near moving bodies)
  bool localFrameEnabled = true;
  math::Vec3d localFrameVelKmS{0,0,0};
  double localFrameBlendTauSec = 1.0; // seconds (real time)

  // "Blue zone" turn assist (max turn rate near ~mid-speed when flight assist is on)
  bool blueZoneTurnAssist = true;

  // Supercruise (Elite-style in-system travel)
  enum class SupercruiseState { Idle, Charging, Active, Cooldown };
  SupercruiseState supercruiseState = SupercruiseState::Idle;

  bool supercruiseAssist = true; // "Nav assist" keeps a safe drop profile
  double supercruiseMaxSpeedKmS = 18000.0;

  double supercruiseChargeRemainingSec = 0.0;   // real seconds
  double supercruiseCooldownRemainingSec = 0.0; // real seconds
  double supercruiseCooldownTotalSec = 0.0;     // real seconds (for UI)
  bool supercruiseDropRequested = false;
  bool supercruiseSafeDropReady = false;
  double supercruiseTtaSec = 0.0;
  double supercruiseDistKm = 0.0;
  double supercruiseClosingKmS = 0.0;

  // FSD / hyperspace (system-to-system)
  enum class FsdState { Idle, Charging, Jumping };
  FsdState fsdState = FsdState::Idle;
  sim::SystemId fsdTargetSystem{0};
  double fsdChargeRemainingSec = 0.0;
  double fsdTravelRemainingSec = 0.0;
  double fsdTravelTotalSec = 0.0;
  double fsdFuelCost = 0.0;
  double fsdJumpDistanceLy = 0.0;

  // Galaxy navigation / route plotting
  sim::SystemId galaxySelectedSystemId = 0;
  std::vector<sim::SystemId> navRoute;
  std::size_t navRouteHop = 0;
  bool navAutoRun = false;

  // Scanner interaction (bounty scan + exploration scans)
  bool scanning = false;
  Target scanLockedTarget{};
  core::u64 scanLockedId = 0;
  std::string scanLabel;
  double scanProgressSec = 0.0;
  double scanDurationSec = 4.0;
  double scanRangeKm = 80000.0;

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

auto getBounty = [&](core::u32 factionId) -> double {
  auto it = bountyByFaction.find(factionId);
  return it == bountyByFaction.end() ? 0.0 : std::max(0.0, it->second);
};

auto addBounty = [&](core::u32 factionId, double deltaCr) {
  if (factionId == 0) return;
  bountyByFaction[factionId] = std::max(0.0, getBounty(factionId) + deltaCr);
};

auto clearBounty = [&](core::u32 factionId) {
  if (factionId == 0) return;
  bountyByFaction[factionId] = 0.0;
};

auto commitCrime = [&](core::u32 factionId, double bountyAddCr, double repPenalty, const std::string& reason) {
  if (factionId == 0) return;
  addBounty(factionId, bountyAddCr);
  addRep(factionId, repPenalty);
  policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (120.0 / 86400.0)); // 2 minutes
  nextPoliceSpawnDays = std::min(nextPoliceSpawnDays, timeDays + (6.0 / 86400.0));
  toast(toasts, "Crime (" + factionName(factionId) + "): " + reason, 3.0);
};

// Scan/discovery keys (player-local)
const auto scanKeyStar = [&](sim::SystemId sysId) -> core::u64 {
  return core::hashCombine((core::u64)sysId, 0x53544152ULL); // 'STAR'
};
const auto scanKeyPlanet = [&](sim::SystemId sysId, std::size_t planetIndex) -> core::u64 {
  return core::hashCombine((core::u64)sysId, core::hashCombine(0x504C414EULL, (core::u64)planetIndex)); // 'PLAN'
};
const auto scanKeyStation = [&](sim::StationId stId) -> core::u64 {
  return core::hashCombine((core::u64)stId, 0x53544154ULL); // 'STAT'
};
const auto scanKeySystemComplete = [&](sim::SystemId sysId) -> core::u64 {
  return core::hashCombine((core::u64)sysId, 0x434F4D50ULL); // 'COMP'
};

  auto effectiveFeeRate = [&](const sim::Station& st) -> double {
    return applyRepToFee(st.feeRate, getRep(st.factionId));
  };

  // FSD / jump parameters
  const double kFsdFuelBase = 2.0;
  const double kFsdFuelPerLy = 0.5;
  const double kFsdChargeSec = 4.0;
  const double kFsdCooldownSec = 25.0;

  // Supercruise parameters (in-system travel mode)
  const double kSupercruiseChargeSec = 3.0;
  const double kSupercruiseCooldownSec = 6.0;
  const double kSupercruiseEmergencyCooldownSec = 14.0;
  const double kSupercruiseSafeTtaSec = 7.0; // the classic "7-second rule"

  auto fsdBaseRangeLy = [&]() -> double {
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
    return std::min(fsdBaseRangeLy(), fsdFuelLimitedRangeLy());
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
    if (supercruiseState != SupercruiseState::Idle) {
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

    const double rangeLy = fsdBaseRangeLy();
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
    fsdTravelTotalSec = 0.0;
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

  while (running) {
    // Timing
    auto now = std::chrono::high_resolution_clock::now();
    const double dtReal = std::chrono::duration<double>(now - last).count();
    last = now;

    // Keep relative mouse mode in sync with our control mode.
    // (Automatically release the mouse when docked so UI interaction is painless.)
    if (docked && mouseSteer) {
      mouseSteer = false;
      SDL_SetRelativeMouseMode(SDL_FALSE);
      SDL_GetRelativeMouseState(nullptr, nullptr); // flush deltas
      toast(toasts, "Mouse steer disabled while docked.", 2.0);
    }

    const SDL_bool wantRelMouse = mouseSteer ? SDL_TRUE : SDL_FALSE;
    if (SDL_GetRelativeMouseMode() != wantRelMouse) {
      SDL_SetRelativeMouseMode(wantRelMouse);
      SDL_GetRelativeMouseState(nullptr, nullptr); // flush deltas
    }


    // Per-frame helper: apply damage from the player to a contact (laser / cannon / projectiles).
    auto playerDamageContact = [&](int idx, double dmg) {
      if (idx < 0 || idx >= (int)contacts.size()) return;
      auto& hit = contacts[(std::size_t)idx];
      if (!hit.alive) return;

      // Apply damage
      applyDamage(dmg, hit.shield, hit.hull);

      // Crimes / reactions (on hit)
      if (hit.alive && hit.role == ContactRole::Trader) {
        hit.fleeUntilDays = timeDays + (180.0 / 86400.0); // flee ~3 minutes
        commitCrime(hit.factionId, 250.0, -5.0, "Assault on trader");
      }
      if (hit.alive && hit.role == ContactRole::Police) {
        hit.hostileToPlayer = true;
        commitCrime(hit.factionId, 600.0, -10.0, "Assault on security");
      }

      // If destroyed
      if (hit.hull <= 0.0) {
        const core::u64 deadId = hit.id;
        const ContactRole deadRole = hit.role;
        const core::u32 deadFaction = hit.factionId;
        const double lootCr = hit.cargoValueCr;

        hit.alive = false;

        if (deadRole == ContactRole::Pirate) {
          const double bountyCr = 450.0;
          credits += bountyCr;
          toast(toasts, "Pirate destroyed. +" + std::to_string((int)bountyCr) + " cr", 2.5);

          const core::u32 lf = currentSystem ? currentSystem->stub.factionId : 0;
          if (lf != 0) addRep(lf, +0.5);
        } else if (deadRole == ContactRole::Trader) {
          // Piracy payout (simplified as immediate credits)
          const double payout = std::max(0.0, lootCr);
          credits += payout;
          toast(toasts, "Trader destroyed. Loot +" + std::to_string((int)payout) + " cr (WANTED!)", 3.0);

          commitCrime(deadFaction, 1200.0, -18.0, "Murder of trader");
        } else if (deadRole == ContactRole::Police) {
          toast(toasts, "Security destroyed. (WANTED!)", 3.0);
          commitCrime(deadFaction, 2500.0, -35.0, "Murder of security");
        }

        // Bounty kill missions (pirate targets)
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
    };

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
// Exploration / law
s.explorationDataCr = explorationDataCr;
s.scannedKeys.assign(scannedKeys.begin(), scannedKeys.end());
s.bounties.clear();
for (const auto& [fid, b] : bountyByFaction) {
  if (b > 0.0) s.bounties.push_back({fid, b});
}


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

// Exploration / law
explorationDataCr = s.explorationDataCr;
scannedKeys.clear();
for (core::u64 k : s.scannedKeys) scannedKeys.insert(k);

bountyByFaction.clear();
for (const auto& b : s.bounties) bountyByFaction[b.factionId] = b.bountyCr;
policeAlertUntilDays = 0.0;

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
          nextPirateSpawnDays = timeDays + (rng.range(25.0, 55.0) / 86400.0);
          nextTraderSpawnDays = timeDays + (rng.range(15.0, 35.0) / 86400.0);
          nextPoliceSpawnDays = timeDays + (rng.range(10.0, 25.0) / 86400.0);
            autopilot = false;
            supercruiseState = SupercruiseState::Idle;
            supercruiseChargeRemainingSec = 0.0;
            supercruiseCooldownRemainingSec = 0.0;
      supercruiseCooldownTotalSec = 0.0;
            supercruiseDropRequested = false;
            supercruiseSafeDropReady = false;
            supercruiseTtaSec = 0.0;
            supercruiseDistKm = 0.0;
            supercruiseClosingKmS = 0.0;
            fsdState = FsdState::Idle;
            fsdTargetSystem = 0;
            navRoute.clear();
            navRouteHop = 0;
            navAutoRun = false;
            scanning = false;
            scanProgressSec = 0.0;
            scanLockedId = 0;
            scanLabel.clear();
            scanLockedTarget = Target{};
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
        if (event.key.keysym.sym == SDLK_F6) showScanner = !showScanner;

        if (event.key.keysym.sym == SDLK_SPACE) paused = !paused;

        if (event.key.keysym.sym == SDLK_p) {
          autopilot = !autopilot;
          autopilotPhase = 0;
        }

        if (event.key.keysym.sym == SDLK_m) {
          if (!io.WantCaptureKeyboard) {
            mouseSteer = !mouseSteer;
            SDL_SetRelativeMouseMode(mouseSteer ? SDL_TRUE : SDL_FALSE);
            SDL_GetRelativeMouseState(nullptr, nullptr); // flush deltas
            toast(toasts, std::string("Mouse steer ") + (mouseSteer ? "ON" : "OFF") + " (M)", 1.8);
          }
        }

        if (event.key.keysym.sym == SDLK_h) {
          if (!io.WantCaptureKeyboard && !docked && fsdState == FsdState::Idle) {
            // Interdiction-lite: if pirates are nearby, block *engage* and force emergency drops while active.
            double nearestPirateKm = 1e99;
            for (const auto& c : contacts) {
              if (!c.alive || c.role != ContactRole::Pirate) continue;
              nearestPirateKm = std::min(nearestPirateKm, (c.ship.positionKm() - ship.positionKm()).length());
            }
            const bool pirateNearby = (nearestPirateKm < 120000.0);

            if (supercruiseState == SupercruiseState::Idle) {
              if (target.kind == TargetKind::None) {
                toast(toasts, "No nav target set.", 1.8);
              } else if (pirateNearby) {
                toast(toasts, "Interdiction risk: can't engage supercruise with pirates nearby.", 2.3);
              } else {
                supercruiseState = SupercruiseState::Charging;
                supercruiseChargeRemainingSec = kSupercruiseChargeSec;
                supercruiseDropRequested = false;
                autopilot = false;
                autopilotPhase = 0;
                scanning = false;
                scanProgressSec = 0.0;
                navAutoRun = false;
                toast(toasts, "Supercruise charging...", 1.6);
              }
            } else if (supercruiseState == SupercruiseState::Charging) {
              supercruiseState = SupercruiseState::Idle;
              supercruiseChargeRemainingSec = 0.0;
              toast(toasts, "Supercruise canceled.", 1.4);
            } else if (supercruiseState == SupercruiseState::Active) {
              supercruiseDropRequested = true;
              toast(toasts, "Drop requested.", 1.2);
            } else if (supercruiseState == SupercruiseState::Cooldown) {
              toast(toasts, "Supercruise cooling down...", 1.6);
            }
          }
        }

        if (event.key.keysym.sym == SDLK_k) {
            // Scanner action (mission scans + exploration scans)
            if (scanning) {
              scanning = false;
              scanProgressSec = 0.0;
              scanLockedId = 0;
              scanLabel.clear();
              toast(toasts, "Scan cancelled.", 1.5);
            } else {
              if (docked) {
                toast(toasts, "Undock to scan.", 1.8);
              } else if (supercruiseState != SupercruiseState::Idle || fsdState != FsdState::Idle) {
                toast(toasts, "Scanning unavailable in supercruise / hyperspace.", 2.0);
              } else if (target.kind == TargetKind::None) {
                toast(toasts, "No scan target selected. (Use T/B/N/U or System Scanner)", 2.5);
              } else {
                // Lock scan to current target so cycling targets cancels the scan.
                scanLockedTarget = target;
                scanLockedId = 0;
                scanLabel.clear();
                scanProgressSec = 0.0;

                bool ok = false;

                if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
                  const auto& c = contacts[target.index];
                  if (c.alive) {
                    scanLockedId = c.id;
                    scanDurationSec = 4.0;
                    scanRangeKm = 85000.0;
                    scanLabel = "Contact scan: " + c.name;
                    ok = true;
                  }
                } else if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
                  const auto& st = currentSystem->stations[target.index];
                  scanLockedId = st.id;
                  scanDurationSec = 3.0;
                  scanRangeKm = std::max(25000.0, st.commsRangeKm * 0.9);
                  scanLabel = "Station scan: " + st.name;
                  ok = true;
                } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
                  const auto& p = currentSystem->planets[target.index];
                  scanLockedId = (core::u64)target.index;
                  scanDurationSec = 5.0;
                  const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
                  scanRangeKm = std::max(200000.0, rKm * 45.0);
                  scanLabel = "Planet scan: " + p.name;
                  ok = true;
                } else if (target.kind == TargetKind::Star) {
                  scanLockedId = currentSystem->stub.id;
                  scanDurationSec = 2.5;
                  scanRangeKm = 1.0e18; // anywhere in-system for now
                  scanLabel = std::string("Star scan: ") + starClassName(currentSystem->star.cls);
                  ok = true;
                }

                if (ok) {
                  scanning = true;
                  toast(toasts, scanLabel + " (hold steady)...", 2.0);
                } else {
                  toast(toasts, "Invalid scan target.", 2.0);
                }
              }
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
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
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
  // cycle planet targets
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
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
  // cycle contact targets (alive only)
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  if (!contacts.empty()) {
    std::size_t start = 0;
    if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
      start = (target.index + 1) % contacts.size();
    }
    std::size_t idx = start;
    for (std::size_t step = 0; step < contacts.size(); ++step) {
      if (contacts[idx].alive) {
        target.kind = TargetKind::Contact;
        target.index = idx;
        break;
      }
      idx = (idx + 1) % contacts.size();
    }
  }
}

if (event.key.keysym.sym == SDLK_u) {
  // target system primary star (for scanning / nav reference)
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  target.kind = TargetKind::Star;
  target.index = 0;
}

if (event.key.keysym.sym == SDLK_y) {
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  target = Target{};
}

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
          if (supercruiseState != SupercruiseState::Idle) {
            toast(toasts, "Cannot dock while in supercruise.", 2.0);
          } else if (docked) {
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
          if (playerLaserCooldown <= 0.0 && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
            playerLaserCooldown = 0.18; // rate of fire
            heat = std::min(120.0, heat + 2.5);
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

  // Apply damage
  applyDamage(dmg, hit.shield, hit.hull);

  // Crimes / reactions (on hit)
  if (hit.alive && hit.role == ContactRole::Trader) {
    hit.fleeUntilDays = timeDays + (180.0 / 86400.0); // flee ~3 minutes
    commitCrime(hit.factionId, 250.0, -5.0, "Assault on trader");
  }
  if (hit.alive && hit.role == ContactRole::Police) {
    hit.hostileToPlayer = true;
    commitCrime(hit.factionId, 600.0, -10.0, "Assault on security");
  }

  // If destroyed
  if (hit.hull <= 0.0) {
    const core::u64 deadId = hit.id;
    const ContactRole deadRole = hit.role;
    const core::u32 deadFaction = hit.factionId;
    const double lootCr = hit.cargoValueCr;

    hit.alive = false;

    if (deadRole == ContactRole::Pirate) {
      const double bountyCr = 450.0;
      credits += bountyCr;
      toast(toasts, "Pirate destroyed. +" + std::to_string((int)bountyCr) + " cr", 2.5);

      const core::u32 lf = currentSystem ? currentSystem->stub.factionId : 0;
      if (lf != 0) addRep(lf, +0.5);
    } else if (deadRole == ContactRole::Trader) {
      // Piracy payout (simplified as immediate credits)
      const double payout = std::max(0.0, lootCr);
      credits += payout;
      toast(toasts, "Trader destroyed. Loot +" + std::to_string((int)payout) + " cr (WANTED!)", 3.0);

      commitCrime(deadFaction, 1200.0, -18.0, "Murder of trader");
    } else if (deadRole == ContactRole::Police) {
      toast(toasts, "Security destroyed. (WANTED!)", 3.0);
      commitCrime(deadFaction, 2500.0, -35.0, "Murder of security");
    }

    // Bounty kill missions (pirate targets)
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

      if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_RIGHT) {
        if (!io.WantCaptureMouse) {
          // Fire cannon projectile (kinetic)
          if (playerCannonCooldown <= 0.0 && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
            // NOTE: cooldown is in *sim* seconds, so it naturally scales with timeScale.
            playerCannonCooldown = 0.9;
            heat = std::min(120.0, heat + 4.5);

            const double muzzleSpeedKmS = 250.0;
            const double rangeKm = 130000.0;
            const double dmg = 38.0;

            const double ttlSim = rangeKm / muzzleSpeedKmS;

            const math::Vec3d fwd = ship.forward().normalized();
            const math::Vec3d spawnKm = ship.positionKm() + fwd * 400.0; // start slightly ahead of the ship

            Projectile p{};
            p.prevKm = spawnKm;
            p.posKm = spawnKm;
            p.velKmS = ship.velocityKmS() + fwd * muzzleSpeedKmS;
            p.r = 1.0f; p.g = 0.95f; p.b = 0.65f;
            p.ttlSim = ttlSim;
            p.radiusKm = 700.0;
            p.dmg = dmg;
            p.fromPlayer = true;
            p.shooterId = 0;
            projectiles.push_back(p);

            // Tiny recoil impulse
            ship.setVelocityKmS(ship.velocityKmS() - fwd * 0.002);
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

      // Mouse steering (relative mode). Uses local pitch/yaw.
      if (mouseSteer && !io.WantCaptureMouse) {
        int dx = 0, dy = 0;
        SDL_GetRelativeMouseState(&dx, &dy);
        const double sx = (double)mouseSensitivity;
        const double yaw = (double)dx * sx;
        const double pitch = (double)(mouseInvertY ? dy : -dy) * sx;
        input.torqueLocal.y += yaw;
        input.torqueLocal.x += pitch;
      }

      input.boost = keys[SDL_SCANCODE_LSHIFT] != 0;
      input.brake = keys[SDL_SCANCODE_X] != 0;

      static bool dampers = true;
      if (keys[SDL_SCANCODE_Z]) dampers = true;
      if (keys[SDL_SCANCODE_C]) dampers = false;
      input.dampers = dampers;
    }

    // Autopilot: station approach assist (staging + docking corridor guidance).
    if (autopilot && !docked && currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size() && !captureKeys) {
      const auto& st = currentSystem->stations[target.index];

      const math::Vec3d stPos = stationPosKm(st, timeDays);
      const math::Quatd stQ = stationOrient(st, stPos, timeDays);
      const math::Vec3d stV = stationVelKmS(st, timeDays);

      const double zEntrance = st.radiusKm * 1.10;
      const double zCorridorEnd = zEntrance + st.approachLengthKm;
      const double zStage = zCorridorEnd + st.radiusKm * 0.8;

      const math::Vec3d relLocal = stQ.conjugate().rotate(ship.positionKm() - stPos);
      const double lateralDist = math::Vec3d{relLocal.x, relLocal.y, 0.0}.length();

      // Phase switching: start by staging at the corridor end, then commit through the corridor.
      if (autopilotPhase == 0) {
        if (relLocal.z < zCorridorEnd * 1.05 && lateralDist < st.approachRadiusKm * 0.65) autopilotPhase = 1;
      } else {
        // If we drift far off-axis, go back to staging.
        if (relLocal.z > zStage || lateralDist > st.approachRadiusKm * 1.25) autopilotPhase = 0;
      }

      math::Vec3d desiredLocal{0,0,0};
      if (autopilotPhase == 0) desiredLocal = {0,0, zStage};
      else desiredLocal = {0,0, zEntrance + st.radiusKm * 0.7};

      const math::Vec3d desiredPoint = stPos + stQ.rotate(desiredLocal);
      const math::Vec3d rel = desiredPoint - ship.positionKm();
      const double dist = rel.length();
      const math::Vec3d dir = (dist > 1e-6) ? (rel / dist) : math::Vec3d{0,0,0};

      const double maxV = (autopilotPhase == 0) ? std::max(st.maxApproachSpeedKmS * 2.5, 0.45)
                                                : (st.maxApproachSpeedKmS * 0.85);
      double vMag = std::min(maxV, 0.004 * dist);

      // Slow down more if we're not centered when committing.
      if (autopilotPhase == 1) {
        const double frac = std::clamp(lateralDist / std::max(1.0, st.approachRadiusKm), 0.0, 1.0);
        vMag *= (1.0 - 0.55 * frac);
      }

      const math::Vec3d desiredVel = stV + dir * vMag;
      const math::Vec3d dv = desiredVel - ship.velocityKmS();

      if (dv.lengthSq() > 1e-12) input.thrustLocal = ship.orientation().conjugate().rotate(dv.normalized());
      else input.thrustLocal = {0,0,0};

      input.dampers = true;

      // Align forward into station (towards -axisOut).
      const math::Vec3d axisOut = stQ.rotate({0,0,1});
      const math::Vec3d desiredFwdWorld = -axisOut;

      const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(desiredFwdWorld);
      const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
      const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

      // Roll: keep ship "up" aligned to station +Y so you fly the slot level.
      const math::Vec3d desiredUpWorld = stQ.rotate({0,1,0});
      const math::Vec3d desiredUpLocal = ship.orientation().conjugate().rotate(desiredUpWorld);
      const double rollErr = std::atan2(desiredUpLocal.x, desiredUpLocal.y);

      input.torqueLocal.x = (float)std::clamp(pitchErr * 1.8, -1.0, 1.0);
      input.torqueLocal.y = (float)std::clamp(yawErr * 1.8, -1.0, 1.0);
      input.torqueLocal.z = (float)std::clamp(rollErr * 1.6, -1.0, 1.0);
    }

    // Supercruise: fast in-system travel assist to current target.
    // Also: update the "local reference frame" velocity that Flight Assist uses for dampers.
    if (localFrameEnabled && currentSystem) {
      math::Vec3d desiredVel{0,0,0};
      double best = 1e99;

      // Prefer stations when nearby (keeps you "riding along" with their orbital motion).
      for (const auto& st : currentSystem->stations) {
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const double d = (stPos - ship.positionKm()).length();
        const double influence = std::max(st.commsRangeKm * 1.5, st.approachLengthKm * 2.0 + st.radiusKm * 12.0);
        if (d < influence && d < best) {
          best = d;
          desiredVel = stationVelKmS(st, timeDays);
        }
      }

      // Fall back to planets if no station "wins".
      for (const auto& p : currentSystem->planets) {
        const math::Vec3d pPos = planetPosKm(p, timeDays);
        const double d = (pPos - ship.positionKm()).length();
        const double rKm = p.radiusEarth * 6371.0;
        const double influence = std::max(250000.0, rKm * 80.0);
        if (d < influence && d < best) {
          best = d;
          desiredVel = planetVelKmS(p, timeDays);
        }
      }

      const double tau = std::max(0.05, localFrameBlendTauSec);
      const double alpha = 1.0 - std::exp(-dtReal / tau);
      localFrameVelKmS = localFrameVelKmS + (desiredVel - localFrameVelKmS) * alpha;
    } else {
      localFrameVelKmS = {0,0,0};
    }
    ship.setDampingFrameVelocityKmS(localFrameVelKmS);

    // Supercruise state timers (real time)
    if (supercruiseState == SupercruiseState::Charging) {
      supercruiseChargeRemainingSec = std::max(0.0, supercruiseChargeRemainingSec - dtReal);
      if (supercruiseChargeRemainingSec <= 0.0) {
        supercruiseState = SupercruiseState::Active;
        supercruiseDropRequested = false;
        toast(toasts, "Supercruise engaged.", 1.6);
      }
    } else if (supercruiseState == SupercruiseState::Cooldown) {
      supercruiseCooldownRemainingSec = std::max(0.0, supercruiseCooldownRemainingSec - dtReal);
      if (supercruiseCooldownRemainingSec <= 0.0) {
        supercruiseState = SupercruiseState::Idle;
      }
    }

    // Reset per-frame HUD values
    supercruiseSafeDropReady = false;
    supercruiseTtaSec = 0.0;
    supercruiseDistKm = 0.0;
    supercruiseClosingKmS = 0.0;

    auto endSupercruise = [&](bool emergency, const math::Vec3d& destVelKmS, const math::Vec3d& dirToDest, const char* msg) {
      supercruiseState = SupercruiseState::Cooldown;
      supercruiseCooldownRemainingSec = emergency ? kSupercruiseEmergencyCooldownSec : kSupercruiseCooldownSec;
      supercruiseCooldownTotalSec = supercruiseCooldownRemainingSec;
      supercruiseDropRequested = false;

      const double approachSpeed = emergency ? 0.45 : 0.16;
      ship.setVelocityKmS(destVelKmS + dirToDest * approachSpeed);

      if (emergency) {
        heat = std::min(120.0, heat + 25.0);
        playerShield = std::max(0.0, playerShield - 12.0);
        playerHull = std::max(0.0, playerHull - 4.0);
        ship.setAngularVelocityRadS({rng.range(-0.6, 0.6), rng.range(-0.6, 0.6), rng.range(-0.4, 0.4)});
      } else {
        ship.setAngularVelocityRadS({0,0,0});
      }

      toast(toasts, msg, 2.2);
    };

    if (supercruiseState == SupercruiseState::Active && !docked && !captureKeys && fsdState == FsdState::Idle) {
      // Determine destination (station or planet)
      bool hasDest = false;
      math::Vec3d destPosKm{0,0,0};
      math::Vec3d destVelKmS{0,0,0};
      double dropKm = 0.0;

      if (currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
        const auto& st = currentSystem->stations[target.index];
        destPosKm = stationPosKm(st, timeDays);
        destVelKmS = stationVelKmS(st, timeDays);
        dropKm = std::max(15000.0, st.radiusKm * 3.0);
        hasDest = true;
      } else if (currentSystem && target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
        const auto& p = currentSystem->planets[target.index];
        destPosKm = planetPosKm(p, timeDays);
        destVelKmS = planetVelKmS(p, timeDays);
        const double rKm = p.radiusEarth * 6371.0;
        dropKm = std::max(60000.0, rKm * 12.0);
        hasDest = true;
      }

      math::Vec3d dirToDest = ship.forward().normalized();
      if (hasDest) {
        const math::Vec3d rel = destPosKm - ship.positionKm();
        const double dist = rel.length();
        if (dist > 1e-6) dirToDest = rel / dist;
      }

      // Interdiction-lite: pirates nearby can force an emergency drop
      double nearestPirateKm = 1e99;
      for (const auto& c : contacts) {
        if (!c.alive || c.role != ContactRole::Pirate) continue;
        nearestPirateKm = std::min(nearestPirateKm, (c.ship.positionKm() - ship.positionKm()).length());
      }
      if (nearestPirateKm < 120000.0) {
        endSupercruise(true, destVelKmS, dirToDest, "Interdicted! Emergency drop.");
      } else if (!hasDest) {
        endSupercruise(true, localFrameVelKmS, dirToDest, "Supercruise lost target. Emergency drop.");
      } else {
        // Compute safe-drop window ("7-second rule")
        const math::Vec3d rel = destPosKm - ship.positionKm();
        const double dist = rel.length();
        supercruiseDistKm = dist;

        if (dist < 1e-6) {
          endSupercruise(false, destVelKmS, dirToDest, "Supercruise drop.");
        } else {
          const math::Vec3d dir = rel / dist;
          const math::Vec3d vRel = ship.velocityKmS() - destVelKmS;
          const double closing = math::dot(vRel, dir);
          supercruiseClosingKmS = closing;

          const double tta = (closing > 1e-3) ? (dist / closing) : 1e9;
          supercruiseTtaSec = tta;

          const double safeMin = kSupercruiseSafeTtaSec - 2.0;
          const double safeMax = kSupercruiseSafeTtaSec + 2.0;
          const bool safeWindow = (dist < dropKm) && (tta > safeMin) && (tta < safeMax) && (closing > 0.05);
          supercruiseSafeDropReady = safeWindow;

          // Manual drop request (H while in supercruise)
          if (supercruiseDropRequested) {
            if (safeWindow) endSupercruise(false, destVelKmS, dir, "Supercruise drop.");
            else endSupercruise(true, destVelKmS, dir, "EMERGENCY DROP!");
          } else if (supercruiseAssist && safeWindow) {
            // Nav assist auto-drop when safe
            endSupercruise(false, destVelKmS, dir, "Supercruise drop.");
          } else {
            // Follow a speed profile toward the target
            const double desiredSpeed =
              supercruiseAssist ? std::clamp(dist / kSupercruiseSafeTtaSec, 60.0, supercruiseMaxSpeedKmS)
                               : std::clamp(dist * 0.0008, 90.0, supercruiseMaxSpeedKmS);

            const math::Vec3d desiredVel = destVelKmS + dir * desiredSpeed;
            const math::Vec3d dv = desiredVel - ship.velocityKmS();

            if (dv.lengthSq() > 1e-12) {
              const math::Vec3d accelDir = dv.normalized();
              input.thrustLocal = ship.orientation().conjugate().rotate(accelDir);
            } else {
              input.thrustLocal = {0,0,0};
            }

            // Always face the direction of travel while in supercruise
            const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(dir);
            const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
            const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

            input.torqueLocal.x = (float)std::clamp(pitchErr * 1.6, -1.0, 1.0);
            input.torqueLocal.y = (float)std::clamp(yawErr * 1.6, -1.0, 1.0);
            input.torqueLocal.z = 0.0f;

            input.dampers = true;
            input.brake = false;
            input.boost = false;

            ship.setMaxLinearAccelKmS2(6.0);
            ship.setMaxAngularAccelRadS2(1.2);
          }
        }
      }
    }

    if (supercruiseState == SupercruiseState::Active) {
      ship.setMaxLinearAccelKmS2(6.0);
      ship.setMaxAngularAccelRadS2(1.2);
    }

    // Restore handling caps when not in active supercruise (also applies "blue zone" turn assist).
    if (supercruiseState != SupercruiseState::Active) {
      ship.setMaxLinearAccelKmS2(kPlayerBaseLinAccelKmS2);

      double ang = kPlayerBaseAngAccelRadS2;
      if (blueZoneTurnAssist && input.dampers && !docked && fsdState == FsdState::Idle) {
        const double v = (ship.velocityKmS() - localFrameVelKmS).length();
        const double vTerm = std::max(0.001, ship.maxLinearAccelKmS2() / std::max(0.05, ship.dampingLinear()));
        const double x = std::clamp(v / vTerm, 0.0, 1.0);

        // Gaussian-ish bump centered around mid-speed.
        const double center = 0.55;
        const double sigma = 0.25;
        const double bump = std::exp(-std::pow((x - center) / sigma, 2.0));
        const double factor = 0.65 + 0.35 * bump;
        ang *= factor;
      }

      ship.setMaxAngularAccelRadS2(ang);
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
          fsdState = FsdState::Jumping;
          fsdTravelTotalSec = 1.25 + fsdJumpDistanceLy * 0.08;
          fsdTravelRemainingSec = fsdTravelTotalSec;
          toast(toasts, "FSD: entering hyperspace...", 1.8);
        }
      } else if (fsdState == FsdState::Jumping) {
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
          scanLockedId = 0;
          scanLabel.clear();
          scanLockedTarget = Target{};
          target = Target{};

          // Spawn near the first station.
          respawnNearStation(*currentSystem, 0);
          galaxySelectedSystemId = currentSystem->stub.id;

          // Cooldown (sim time)
          fsdReadyDay = timeDays + (kFsdCooldownSec / 86400.0);

          // Advance route cursor if this jump matches the plotted route.
          if (!navRoute.empty() && navRouteHop + 1 < navRoute.size() && navRoute[navRouteHop + 1] == currentSystem->stub.id) {
            navRouteHop++;
            if (navRouteHop + 1 >= navRoute.size()) {
              navAutoRun = false;
              toast(toasts, "Route complete.", 2.0);
            }
          }
// QoL: auto-select the next hop in the Galaxy UI after each jump.
if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
  galaxySelectedSystemId = navRoute[navRouteHop + 1];
} else {
  galaxySelectedSystemId = currentSystem->stub.id;
}


          toast(toasts, std::string("Arrived in ") + currentSystem->stub.name + ".", 2.0);

          fsdState = FsdState::Idle;
          fsdTargetSystem = 0;
          fsdJumpDistanceLy = 0.0;
          fsdFuelCost = 0.0;
          fsdTravelTotalSec = 0.0;
          fsdJustArrived = true;
        }
      }

      // Auto-run route (hands-free multi-jump). We only auto-trigger when safe.
      if (navAutoRun && fsdState == FsdState::Idle && timeDays >= fsdReadyDay && !docked && supercruiseState == SupercruiseState::Idle) {
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
        if (supercruiseState == SupercruiseState::Active) {
          const double v = ship.velocityKmS().length();
          fuel -= dtSim * (0.0020 + v * 0.00000015);
        }

        if (fuel < 0.0) fuel = 0.0;
        if (fuel > fuelMax) fuel = fuelMax;

        if (fuel <= 0.0 && supercruiseState == SupercruiseState::Active) {
          supercruiseState = SupercruiseState::Cooldown;
          supercruiseCooldownRemainingSec = kSupercruiseEmergencyCooldownSec;
          supercruiseCooldownTotalSec = supercruiseCooldownRemainingSec;
          supercruiseDropRequested = false;

          heat = std::min(120.0, heat + 25.0);
          playerShield = std::max(0.0, playerShield - 10.0);
          playerHull = std::max(0.0, playerHull - 3.0);
          ship.setVelocityKmS(localFrameVelKmS + ship.forward().normalized() * 0.45);
          ship.setAngularVelocityRadS({rng.range(-0.6, 0.6), rng.range(-0.6, 0.6), rng.range(-0.4, 0.4)});

          toast(toasts, "Fuel depleted: EMERGENCY DROP!", 2.4);
        }
      }

const bool combatSimEnabled = (fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle);
if (combatSimEnabled) {
  // ---- Spawns / combat / traffic ----

  int alivePirates = 0;
int aliveTraders = 0;
int alivePolice = 0;
int aliveTotal = 0;
for (const auto& c : contacts) {
  if (!c.alive) continue;
  ++aliveTotal;
  if (c.role == ContactRole::Pirate) ++alivePirates;
  if (c.role == ContactRole::Trader) ++aliveTraders;
  if (c.role == ContactRole::Police) ++alivePolice;
}

  const core::u32 localFaction = currentSystem ? currentSystem->stub.factionId : 0;
  const double localRep = getRep(localFaction);
  const double localBounty = getBounty(localFaction);
  const bool playerWantedHere = (localFaction != 0) && (localBounty > 0.0);

  // Simple helper for random directions.
  auto randDir = [&]() -> math::Vec3d {
    return math::Vec3d{rng.range(-1.0,1.0), rng.range(-0.3,0.3), rng.range(-1.0,1.0)}.normalized();
  };

  auto chooseHomeStation = [&]() -> std::size_t {
    if (!currentSystem || currentSystem->stations.empty()) return 0;
    return (std::size_t)(rng.nextU64() % (core::u64)currentSystem->stations.size());
  };

  auto spawnNearStation = [&](std::size_t stIdx, double minDistMul, double maxDistMul, sim::Ship& outShip) {
    if (!currentSystem || currentSystem->stations.empty()) {
      // fallback: near player
      const math::Vec3d d = randDir();
      const double distKm = rng.range(60000.0, 130000.0);
      outShip.setPositionKm(ship.positionKm() + d * distKm);
      outShip.setVelocityKmS(ship.velocityKmS());
      outShip.setOrientation(quatFromTo({0,0,1}, (-d).normalized()));
      return;
    }
    stIdx = std::min(stIdx, currentSystem->stations.size() - 1);
    const auto& st = currentSystem->stations[stIdx];
    const math::Vec3d stPos = stationPosKm(st, timeDays);
    const math::Quatd stQ = stationOrient(st, stPos, timeDays);
    const math::Vec3d axis = stQ.rotate({0,0,1});
    const math::Vec3d tangent = stQ.rotate({1,0,0});

    const double distKm = rng.range(st.radiusKm * minDistMul, st.radiusKm * maxDistMul);
    const math::Vec3d p = stPos + axis * distKm + tangent * rng.range(-st.radiusKm*2.0, st.radiusKm*2.0);
    outShip.setPositionKm(p);
    outShip.setVelocityKmS(stationVelKmS(st, timeDays));
    outShip.setOrientation(quatFromTo({0,0,1}, (stPos - p).normalized()));
  };

  auto spawnPirate = [&]() {
    Contact p{};
    p.id = std::max<core::u64>(1, rng.nextU64());
    p.role = ContactRole::Pirate;
    p.name = "Pirate " + std::to_string((int)contacts.size() + 1);
    p.ship.setMaxLinearAccelKmS2(0.06);
    p.ship.setMaxAngularAccelRadS2(0.9);

    // Spawn somewhere near the player, but not too close.
    const math::Vec3d d = randDir();
    const double distKm = rng.range(55000.0, 120000.0);
    p.ship.setPositionKm(ship.positionKm() + d * distKm);
    p.ship.setVelocityKmS(ship.velocityKmS());
    p.ship.setOrientation(quatFromTo({0,0,1}, (-d).normalized()));

    contacts.push_back(std::move(p));
    toast(toasts, "Contact: pirate detected!", 3.0);
  };

  auto spawnTrader = [&]() {
    Contact t{};
    t.id = std::max<core::u64>(1, rng.nextU64());
    t.role = ContactRole::Trader;
    t.factionId = localFaction;
    t.homeStationIndex = chooseHomeStation();
    t.name = "Trader " + std::to_string((int)contacts.size() + 1);
    t.shield = 35.0;
    t.hull = 55.0;
    t.cargoValueCr = rng.range(350.0, 900.0);

    t.ship.setMaxLinearAccelKmS2(0.05);
    t.ship.setMaxAngularAccelRadS2(0.6);
    spawnNearStation(t.homeStationIndex, 18.0, 30.0, t.ship);

    contacts.push_back(std::move(t));
  };

  auto spawnPolice = [&]() {
    Contact p{};
    p.id = std::max<core::u64>(1, rng.nextU64());
    p.role = ContactRole::Police;
    p.factionId = localFaction;
    p.homeStationIndex = chooseHomeStation();
    p.name = "Security " + std::to_string((int)contacts.size() + 1);
    p.shield = 90.0;
    p.hull = 90.0;

    p.ship.setMaxLinearAccelKmS2(0.075);
    p.ship.setMaxAngularAccelRadS2(1.05);
    spawnNearStation(p.homeStationIndex, 16.0, 22.0, p.ship);

    contacts.push_back(std::move(p));
    toast(toasts, "Local security on patrol.", 2.0);
  };

  // Pirates occasionally (baseline threat).
  if (timeDays >= nextPirateSpawnDays && alivePirates < 4 && aliveTotal < 14) {
    nextPirateSpawnDays = timeDays + (rng.range(120.0, 220.0) / 86400.0); // every ~2-4 minutes
    spawnPirate();
    ++alivePirates;
    ++aliveTotal;
  }

  // Traders / traffic: gives you something to pirate (but doing so triggers police).
  if (timeDays >= nextTraderSpawnDays && aliveTraders < 3 && currentSystem && !currentSystem->stations.empty() && aliveTotal < 14) {
    nextTraderSpawnDays = timeDays + (rng.range(70.0, 140.0) / 86400.0);
    spawnTrader();
    ++aliveTraders;
    ++aliveTotal;
  }

  // Police / patrols: scale with threat + your legal status.
  if (localFaction != 0) {
    int desiredPolice = 1;
    if (playerWantedHere) desiredPolice += 2;
    if (alivePirates > 0) desiredPolice += 1;
    if (localRep < -25.0) desiredPolice += 1;
    desiredPolice = std::clamp(desiredPolice, 0, 5);

    // If recently alerted by a crime, tighten the spawn interval.
    const double spawnMinSec = (policeAlertUntilDays > timeDays) ? 12.0 : 55.0;
    const double spawnMaxSec = (policeAlertUntilDays > timeDays) ? 24.0 : 95.0;

    if (timeDays >= nextPoliceSpawnDays && alivePolice < desiredPolice && aliveTotal < 16) {
      nextPoliceSpawnDays = timeDays + (rng.range(spawnMinSec, spawnMaxSec) / 86400.0);
      spawnPolice();
      ++alivePolice;
      ++aliveTotal;
    }
  }

  // Ensure mission bounty targets exist in their target system.
  if (!docked && !missions.empty()) {
    for (const auto& m : missions) {
      if (m.completed || m.failed) continue;
      if (!((m.type == sim::MissionType::BountyScan) || (m.type == sim::MissionType::BountyKill))) continue;
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
      tgt.role = ContactRole::Pirate;
      tgt.missionTarget = true;
      tgt.name = "Bounty Target";
      tgt.shield = 80.0;
      tgt.hull = 90.0;
      tgt.ship.setMaxLinearAccelKmS2(0.07);
      tgt.ship.setMaxAngularAccelRadS2(1.0);

      const math::Vec3d d = randDir();
      const double distKm = rng.range(65000.0, 150000.0);
      tgt.ship.setPositionKm(ship.positionKm() + d * distKm);
      tgt.ship.setVelocityKmS(ship.velocityKmS());
      tgt.ship.setOrientation(quatFromTo({0,0,1}, (-d).normalized()));

      contacts.push_back(std::move(tgt));
      break; // spawn at most one target per frame
    }
  }

  // Find a nearby pirate index (for police).
  auto nearestPirateIndex = [&](const math::Vec3d& fromKm, double maxDistKm) -> std::optional<std::size_t> {
    double bestD = maxDistKm;
    std::optional<std::size_t> best{};
    for (std::size_t i = 0; i < contacts.size(); ++i) {
      const auto& c = contacts[i];
      if (!c.alive || c.role != ContactRole::Pirate) continue;
      const double d = (c.ship.positionKm() - fromKm).length();
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  };

  auto chaseTarget = [&](sim::Ship& selfShip,
                         sim::ShipInput& ai,
                         const math::Vec3d& targetPosKm,
                         const math::Vec3d& targetVelKmS,
                         double desiredDistKm,
                         double maxSpeedKmS,
                         double faceGain) {
    ai.dampers = true;

    const math::Vec3d to = targetPosKm - selfShip.positionKm();
    const double dist = to.length();
    const math::Vec3d toN = (dist > 1e-6) ? (to / dist) : math::Vec3d{0,0,1};

    // Face target
    const math::Vec3d desiredFwdWorld = toN;
    const math::Vec3d desiredFwdLocal = selfShip.orientation().conjugate().rotate(desiredFwdWorld);
    const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
    const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);
    ai.torqueLocal.x = std::clamp(pitchErr * faceGain, -1.0, 1.0);
    ai.torqueLocal.y = std::clamp(yawErr * faceGain, -1.0, 1.0);
    ai.torqueLocal.z = 0.0;

    // Speed control
    double vAim = 0.0;
    if (dist > desiredDistKm) vAim = std::min(maxSpeedKmS, 0.000004 * (dist - desiredDistKm));
    if (dist < desiredDistKm * 0.60) vAim = -0.10;

    const math::Vec3d desiredVel = targetVelKmS + toN * vAim;
    const math::Vec3d dv = desiredVel - selfShip.velocityKmS();
    const math::Vec3d accelWorldDir = (dv.lengthSq() > 1e-9) ? dv.normalized() : math::Vec3d{0,0,0};
    ai.thrustLocal = selfShip.orientation().conjugate().rotate(accelWorldDir);
  };

  // Contacts AI + combat
  for (auto& c : contacts) {
    if (!c.alive) continue;

    // Keep contact dampers in the same local reference frame as the player.
    c.ship.setDampingFrameVelocityKmS(localFrameVelKmS);

    // cooldowns
    c.fireCooldown = std::max(0.0, c.fireCooldown - dtSim);

    // ---- PIRATES ----
    if (c.role == ContactRole::Pirate) {
      sim::ShipInput ai{};
      chaseTarget(c.ship, ai, ship.positionKm(), ship.velocityKmS(), 35000.0, 0.22, 1.8);
      c.ship.step(dtSim, ai);

      // Fire if aligned
      const math::Vec3d to = ship.positionKm() - c.ship.positionKm();
      const double dist = to.length();
      if (c.fireCooldown <= 0.0 && dist < 90000.0) {
        const math::Vec3d toN = (dist > 1e-6) ? (to / dist) : math::Vec3d{0,0,1};
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
      continue;
    }

    // ---- TRADERS ----
    if (c.role == ContactRole::Trader) {
      sim::ShipInput ai{};
      ai.dampers = true;

      const bool fleeing = (timeDays < c.fleeUntilDays);
      if (fleeing) {
        // Run away from player.
        const math::Vec3d away = (c.ship.positionKm() - ship.positionKm());
        const math::Vec3d dir = (away.lengthSq() > 1e-6) ? away.normalized() : math::Vec3d{0,0,1};
        chaseTarget(c.ship, ai, c.ship.positionKm() + dir * 200000.0, ship.velocityKmS(), 0.0, 0.25, 1.4);
      } else if (currentSystem && !currentSystem->stations.empty()) {
        // Lazy orbit/patrol around their home station.
        const auto& st = currentSystem->stations[std::min(c.homeStationIndex, currentSystem->stations.size()-1)];
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const double baseR = st.radiusKm * 26.0;
        const double ang = std::fmod((double)(c.id % 1000) * 0.01 + timeDays * 0.75, 2.0 * math::kPi);
        const math::Vec3d offset = math::Vec3d{std::cos(ang), 0.15*std::sin(ang*0.7), std::sin(ang)} * baseR;
        chaseTarget(c.ship, ai, stPos + offset, stationVelKmS(st, timeDays), baseR * 0.6, 0.18, 1.2);
      }

      c.ship.step(dtSim, ai);
      continue;
    }

    // ---- POLICE ----
    if (c.role == ContactRole::Police) {
      // Police engage if you're wanted here, or if you attacked security.
      const bool hostile = c.hostileToPlayer || (c.factionId != 0 && getBounty(c.factionId) > 0.0) || (c.factionId != 0 && getRep(c.factionId) < -45.0);

      // If not hostile to player, try to engage pirates near the player.
      std::optional<std::size_t> pirateIdx{};
      if (!hostile) {
        pirateIdx = nearestPirateIndex(c.ship.positionKm(), 140000.0);
      }

      sim::ShipInput ai{};
      if (hostile) {
        chaseTarget(c.ship, ai, ship.positionKm(), ship.velocityKmS(), 42000.0, 0.26, 2.0);
      } else if (pirateIdx) {
        const auto& p = contacts[*pirateIdx];
        chaseTarget(c.ship, ai, p.ship.positionKm(), p.ship.velocityKmS(), 42000.0, 0.26, 2.0);
      } else if (currentSystem && !currentSystem->stations.empty()) {
        // Patrol around home station.
        const auto& st = currentSystem->stations[std::min(c.homeStationIndex, currentSystem->stations.size()-1)];
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const double baseR = st.radiusKm * 22.0;
        const double ang = std::fmod((double)(c.id % 1000) * 0.012 + timeDays * 1.10, 2.0 * math::kPi);
        const math::Vec3d offset = math::Vec3d{std::cos(ang), 0.10*std::sin(ang*0.6), std::sin(ang)} * baseR;
        chaseTarget(c.ship, ai, stPos + offset, stationVelKmS(st, timeDays), baseR * 0.6, 0.20, 1.6);
      }

      c.ship.step(dtSim, ai);

      // Fire if aligned (at player OR pirate target).
      const math::Vec3d targetPos = hostile ? ship.positionKm()
                                            : (pirateIdx ? contacts[*pirateIdx].ship.positionKm() : math::Vec3d{0,0,0});
      const bool hasTarget = hostile || (bool)pirateIdx;
      if (hasTarget) {
        const math::Vec3d to = targetPos - c.ship.positionKm();
        const double dist = to.length();
        if (c.fireCooldown <= 0.0 && dist < 95000.0) {
          const math::Vec3d toN = (dist > 1e-6) ? (to / dist) : math::Vec3d{0,0,1};
          const double aim = math::dot(c.ship.forward().normalized(), toN);
          if (aim > 0.993) {
            c.fireCooldown = 0.28;
            const double dmg = 12.0;

            if (hostile) {
              applyDamage(dmg, playerShield, playerHull);
              const math::Vec3d aKm = c.ship.positionKm();
              const math::Vec3d bKm = ship.positionKm();
              beams.push_back({toRenderU(aKm), toRenderU(bKm), 0.35f, 0.75f, 1.0f, 0.07});
            } else if (pirateIdx) {
              auto& p = contacts[*pirateIdx];
              applyDamage(dmg, p.shield, p.hull);
              const math::Vec3d aKm = c.ship.positionKm();
              const math::Vec3d bKm = p.ship.positionKm();
              beams.push_back({toRenderU(aKm), toRenderU(bKm), 0.35f, 0.75f, 1.0f, 0.07});
              if (p.hull <= 0.0) {
                p.alive = false;
                credits += 180.0;
                toast(toasts, "Security destroyed a pirate (+180).", 2.0);
              }
            }
          }
        }
      }
      continue;
    }
  }

  // Station turrets:
  // - Help against pirates
  // - If you are WANTED with the station's faction, the station will also chip at you near the no-fire zone.
  for (const auto& st : currentSystem->stations) {
    const math::Vec3d stPos = stationPosKm(st, timeDays);
    const double zoneKm = st.radiusKm * 25.0;
    const double distShip = (ship.positionKm() - stPos).length();

    // Only if player is nearby.
    if (distShip > zoneKm) continue;

    // Shoot pirates in zone
    for (auto& c : contacts) {
      if (!c.alive || c.role != ContactRole::Pirate) continue;
      const double d = (c.ship.positionKm() - stPos).length();
      if (d < zoneKm) {
        applyDamage(4.0 * dtSim, c.shield, c.hull);
        if (c.hull <= 0.0) {
          c.alive = false;
          credits += 250.0;
          toast(toasts, "Station defenses destroyed a pirate (+250).", 2.5);
        }
      }
    }

    // Station security vs wanted player
    if (st.factionId != 0 && getBounty(st.factionId) > 0.0) {
      // very light pressure - enough to create urgency without insta-kill
      applyDamage(1.5 * dtSim, playerShield, playerHull);
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

// --- Scanner progress (missions + exploration) ---
if (scanning && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
  bool valid = false;

  // If the player changes target while scanning, cancel.
  if (target.kind != scanLockedTarget.kind || target.index != scanLockedTarget.index) {
    scanning = false;
    scanProgressSec = 0.0;
    scanLockedId = 0;
    scanLabel.clear();
  } else {
    auto completeScan = [&](const std::string& msg, double toastSec) {
      scanning = false;
      scanProgressSec = 0.0;
      scanLockedId = 0;
      scanLabel.clear();
      scanLockedTarget = Target{};
      scanLockedId = 0;
      scanLabel.clear();
      toast(toasts, msg, toastSec);
    };

    // --- CONTACT SCAN (primarily for bounty-scan missions) ---
    if (scanLockedTarget.kind == TargetKind::Contact && scanLockedTarget.index < contacts.size()) {
      auto& c = contacts[scanLockedTarget.index];
      if (c.alive && c.id == scanLockedId) {
        const double dist = (c.ship.positionKm() - ship.positionKm()).length();
        if (dist <= scanRangeKm) {
          valid = true;
          scanProgressSec += dtReal;

          if (scanProgressSec >= scanDurationSec) {
            bool completedAnyMission = false;

            // Find an active bounty-scan mission for this target.
            for (auto& m : missions) {
              if (m.completed || m.failed) continue;
              if (m.type != sim::MissionType::BountyScan) continue;
              if (m.toSystem != currentSystem->stub.id) continue;
              if (m.targetNpcId != c.id) continue;

              m.scanned = true;
              m.completed = true;
              credits += m.reward;
              addRep(m.factionId, +2.0);
              toast(toasts, "Mission complete: bounty scan uploaded! +" + std::to_string((int)m.reward) + " cr", 3.0);
              completedAnyMission = true;
              break;
            }

            if (completedAnyMission) {
              scanning = false;
              scanProgressSec = 0.0;
              scanLockedId = 0;
              scanLabel.clear();
            } else {
              completeScan("Scan complete.", 1.6);
            }
          }
        }
      }
    }

    // --- STATION SCAN ---
    if (!valid && scanLockedTarget.kind == TargetKind::Station && scanLockedTarget.index < currentSystem->stations.size()) {
      const auto& st = currentSystem->stations[scanLockedTarget.index];
      if (st.id == scanLockedId) {
        const double dist = (stationPosKm(st, timeDays) - ship.positionKm()).length();
        if (dist <= scanRangeKm) {
          valid = true;
          scanProgressSec += dtReal;

          if (scanProgressSec >= scanDurationSec) {
            const core::u64 key = scanKeyStation(st.id);
            if (scannedKeys.find(key) == scannedKeys.end()) {
              scannedKeys.insert(key);

              const double value = 180.0 + (double)st.type * 40.0;
              explorationDataCr += value;
              completeScan("Station scan logged (+data " + std::to_string((int)value) + " cr).", 2.5);
            } else {
              completeScan("Station already scanned.", 1.8);
            }
          }
        }
      }
    }

    // --- PLANET SCAN ---
    if (!valid && scanLockedTarget.kind == TargetKind::Planet && scanLockedTarget.index < currentSystem->planets.size()) {
      const auto& p = currentSystem->planets[scanLockedTarget.index];
      const math::Vec3d pPosKm = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
      const double dist = (pPosKm - ship.positionKm()).length();
      if (dist <= scanRangeKm) {
        valid = true;
        scanProgressSec += dtReal;

        if (scanProgressSec >= scanDurationSec) {
          const core::u64 key = scanKeyPlanet(currentSystem->stub.id, scanLockedTarget.index);
          if (scannedKeys.find(key) == scannedKeys.end()) {
            scannedKeys.insert(key);

            double typeMul = 1.0;
            if (p.type == sim::PlanetType::Ocean) typeMul = 1.25;
            if (p.type == sim::PlanetType::Ice) typeMul = 1.15;
            if (p.type == sim::PlanetType::GasGiant) typeMul = 1.45;

            const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
            const double base = 240.0 + radiusKm * 0.03;
            const double value = base * typeMul + (p.orbit.semiMajorAxisAU * 22.0);
            explorationDataCr += value;

            completeScan("Planet scan logged (+data " + std::to_string((int)value) + " cr).", 2.5);
          } else {
            completeScan("Planet already scanned.", 1.8);
          }
        }
      }
    }

    // --- STAR SCAN ---
    if (!valid && scanLockedTarget.kind == TargetKind::Star) {
      valid = true;
      scanProgressSec += dtReal;

      if (scanProgressSec >= scanDurationSec) {
        const core::u64 key = scanKeyStar(currentSystem->stub.id);
        if (scannedKeys.find(key) == scannedKeys.end()) {
          scannedKeys.insert(key);

          const double value = 220.0 + (double)static_cast<int>(currentSystem->star.cls) * 80.0;
          explorationDataCr += value;
          completeScan("Star scan logged (+data " + std::to_string((int)value) + " cr).", 2.2);
        } else {
          completeScan("Star already scanned.", 1.8);
        }
      }
    }

    // If we successfully scanned something, check "system completion" bonus once.
    if (!scanning) {
      const bool starDone = scannedKeys.find(scanKeyStar(currentSystem->stub.id)) != scannedKeys.end();
      bool planetsDone = true;
      for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
        if (scannedKeys.find(scanKeyPlanet(currentSystem->stub.id, i)) == scannedKeys.end()) {
          planetsDone = false;
          break;
        }
      }
      bool stationsDone = true;
      for (const auto& st : currentSystem->stations) {
        if (scannedKeys.find(scanKeyStation(st.id)) == scannedKeys.end()) {
          stationsDone = false;
          break;
        }
      }

      if (starDone && planetsDone && stationsDone) {
        const core::u64 compKey = scanKeySystemComplete(currentSystem->stub.id);
        if (scannedKeys.find(compKey) == scannedKeys.end()) {
          scannedKeys.insert(compKey);
          const double bonus = 600.0 + 90.0 * (double)currentSystem->planets.size() + 120.0 * (double)currentSystem->stations.size();
          explorationDataCr += bonus;
          const core::u32 lf = currentSystem ? currentSystem->stub.factionId : 0;
          if (lf != 0) addRep(lf, +1.0);
          toast(toasts, "System survey complete! Bonus data +" + std::to_string((int)bonus) + " cr.", 3.0);
        }
      }
    }

    if (!valid) {
      // Cancel silently if the target became invalid or we drifted out of range.
      scanning = false;
      scanProgressSec = 0.0;
      scanLockedId = 0;
      scanLabel.clear();
    }
  }
} else if (!scanning) {
  scanProgressSec = 0.0;
}

      // --- Heat model (real-time) ---
      {
        double heatIn = 0.0;
        if (!docked) {
          if (input.boost) heatIn += 18.0;
          if (supercruiseState == SupercruiseState::Active) heatIn += 6.0;
          if (fsdState == FsdState::Charging) heatIn += 12.0;
          if (fsdState == FsdState::Jumping) heatIn += 16.0;
        }

        const double coolRate = docked ? 20.0 : 10.0;
        heat += (heatIn - coolRate) * dtReal;
        heat = std::clamp(heat, 0.0, 120.0);

        // Very simple overheat consequence for now: take slow hull damage.
        if (heat > 100.0) {
          const double over = heat - 100.0;
          playerHull = std::max(0.0, playerHull - over * 0.08 * dtReal);
        }
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
                const econ::CommodityId cid = m.commodity;
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

      // Projectiles update (simple ballistic slugs)
      for (auto& p : projectiles) {
        if (p.ttlSim <= 0.0) continue;

        const math::Vec3d a = p.posKm;
        const math::Vec3d b = p.posKm + p.velKmS * dtSim;

        p.prevKm = a;
        p.posKm = b;
        p.ttlSim -= dtSim;

        // Collisions
        if (p.fromPlayer) {
          for (int i = 0; i < (int)contacts.size(); ++i) {
            const auto& c = contacts[(std::size_t)i];
            if (!c.alive) continue;
            if (p.shooterId != 0 && c.id == p.shooterId) continue;

            const double hitRadiusKm = 900.0 + p.radiusKm;
            if (segmentHitsSphere(a, b, c.ship.positionKm(), hitRadiusKm)) {
              playerDamageContact(i, p.dmg);
              p.ttlSim = 0.0;
              break;
            }
          }
        } else {
          if (!docked) {
            const double playerHitRadiusKm = 900.0 + p.radiusKm;
            if (segmentHitsSphere(a, b, ship.positionKm(), playerHitRadiusKm)) {
              applyDamage(p.dmg, playerShield, playerHull);
              p.ttlSim = 0.0;
            }
          }
        }
      }

      projectiles.erase(
        std::remove_if(projectiles.begin(), projectiles.end(), [](const Projectile& p) { return p.ttlSim <= 0.0; }),
        projectiles.end());


      playerLaserCooldown = std::max(0.0, playerLaserCooldown - dtSim);
      playerCannonCooldown = std::max(0.0, playerCannonCooldown - dtSim);

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
      supercruiseState = SupercruiseState::Idle;
      supercruiseChargeRemainingSec = 0.0;
      supercruiseCooldownRemainingSec = 0.0;
      supercruiseDropRequested = false;
      scanning = false;
      scanProgressSec = 0.0;
      fsdState = FsdState::Idle;
      fsdTargetSystem = 0;
      fsdChargeRemainingSec = 0.0;
      fsdTravelRemainingSec = 0.0;
      fsdTravelTotalSec = 0.0;
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

      // docking corridor guidance for targeted station
      if (target.kind == TargetKind::Station && target.index == i) {
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const math::Quatd stQ = stationOrient(st, stPos, timeDays);
        const math::Vec3d axis = stQ.rotate({0,0,1});
        const math::Vec3d right = stQ.rotate({1,0,0});
        const math::Vec3d up = stQ.rotate({0,1,0});

        const double zEntrance = st.radiusKm * 1.10;
        const double zEnd = zEntrance + st.approachLengthKm;

        const math::Vec3d startKm = stPos + axis * zEntrance;
        const math::Vec3d endKm = stPos + axis * zEnd;

        bool clearanceGranted = false;
        if (auto it = clearances.find(st.id); it != clearances.end()) {
          clearanceGranted = it->second.granted && (timeDays <= it->second.expiresDays);
        }

        const float cr = clearanceGranted ? 0.35f : 0.85f;
        const float cg = clearanceGranted ? 0.85f : 0.55f;
        const float cb = clearanceGranted ? 0.45f : 0.25f;

        auto addLineKm = [&](const math::Vec3d& aKm, const math::Vec3d& bKm, float r, float g, float b) {
          const math::Vec3d aU = toRenderU(aKm);
          const math::Vec3d bU = toRenderU(bKm);
          lines.push_back({(float)aU.x,(float)aU.y,(float)aU.z, r,g,b});
          lines.push_back({(float)bU.x,(float)bU.y,(float)bU.z, r,g,b});
        };

        // Centerline
        addLineKm(startKm, endKm, cr,cg,cb);

        // Corridor wireframe box
        const double hx = st.approachRadiusKm;
        const double hy = st.approachRadiusKm;

        math::Vec3d c0[4] = {
          startKm + right*hx + up*hy,
          startKm + right*hx - up*hy,
          startKm - right*hx - up*hy,
          startKm - right*hx + up*hy,
        };
        math::Vec3d c1[4] = {
          endKm + right*hx + up*hy,
          endKm + right*hx - up*hy,
          endKm - right*hx - up*hy,
          endKm - right*hx + up*hy,
        };

        for (int k = 0; k < 4; ++k) {
          const int k2 = (k + 1) % 4;
          addLineKm(c0[k], c0[k2], cr,cg,cb);
          addLineKm(c1[k], c1[k2], cr,cg,cb);
          addLineKm(c0[k], c1[k], cr,cg,cb);
        }

        // Slot frame at the corridor entrance (yellow)
        const double sw = st.slotWidthKm * 0.5;
        const double sh = st.slotHeightKm * 0.5;
        math::Vec3d s[4] = {
          startKm + right*sw + up*sh,
          startKm + right*sw - up*sh,
          startKm - right*sw - up*sh,
          startKm - right*sw + up*sh,
        };
        for (int k = 0; k < 4; ++k) {
          addLineKm(s[k], s[(k + 1) % 4], 0.95f,0.9f,0.15f);
        }

        // Small staging cross at the corridor end (cyan)
        const double cross = std::max(250.0, st.approachRadiusKm * 0.15);
        addLineKm(endKm - right*cross, endKm + right*cross, 0.25f,0.7f,1.0f);
        addLineKm(endKm - up*cross, endKm + up*cross, 0.25f,0.7f,1.0f);
      }
    }

    // Laser beams
    for (const auto& b : beams) {
      lines.push_back({(float)b.aU.x,(float)b.aU.y,(float)b.aU.z, b.r,b.g,b.b});
      lines.push_back({(float)b.bU.x,(float)b.bU.y,(float)b.bU.z, b.r,b.g,b.b});
    }

    // Projectile tracers (draw a fixed-length tail so they're visible at astronomical scales)
    for (const auto& p : projectiles) {
      math::Vec3d tailKm = p.prevKm;
      if (p.velKmS.lengthSq() > 1e-12) {
        const double tracerLenKm = 15000.0;
        tailKm = p.posKm - p.velKmS.normalized() * tracerLenKm;
      }
      const math::Vec3d aU = toRenderU(tailKm);
      const math::Vec3d bU = toRenderU(p.posKm);
      lines.push_back({(float)aU.x,(float)aU.y,(float)aU.z, p.r,p.g,p.b});
      lines.push_back({(float)bU.x,(float)bU.y,(float)bU.z, p.r,p.g,p.b});
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
      float r = 0.65f, g = 0.75f, b = 0.85f;
      if (c.role == ContactRole::Pirate) { r = 1.0f; g = 0.25f; b = 0.25f; }
      if (c.role == ContactRole::Police) { r = 0.35f; g = 0.75f; b = 1.0f; }
      if (c.role == ContactRole::Trader) { r = 0.45f; g = 0.95f; b = 0.45f; }
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
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();

    // DockSpace is only available in Dear ImGui's docking branch; keep UI
    // functional without it (windows will simply be floating).

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
      } else if (target.kind == TargetKind::Star) {
        tgtKm = math::Vec3d{0,0,0};
        tgtLabel = std::string("Star ") + starClassName(currentSystem->star.cls);
      } else if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
        const auto& c = contacts[target.index];
        if (c.alive) {
          tgtKm = c.ship.positionKm();
          tgtLabel = c.name + std::string(" [") + contactRoleName(c.role) + "]";
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

      // Docking corridor HUD for targeted station
      if (currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
        const auto& st = currentSystem->stations[target.index];

        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const math::Quatd stQ = stationOrient(st, stPos, timeDays);
        const math::Vec3d stV = stationVelKmS(st, timeDays);

        const math::Vec3d relLocal = stQ.conjugate().rotate(ship.positionKm() - stPos);
        const math::Vec3d vRelLocal = stQ.conjugate().rotate(ship.velocityKmS() - stV);

        const double zEntrance = st.radiusKm * 1.10;
        const double zEnd = zEntrance + st.approachLengthKm;

        const double hx = st.approachRadiusKm;
        const double hy = st.approachRadiusKm;

        const bool inZ = (relLocal.z >= zEntrance && relLocal.z <= zEnd);
        const bool inXY = (std::abs(relLocal.x) <= hx && std::abs(relLocal.y) <= hy);
        const bool inCorridor = inZ && inXY;

        bool clearanceGranted = false;
        if (auto it = clearances.find(st.id); it != clearances.end()) {
          clearanceGranted = it->second.granted && (timeDays <= it->second.expiresDays);
        }

        const math::Vec3d axisOut = stQ.rotate({0,0,1});
        const math::Vec3d slotFwd = (-axisOut).normalized();
        const math::Vec3d slotUp = stQ.rotate({0,1,0}).normalized();

        const double fwdAlign = math::dot(ship.forward().normalized(), slotFwd);
        const double upAlign = math::dot(ship.up().normalized(), slotUp);

        const double relSpeed = (ship.velocityKmS() - stV).length();

        ImVec2 ds = io.DisplaySize;
        float y0 = ds.y - 78.0f;

        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "Corridor: %s | offX %.0f km offY %.0f km | z %.0f km | vrel %.2f km/s",
                      inCorridor ? "IN" : "OUT",
                      relLocal.x, relLocal.y, relLocal.z, relSpeed);

        const ImU32 col = clearanceGranted ? IM_COL32(160,255,180,220) : IM_COL32(255,220,140,220);
        draw->AddText({18.0f, y0}, col, buf);

        std::snprintf(buf, sizeof(buf),
                      "Align: fwd %.2f up %.2f | speed limit %.2f km/s | %s",
                      fwdAlign, upAlign, st.maxApproachSpeedKmS,
                      clearanceGranted ? "CLEARANCE GRANTED" : "NO CLEARANCE (press L to request)");
        draw->AddText({18.0f, y0 + 18.0f}, col, buf);
      }

      // Toasts
      float y = 18.0f;
      for (const auto& t : toasts) {
        draw->AddText({18.0f, y}, IM_COL32(240,240,240,220), t.text.c_str());
        y += 18.0f;
      }
    }

if (showShip) {
  ImGui::Begin("Ship / Status");

  ImGui::Text("System: %s", currentSystem->stub.name.c_str());

  const core::u32 localFaction = currentSystem ? currentSystem->stub.factionId : 0;
  if (localFaction != 0) {
    ImGui::Text("Local faction: %s | Rep %.1f | Bounty %.0f",
                factionName(localFaction).c_str(),
                getRep(localFaction),
                getBounty(localFaction));
  } else {
    ImGui::Text("Local faction: (none)");
  }

  ImGui::Separator();

  ImGui::Text("Credits: %.0f | Exploration data: %.0f cr", credits, explorationDataCr);
  ImGui::Text("Fuel: %.1f | Heat: %.0f", fuel, heat);
  ImGui::Text("Cargo: %.0f / %.0f kg", cargoMassKg(cargo), cargoCapacityKg);
  ImGui::Text("Shield: %.0f Hull: %.0f", playerShield, playerHull);

  ImGui::Text("Ship pos: (%.0f,%.0f,%.0f) km", ship.positionKm().x, ship.positionKm().y, ship.positionKm().z);
  ImGui::Text("Vel: %.3f km/s", ship.velocityKmS().length());

  ImGui::Separator();
  ImGui::Text("Controls");
  ImGui::Checkbox("Mouse steer (M)", &mouseSteer);
  ImGui::Checkbox("Invert mouse Y", &mouseInvertY);
  ImGui::SliderFloat("Mouse sensitivity", &mouseSensitivity, 0.0006f, 0.0080f, "%.4f");
  ImGui::TextDisabled("Mouse steer captures the cursor (relative mouse mode).");

  if (scanning) {
    ImGui::Separator();
    ImGui::TextColored(ImVec4(0.9f, 0.95f, 1.0f, 1.0f), "%s", scanLabel.empty() ? "Scanning..." : scanLabel.c_str());
    const float frac = (float)std::clamp(scanProgressSec / std::max(0.01, scanDurationSec), 0.0, 1.0);
    ImGui::ProgressBar(frac, ImVec2(-1, 0));
    ImGui::TextDisabled("K cancels scan.");
  }

  // FSD status (Charging / Jumping / Cooling)
  {
    const bool cooling = (fsdState == FsdState::Idle && timeDays < fsdReadyDay);
    if (fsdState != FsdState::Idle || cooling) {
      ImGui::Separator();

      const char* label = "Idle";
      double timer = 0.0;
      double total = 1.0;

      if (fsdState == FsdState::Charging) {
        label = "Charging";
        total = kFsdChargeSec;
        timer = std::clamp(total - fsdChargeRemainingSec, 0.0, total);
      } else if (fsdState == FsdState::Jumping) {
        label = "Jumping";
        total = std::max(0.001, fsdTravelTotalSec);
        timer = std::clamp(total - fsdTravelRemainingSec, 0.0, total);
      } else {
        // Cooldown uses simulated time.
        label = "Cooling";
        total = kFsdCooldownSec;
        const double remain = std::max(0.0, (fsdReadyDay - timeDays) * 86400.0);
        timer = std::clamp(total - remain, 0.0, total);
      }

      ImGui::Text("FSD: %s", label);
      ImGui::ProgressBar((float)std::clamp(timer / total, 0.0, 1.0), ImVec2(-1, 0));
    }
  }

  if (supercruiseState != SupercruiseState::Idle) {
    ImGui::Separator();

    if (supercruiseState == SupercruiseState::Charging) {
      ImGui::TextColored(ImVec4(0.6f, 0.9f, 1.0f, 1.0f), "Supercruise: CHARGING");
      const float t = (float)std::clamp(1.0 - (supercruiseChargeRemainingSec / std::max(0.001, kSupercruiseChargeSec)), 0.0, 1.0);
      ImGui::ProgressBar(t, ImVec2(-1, 0));
    } else if (supercruiseState == SupercruiseState::Active) {
      ImGui::TextColored(ImVec4(0.6f, 0.9f, 1.0f, 1.0f), "Supercruise: ACTIVE");
      if (supercruiseDistKm > 0.0) ImGui::Text("Range: %.0f km", supercruiseDistKm);
      if (supercruiseTtaSec > 0.0 && supercruiseTtaSec < 1e8) ImGui::Text("ETA: %.1f s", supercruiseTtaSec);

      if (supercruiseSafeDropReady) {
        ImGui::TextColored(ImVec4(0.2f, 1.0f, 0.2f, 1.0f), "SAFE DROP");
      } else {
        ImGui::TextColored(ImVec4(1.0f, 0.9f, 0.2f, 1.0f), "Keep ~7s ETA for safe drop");
      }

      if (!supercruiseAssist) {
        ImGui::Text("Manual drop: press H");
      } else {
        ImGui::Text("Nav assist: auto-drop when safe");
      }
    } else if (supercruiseState == SupercruiseState::Cooldown) {
      ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Supercruise: COOLDOWN");
      const float t = (float)std::clamp(1.0 - (supercruiseCooldownRemainingSec / std::max(0.001, supercruiseCooldownTotalSec)), 0.0, 1.0);
      ImGui::ProgressBar(t, ImVec2(-1, 0));
    }
  }

  // Navigation / route status
  if (navRoute.size() >= 2) {
    ImGui::Separator();
    const int totalJumps = (int)navRoute.size() - 1;
    const int remaining = std::max(0, totalJumps - (int)navRouteHop);

    double remDist = 0.0;
    double remFuel = 0.0;
    for (std::size_t i = navRouteHop; i + 1 < navRoute.size(); ++i) {
      const auto& a = universe.getSystem(navRoute[i]).stub;
      const auto& b = universe.getSystem(navRoute[i + 1]).stub;
      const double d = (b.posLy - a.posLy).length();
      remDist += d;
      remFuel += fsdFuelCostFor(d);
    }

    ImGui::Text("Route: %d jumps (%d remaining) | est fuel %.1f", totalJumps, remaining, remFuel);
    ImGui::TextDisabled("Auto-run: %s", navAutoRun ? "ON" : "OFF");

    if (remFuel > fuel + 1e-6) {
      ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.45f, 1.0f), "Fuel shortfall: +%.1f", remFuel - fuel);
    }
  }

  // Target summary
  ImGui::Separator();
  if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
    const auto& st = currentSystem->stations[target.index];
    ImGui::Text("Target: %s (%.0f km)", st.name.c_str(), (stationPosKm(st, timeDays) - ship.positionKm()).length());
  } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
    const auto& p = currentSystem->planets[target.index];
    const math::Vec3d pPos = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
    ImGui::Text("Target: %s (%.0f km)", p.name.c_str(), (pPos - ship.positionKm()).length());
  } else if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
    const auto& c = contacts[target.index];
    ImGui::Text("Target: %s [%s] (%.0f km)", c.name.c_str(), contactRoleName(c.role), (c.ship.positionKm() - ship.positionKm()).length());
  } else if (target.kind == TargetKind::Star) {
    ImGui::Text("Target: Star (%s)", starClassName(currentSystem->star.cls));
  } else {
    ImGui::Text("Target: (none)  [T station / B planet / N contact / U star]");
  }

  ImGui::Separator();
  ImGui::TextDisabled("Controls:");
  ImGui::BulletText("WASD / Space/Ctrl: translate | Arrows: pitch/yaw | Q/E: roll");
  ImGui::BulletText("Shift: boost | X: brake | Left click: laser | Right click: cannon");
  ImGui::BulletText("H: supercruise (charge/engage). While active: H requests drop (~7s safe)");
  ImGui::BulletText("J engage FSD jump (uses plotted route if present)");
  ImGui::BulletText("P: autopilot to station (staging + corridor guidance)");
  ImGui::BulletText("T/B/N/U cycle targets, Y clear target");
  ImGui::BulletText("K scan target (missions + exploration scans)");
  ImGui::BulletText("L request docking clearance");
  ImGui::BulletText("TAB Galaxy map, F1 Ship, F2 Market, F3 Contacts, F4 Missions, F6 Scanner");
  ImGui::BulletText("F5 quicksave, F9 quickload");

  ImGui::Separator();
  float scMax = (float)supercruiseMaxSpeedKmS;
  if (ImGui::SliderFloat("Supercruise max speed (km/s)", &scMax, 4000.0f, 60000.0f, "%.0f")) {
    supercruiseMaxSpeedKmS = (double)scMax;
  }
  ImGui::Checkbox("Supercruise assist", &supercruiseAssist);

  ImGui::Separator();
  ImGui::TextDisabled("Flight Assist:");
  ImGui::Checkbox("Local reference frame", &localFrameEnabled);
  float lfTau = (float)localFrameBlendTauSec;
  if (ImGui::SliderFloat("Local frame blend (s)", &lfTau, 0.05f, 6.0f, "%.2f")) {
    localFrameBlendTauSec = (double)lfTau;
  }
  ImGui::Checkbox("Blue-zone turn assist", &blueZoneTurnAssist);

  ImGui::End();
}

if (showScanner) {
  ImGui::Begin("System Scanner");

  ImGui::Text("System: %s", currentSystem->stub.name.c_str());
  ImGui::Text("Star: %s | Mass %.2f | Radius %.2f | Temp %.0fK",
              starClassName(currentSystem->star.cls),
              currentSystem->star.massSol,
              currentSystem->star.radiusSol,
              currentSystem->star.temperatureK);

  const bool starScanned = (scannedKeys.find(scanKeyStar(currentSystem->stub.id)) != scannedKeys.end());

  // Quick UI helper: start a scan for a given target.
  auto uiStartScan = [&](TargetKind kind, std::size_t idx) {
    if (scanning) {
      toast(toasts, "Already scanning. (K cancels)", 2.0);
      return;
    }
    if (docked) {
      toast(toasts, "Undock to scan.", 1.8);
      return;
    }
    if (supercruiseState != SupercruiseState::Idle || fsdState != FsdState::Idle) {
      toast(toasts, "Scanning unavailable in supercruise / hyperspace.", 2.0);
      return;
    }

    target.kind = kind;
    target.index = idx;

    scanLockedTarget = target;
    scanLockedId = 0;
    scanLabel.clear();
    scanProgressSec = 0.0;

    bool ok = false;

    if (kind == TargetKind::Star) {
      scanLockedId = currentSystem->stub.id;
      scanDurationSec = 2.5;
      scanRangeKm = 1.0e18;
      scanLabel = std::string("Star scan: ") + starClassName(currentSystem->star.cls);
      ok = true;
    } else if (kind == TargetKind::Planet && idx < currentSystem->planets.size()) {
      const auto& p = currentSystem->planets[idx];
      scanLockedId = (core::u64)idx;
      scanDurationSec = 5.0;
      const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
      scanRangeKm = std::max(200000.0, rKm * 45.0);
      scanLabel = "Planet scan: " + p.name;
      ok = true;
    } else if (kind == TargetKind::Station && idx < currentSystem->stations.size()) {
      const auto& st = currentSystem->stations[idx];
      scanLockedId = st.id;
      scanDurationSec = 3.0;
      scanRangeKm = std::max(25000.0, st.commsRangeKm * 0.9);
      scanLabel = "Station scan: " + st.name;
      ok = true;
    }

    if (ok) {
      scanning = true;
      toast(toasts, scanLabel + " (hold steady)...", 2.0);
    } else {
      toast(toasts, "Invalid scan target.", 2.0);
    }
  };

  ImGui::Separator();

  // Survey status
  int scannedPlanets = 0;
  for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
    if (scannedKeys.find(scanKeyPlanet(currentSystem->stub.id, i)) != scannedKeys.end()) ++scannedPlanets;
  }
  int scannedStations = 0;
  for (const auto& st : currentSystem->stations) {
    if (scannedKeys.find(scanKeyStation(st.id)) != scannedKeys.end()) ++scannedStations;
  }

  ImGui::Text("Survey: Star %s | Planets %d/%d | Stations %d/%d",
              starScanned ? "OK" : "UNSCANNED",
              scannedPlanets, (int)currentSystem->planets.size(),
              scannedStations, (int)currentSystem->stations.size());

  if (!starScanned) {
    if (ImGui::Button("Target Star")) {
      target.kind = TargetKind::Star;
      target.index = 0;
    }
    ImGui::SameLine();
    if (ImGui::Button("Scan Star")) {
      uiStartScan(TargetKind::Star, 0);
    }
  }

  ImGui::Separator();
  ImGui::Text("Planets");

  for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
    const auto& p = currentSystem->planets[i];
    const bool scanned = (scannedKeys.find(scanKeyPlanet(currentSystem->stub.id, i)) != scannedKeys.end());

    ImGui::PushID((int)i);

    ImGui::Text("%s %s", p.name.c_str(), scanned ? "" : "(unscanned)");
    ImGui::SameLine();
    if (ImGui::SmallButton("Target##p")) {
      target.kind = TargetKind::Planet;
      target.index = i;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Scan##p")) {
      uiStartScan(TargetKind::Planet, i);
    }

    if (scanned) {
      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const double periodDays = p.orbit.periodDays;
      ImGui::TextDisabled("Type: %s | Radius: %.0f km | a: %.2f AU | Period: %.1f d",
                          planetTypeName(p.type), radiusKm, p.orbit.semiMajorAxisAU, periodDays);
    } else {
      ImGui::TextDisabled("Type: ??? | Radius: ??? | Orbit: ???");
    }

    ImGui::Separator();
    ImGui::PopID();
  }

  ImGui::Text("Stations");

  for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
    const auto& st = currentSystem->stations[i];
    const bool scanned = (scannedKeys.find(scanKeyStation(st.id)) != scannedKeys.end());

    ImGui::PushID((int)(1000 + i));

    ImGui::Text("%s %s", st.name.c_str(), scanned ? "" : "(unscanned)");
    ImGui::SameLine();
    if (ImGui::SmallButton("Target##s")) {
      target.kind = TargetKind::Station;
      target.index = i;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Scan##s")) {
      uiStartScan(TargetKind::Station, i);
    }

    if (scanned) {
      ImGui::TextDisabled("Type: %s | Faction: %s | Fee: %.1f%%",
                          stationTypeName(st.type),
                          factionName(st.factionId).c_str(),
                          st.feeRate * 100.0);
    } else {
      ImGui::TextDisabled("Type: ??? | Faction: ??? | Services: ???");
    }

    ImGui::Separator();
    ImGui::PopID();
  }

  ImGui::Separator();
  ImGui::Text("Exploration data bank: %.0f cr", explorationDataCr);
  if (docked && selectedStationIndex >= 0 && selectedStationIndex < (int)currentSystem->stations.size()) {
    if (explorationDataCr > 0.0) {
      if (ImGui::Button("Sell exploration data here")) {
        credits += explorationDataCr;
        toast(toasts, "Sold exploration data +" + std::to_string((int)explorationDataCr) + " cr.", 2.5);
        explorationDataCr = 0.0;
      }
    } else {
      ImGui::TextDisabled("No data to sell.");
    }
  } else {
    ImGui::TextDisabled("Dock at a station to sell exploration data.");
  }

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

// Exploration / legal services
if (canTrade) {
  ImGui::Separator();
  ImGui::Text("Exploration Data");
  ImGui::Text("Bank: %.0f cr", explorationDataCr);
  if (explorationDataCr > 0.0) {
    if (ImGui::Button("Sell all exploration data")) {
      credits += explorationDataCr;
      toast(toasts, "Sold exploration data +" + std::to_string((int)explorationDataCr) + " cr.", 2.5);
      explorationDataCr = 0.0;
      addRep(station.factionId, +0.5);
    }
  } else {
    ImGui::TextDisabled("No scan data to sell.");
  }

  ImGui::Separator();
  ImGui::Text("Legal");
  const double bounty = getBounty(station.factionId);
  ImGui::Text("Outstanding bounty (this faction): %.0f cr", bounty);
  if (bounty > 0.0) {
    if (ImGui::Button("Pay bounty")) {
      if (credits >= bounty) {
        credits -= bounty;
        clearBounty(station.factionId);
        toast(toasts, "Bounty cleared.", 2.0);
      } else {
        toast(toasts, "Not enough credits to pay bounty.", 2.0);
      }
    }
  } else {
    ImGui::TextDisabled("No bounty.");
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
            const econ::CommodityId cid = m.commodity;
            out = "Delivery: Deliver " + std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name)
                + " to " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::MultiDelivery: {
            const econ::CommodityId cid = m.commodity;
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
                    const auto cid = (econ::CommodityId)(mrng.nextU32() % (core::u32)econ::kCommodityCount);
                    m.commodity = cid;
                    m.units = 25 + (mrng.nextU32() % 120);
                    m.reward = 250.0 + distLy * 120.0 + (double)m.units * 6.0;
                    m.cargoProvided = (mrng.nextUnit() < 0.25);
                  } else if (r < 0.85 && dests.size() >= 2) {
                    m.type = sim::MissionType::MultiDelivery;
                    const auto cid = (econ::CommodityId)(mrng.nextU32() % (core::u32)econ::kCommodityCount);
                    m.commodity = cid;
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
                  const econ::CommodityId cid = offer.commodity;
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
                    const econ::CommodityId cid = m.commodity;
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
  int aliveTraders = 0;
  int alivePolice = 0;
  for (const auto& c : contacts) {
    if (!c.alive) continue;
    if (c.role == ContactRole::Pirate) ++alivePirates;
    if (c.role == ContactRole::Trader) ++aliveTraders;
    if (c.role == ContactRole::Police) ++alivePolice;
  }

  const core::u32 localFaction = currentSystem ? currentSystem->stub.factionId : 0;
  if (localFaction != 0) {
    ImGui::Text("Local faction: %s | Rep %.1f | Bounty %.0f",
                factionName(localFaction).c_str(),
                getRep(localFaction),
                getBounty(localFaction));
  } else {
    ImGui::Text("Local faction: (none)");
  }

  ImGui::Text("Pirates: %d   Traders: %d   Police: %d", alivePirates, aliveTraders, alivePolice);

  if (ImGui::Button("Panic: clear contacts")) {
    for (auto& c : contacts) c.alive = false;
    toast(toasts, "Contacts cleared.", 2.0);
  }

  ImGui::Separator();

  for (std::size_t i = 0; i < contacts.size(); ++i) {
    const auto& c = contacts[i];
    if (!c.alive) continue;

    const double distKm = (c.ship.positionKm() - ship.positionKm()).length();

    ImGui::PushID((int)i);

    std::string tag = std::string("[") + contactRoleName(c.role) + "]";
    if (c.missionTarget) tag += " [BOUNTY]";
    if (c.role == ContactRole::Police) {
      const bool hostile = c.hostileToPlayer || (c.factionId != 0 && getBounty(c.factionId) > 0.0);
      if (hostile) tag += " [HOSTILE]";
    }

    ImGui::Text("%s %s", c.name.c_str(), tag.c_str());
    if (c.factionId != 0) {
      ImGui::SameLine();
      ImGui::TextDisabled("(%s)", factionName(c.factionId).c_str());
    }

    ImGui::SameLine();
    if (ImGui::SmallButton("Target")) {
      target.kind = TargetKind::Contact;
      target.index = i;
    }

    ImGui::TextDisabled("Dist %.0f km | Hull %.0f | Shield %.0f", distKm, c.hull, c.shield);

    if (c.role == ContactRole::Trader) {
      ImGui::TextDisabled("Cargo value ~%.0f cr (piracy is illegal).", c.cargoValueCr);
    }

    ImGui::Separator();
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
            const int remaining = std::max(0, totalJumps - (int)navRouteHop);
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
