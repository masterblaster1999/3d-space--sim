cmake_minimum_required(VERSION 3.20)

project(stellar_forge LANGUAGES CXX)

include(FetchContent)
include(CTest)

# Legacy option preserved for older build scripts / README snippets.
# Prefer STELLAR_ENABLE_RENDER + STELLAR_ENABLE_IMGUI going forward.
option(STELLAR_BUILD_GAME "(DEPRECATED) Build the SDL/OpenGL game prototype" ON)

option(STELLAR_ENABLE_RENDER "Build OpenGL/SDL rendering module and the game app" ON)
option(STELLAR_ENABLE_IMGUI  "Enable Dear ImGui UI (requires STELLAR_ENABLE_RENDER)" ON)

if (NOT STELLAR_BUILD_GAME)
  set(STELLAR_ENABLE_RENDER OFF CACHE BOOL "" FORCE)
  set(STELLAR_ENABLE_IMGUI OFF CACHE BOOL "" FORCE)
endif()

option(STELLAR_FETCH_SDL2 "Fetch SDL2 if not found via find_package (render builds only)" ON)
option(STELLAR_FETCH_IMGUI "Fetch Dear ImGui (render+imgui builds only)" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Dependencies (conditional) ----
set(STELLAR_SDL_TARGET "")
set(IMGUI_TARGET "")

if (STELLAR_ENABLE_RENDER)
  find_package(OpenGL QUIET)
  if (NOT OpenGL_FOUND)
    message(WARNING "OpenGL not found. Disabling STELLAR_ENABLE_RENDER/STELLAR_ENABLE_IMGUI. Configure with -DSTELLAR_ENABLE_RENDER=OFF for headless builds, or install OpenGL development packages.")
    set(STELLAR_ENABLE_RENDER OFF CACHE BOOL "" FORCE)
    set(STELLAR_ENABLE_IMGUI OFF CACHE BOOL "" FORCE)
  else()

  # SDL2
  find_package(SDL2 CONFIG QUIET)
  if (NOT SDL2_FOUND AND STELLAR_FETCH_SDL2)
    message(STATUS "SDL2 not found via find_package. Fetching SDL2...")
    set(SDL_TEST OFF CACHE BOOL "" FORCE)
    set(SDL_TESTS OFF CACHE BOOL "" FORCE)
    set(SDL_SHARED OFF CACHE BOOL "" FORCE)
    set(SDL_STATIC ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
      sdl2
      GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
      GIT_TAG release-2.30.9
    )
    FetchContent_MakeAvailable(sdl2)
  endif()

  # Determine SDL target name
  if (TARGET SDL2::SDL2)
    set(STELLAR_SDL_TARGET SDL2::SDL2)
  elseif (TARGET SDL2-static)
    set(STELLAR_SDL_TARGET SDL2-static)
  elseif (TARGET SDL2)
    set(STELLAR_SDL_TARGET SDL2)
  endif()

  if (STELLAR_SDL_TARGET STREQUAL "")
    message(FATAL_ERROR "SDL2 target not found. Install SDL2 or enable STELLAR_FETCH_SDL2.")
  endif()

  # Dear ImGui (optional)
  if (STELLAR_ENABLE_IMGUI AND STELLAR_FETCH_IMGUI)
    FetchContent_Declare(
      imgui
      GIT_REPOSITORY https://github.com/ocornut/imgui.git
      GIT_TAG v1.91.2
    )
    FetchContent_MakeAvailable(imgui)

    add_library(imgui_lib STATIC
      ${imgui_SOURCE_DIR}/imgui.cpp
      ${imgui_SOURCE_DIR}/imgui_demo.cpp
      ${imgui_SOURCE_DIR}/imgui_draw.cpp
      ${imgui_SOURCE_DIR}/imgui_tables.cpp
      ${imgui_SOURCE_DIR}/imgui_widgets.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_sdl2.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
    )
    target_include_directories(imgui_lib PUBLIC
      ${imgui_SOURCE_DIR}
      ${imgui_SOURCE_DIR}/backends
    )
    target_link_libraries(imgui_lib PUBLIC ${STELLAR_SDL_TARGET} OpenGL::GL)
    target_compile_definitions(imgui_lib PUBLIC IMGUI_DISABLE_OBSOLETE_FUNCTIONS)
    set(IMGUI_TARGET imgui_lib)
  endif()
  endif() # OpenGL_FOUND
endif()

# ---- Library ----
set(STELLAR_CORE_SOURCES
  src/core/Assert.cpp
  src/core/Hash.cpp
  src/core/Log.cpp

  src/proc/Noise.cpp
  src/proc/NameGenerator.cpp
  src/proc/GalaxyGenerator.cpp
  src/proc/SystemGenerator.cpp

  src/sim/Faction.cpp
  src/sim/Orbit.cpp
  src/sim/Ship.cpp
  src/sim/SaveGame.cpp
  src/sim/Traffic.cpp
  src/sim/Universe.cpp
  src/sim/NavRoute.cpp
  src/sim/MissionLogic.cpp

  src/econ/Commodity.cpp
  src/econ/Economy.cpp
  src/econ/Market.cpp
  src/econ/RoutePlanner.cpp
)

set(STELLAR_RENDER_SOURCES
  src/render/Gl.cpp
  src/render/Shader.cpp
  src/render/Texture.cpp
  src/render/Mesh.cpp
  src/render/Camera.cpp
  src/render/MeshRenderer.cpp
  src/render/PointRenderer.cpp
  src/render/LineRenderer.cpp
)

add_library(stellar STATIC
  ${STELLAR_CORE_SOURCES}
  $<$<BOOL:${STELLAR_ENABLE_RENDER}>:${STELLAR_RENDER_SOURCES}>
)

add_library(stellar::stellar ALIAS stellar)

target_include_directories(stellar PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(stellar PUBLIC cxx_std_20)

if (STELLAR_ENABLE_RENDER)
  target_link_libraries(stellar PUBLIC
    OpenGL::GL
    ${STELLAR_SDL_TARGET}
  )
  if (WIN32)
    target_link_libraries(stellar PUBLIC gdi32)
  endif()
endif()

# ---- Apps ----
add_executable(stellar_sandbox apps/stellar_sandbox/main.cpp)
target_link_libraries(stellar_sandbox PRIVATE stellar::stellar)
set_target_properties(stellar_sandbox PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps
)

if (STELLAR_ENABLE_RENDER AND STELLAR_ENABLE_IMGUI)
  add_executable(stellar_game apps/stellar_game/main.cpp)
  target_link_libraries(stellar_game PRIVATE stellar::stellar ${IMGUI_TARGET})
  set_target_properties(stellar_game PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps
  )
endif()

# ---- Tests ----
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
