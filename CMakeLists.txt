cmake_minimum_required(VERSION 3.20)

project(StellarForge
  VERSION 0.1.0
  DESCRIPTION "Procedural space simulation starter framework"
  LANGUAGES CXX
)

option(STELLAR_BUILD_SANDBOX "Build the sandbox executable" ON)
option(STELLAR_BUILD_TESTS "Build tests" ON)

# Library
add_library(stellar)
add_library(stellar::stellar ALIAS stellar)

target_compile_features(stellar PUBLIC cxx_std_20)

target_include_directories(stellar
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Warnings / compile options
if (MSVC)
  target_compile_options(stellar PRIVATE /W4 /permissive-)
else()
  target_compile_options(stellar PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
endif()

# Sources
target_sources(stellar PRIVATE
  src/core/Hash.cpp
  src/core/Log.cpp

  src/math/Vec3.cpp

  src/proc/GalaxyGenerator.cpp
  src/proc/NameGenerator.cpp
  src/proc/Noise.cpp
  src/proc/SystemGenerator.cpp

  src/sim/Orbit.cpp
  src/sim/Universe.cpp
)

# Nice folder structure for IDEs
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
  src/core/Hash.cpp
  src/core/Log.cpp
  src/math/Vec3.cpp
  src/proc/GalaxyGenerator.cpp
  src/proc/NameGenerator.cpp
  src/proc/Noise.cpp
  src/proc/SystemGenerator.cpp
  src/sim/Orbit.cpp
  src/sim/Universe.cpp
)

if (STELLAR_BUILD_SANDBOX)
  add_subdirectory(apps/stellar_sandbox)
endif()

if (STELLAR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
