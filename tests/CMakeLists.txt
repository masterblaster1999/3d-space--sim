add_executable(stellar_tests
  main.cpp
  test_orbit.cpp
  test_streaming.cpp
  test_proc.cpp
  test_economy.cpp
  test_route_planner.cpp
  test_savegame.cpp
  test_missions.cpp
  test_nav.cpp
)

target_compile_features(stellar_tests PRIVATE cxx_std_20)
target_link_libraries(stellar_tests PRIVATE stellar::stellar)

add_test(NAME stellar_tests COMMAND stellar_tests)

set_target_properties(stellar_tests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
