add_executable(stellar_tests
  main.cpp
  test_orbit.cpp
  test_proc.cpp
)

target_link_libraries(stellar_tests PRIVATE stellar::stellar)

add_test(NAME stellar_tests COMMAND stellar_tests)

set_target_properties(stellar_tests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
