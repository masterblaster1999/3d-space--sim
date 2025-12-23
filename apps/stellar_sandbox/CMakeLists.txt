add_executable(stellar_sandbox
  main.cpp
)

target_link_libraries(stellar_sandbox PRIVATE stellar::stellar)

# Put the binary in a predictable folder
set_target_properties(stellar_sandbox PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps/stellar_sandbox
)
