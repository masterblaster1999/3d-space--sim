# Convenience Makefile for CI templates that expect ./configure && make
# Delegates to CMake out-of-source builds.

BUILD_DIR ?= build
BUILD_TYPE ?= Release

# Extra args passed through to the CMake configure step.
# Example: make CMAKE_ARGS="-DSTELLAR_ENABLE_RENDER=OFF -DSTELLAR_ENABLE_IMGUI=OFF"
CMAKE_ARGS ?=

.PHONY: all configure clean test check distcheck

all: configure
	cmake --build $(BUILD_DIR) --config $(BUILD_TYPE)

configure:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) $(CMAKE_ARGS)

test: all
	ctest --test-dir $(BUILD_DIR) --output-on-failure

check: test

# "distcheck" is an Autotools convention. Provide a reasonable stand-in so
# CI templates that run `make distcheck` keep working.
distcheck: clean check

clean:
	rm -rf $(BUILD_DIR)
