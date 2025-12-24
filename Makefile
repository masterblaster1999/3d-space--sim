# Convenience Makefile for CI templates that expect ./configure && make
# Delegates to CMake out-of-source builds.

BUILD_DIR ?= build
BUILD_TYPE ?= Release

.PHONY: all configure clean test check

all: configure
	cmake --build $(BUILD_DIR) --config $(BUILD_TYPE)

configure:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE)

test: all
	ctest --test-dir $(BUILD_DIR) --output-on-failure

check: test

clean:
	rm -rf $(BUILD_DIR)
