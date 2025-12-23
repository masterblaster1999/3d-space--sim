#pragma once

#include "stellar/core/Types.h"

#include <string_view>

namespace stellar::core {

// 64-bit FNV-1a hash (stable, fast, good for IDs / seeds).
u64 fnv1a64(std::string_view text);

// Combine two 64-bit hashes into one.
u64 hashCombine(u64 a, u64 b);

} // namespace stellar::core
