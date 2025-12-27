#pragma once

#include <cctype>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace stellar::core {

// Tiny, dependency-free argument parser for CLI tools.
//
// Supports:
//  - Flags:         --flag   -h
//  - KV args:       --key value   --key=value
//  - Positional:    everything else
//
// Notes:
//  - Unknown/duplicate keys are preserved (last() is convenient but values() exists).
//  - This is intentionally small; it's for tools like stellar_sandbox.
class Args {
public:
  Args() = default;
  Args(int argc, char** argv) { parse(argc, argv); }

  // Configure a long option to consume multiple subsequent values.
  // Example: setArity("pos", 3) allows: --pos 0 0 0
  void setArity(std::string_view key, int valueCount) {
    if (valueCount <= 0) return;
    arity_[std::string(key)] = valueCount;
  }

  void parse(int argc, char** argv) {
    program_.clear();
    kv_.clear();
    flags_.clear();
    positional_.clear();

    if (argc > 0 && argv && argv[0]) program_ = argv[0];

    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i] ? std::string(argv[i]) : std::string();
      if (a.empty()) continue;

      // Long options.
      if (startsWith(a, "--")) {
        const auto eq = a.find('=');
        if (eq != std::string::npos) {
          const std::string key = a.substr(2, eq - 2);
          const std::string val = a.substr(eq + 1);
          kv_[key].push_back(val);
          continue;
        }

        const std::string key = a.substr(2);

        const auto ar = arity_.find(key);
        const int need = (ar != arity_.end()) ? ar->second : 1;

        int took = 0;
        while (took < need && i + 1 < argc && argv[i + 1] && !isSwitch(argv[i + 1])) {
          kv_[key].push_back(std::string(argv[++i]));
          ++took;
        }

        // No value(s) -> treat as a flag.
        if (took == 0) flags_.push_back(key);
        continue;
      }

      // Short flags (-h, -v, -abc)
      if (startsWith(a, "-") && a.size() >= 2 && a[1] != '-') {
        for (std::size_t j = 1; j < a.size(); ++j) {
          const char c = a[j];
          if (std::isalnum((unsigned char)c) || c == '_') {
            flags_.push_back(std::string(1, c));
          }
        }
        continue;
      }

      // Positional.
      positional_.push_back(a);
    }
  }

  const std::string& program() const { return program_; }

  bool hasFlag(std::string_view key) const {
    for (const auto& f : flags_) {
      if (f == key) return true;
    }
    return false;
  }

  bool has(std::string_view key) const {
    if (hasFlag(key)) return true;
    return kv_.find(std::string(key)) != kv_.end();
  }

  std::optional<std::string> last(std::string_view key) const {
    const auto it = kv_.find(std::string(key));
    if (it == kv_.end() || it->second.empty()) return std::nullopt;
    return it->second.back();
  }

  std::vector<std::string> values(std::string_view key) const {
    const auto it = kv_.find(std::string(key));
    if (it == kv_.end()) return {};
    return it->second;
  }

  const std::vector<std::string>& flags() const { return flags_; }
  const std::vector<std::string>& positional() const { return positional_; }

  // Typed helpers (return true if provided & parsed).
  bool getU64(std::string_view key, unsigned long long& out) const {
    const auto v = last(key);
    if (!v) return false;
    char* end = nullptr;
    const auto val = std::strtoull(v->c_str(), &end, 10);
    if (end == v->c_str()) return false;
    out = val;
    return true;
  }

  bool getI64(std::string_view key, long long& out) const {
    const auto v = last(key);
    if (!v) return false;
    char* end = nullptr;
    const auto val = std::strtoll(v->c_str(), &end, 10);
    if (end == v->c_str()) return false;
    out = val;
    return true;
  }

  bool getInt(std::string_view key, int& out) const {
    long long v = 0;
    if (!getI64(key, v)) return false;
    out = static_cast<int>(v);
    return true;
  }

  bool getDouble(std::string_view key, double& out) const {
    const auto v = last(key);
    if (!v) return false;
    char* end = nullptr;
    const auto val = std::strtod(v->c_str(), &end);
    if (end == v->c_str()) return false;
    out = val;
    return true;
  }

  bool getString(std::string_view key, std::string& out) const {
    const auto v = last(key);
    if (!v) return false;
    out = *v;
    return true;
  }

private:
  static bool startsWith(const std::string& s, const char* prefix) {
    const std::size_t n = std::char_traits<char>::length(prefix);
    return s.size() >= n && s.compare(0, n, prefix) == 0;
  }

  static bool isSwitch(const char* s) {
    if (!s || !*s) return false;
    return s[0] == '-';
  }

  std::string program_;
  std::unordered_map<std::string, int> arity_;
  std::unordered_map<std::string, std::vector<std::string>> kv_;
  std::vector<std::string> flags_;
  std::vector<std::string> positional_;
};

} // namespace stellar::core
