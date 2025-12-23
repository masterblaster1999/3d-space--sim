#include "stellar/core/Log.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>

namespace stellar::core {
namespace {
  std::mutex g_logMutex;
  LogLevel g_level = LogLevel::Info;

  std::tm localTime(std::time_t t) {
    std::tm out{};
  #if defined(_WIN32)
    localtime_s(&out, &t);
  #else
    localtime_r(&t, &out);
  #endif
    return out;
  }

  std::string nowTimeString() {
    using namespace std::chrono;
    const auto now = system_clock::now();
    const auto tt = system_clock::to_time_t(now);
    const auto tm = localTime(tt);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%H:%M:%S");
    return oss.str();
  }
} // namespace

void setLogLevel(LogLevel level) { g_level = level; }
LogLevel getLogLevel() { return g_level; }

std::string_view toString(LogLevel level) {
  switch (level) {
    case LogLevel::Trace: return "TRACE";
    case LogLevel::Debug: return "DEBUG";
    case LogLevel::Info:  return "INFO";
    case LogLevel::Warn:  return "WARN";
    case LogLevel::Error: return "ERROR";
    case LogLevel::Off:   return "OFF";
  }
  return "UNKNOWN";
}

void log(LogLevel level, std::string_view message) {
  if (level < g_level || g_level == LogLevel::Off) {
    return;
  }

  std::lock_guard<std::mutex> lock(g_logMutex);

  std::ostream& os = (level >= LogLevel::Warn) ? std::cerr : std::cout;
  os << "[" << nowTimeString() << "]"
     << "[" << toString(level) << "] "
     << message
     << "\n";
}

} // namespace stellar::core
