#pragma once

#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::core {

// Minimal JSON writer with proper string escaping + pretty printing.
// Designed for headless tooling output (stellar_sandbox).
class JsonWriter {
public:
  explicit JsonWriter(std::ostream& out, bool pretty = true, int indentSpaces = 2)
      : out_(out), pretty_(pretty), indentSpaces_(indentSpaces) {}

  void beginObject() { beginScope(Scope::Object); }
  void endObject() { endScope(Scope::Object); }

  void beginArray() { beginScope(Scope::Array); }
  void endArray() { endScope(Scope::Array); }

  void key(std::string_view k) {
    // Keys only valid inside object.
    preValue();
    writeString(k);
    out_ << (pretty_ ? ": " : ":");
    afterKey_ = true;
  }

  void value(std::string_view s) {
    if (!afterKey_) preValue();
    writeString(s);
    postValue();
  }
  void value(const char* s) {
    if (!s) {
      nullValue();
      return;
    }
    value(std::string_view(s));
  }
  void value(const std::string& s) { value(std::string_view(s)); }

  void value(double v) {
    if (!afterKey_) preValue();
    // Avoid scientific notation for small tooling outputs.
    out_ << v;
    postValue();
  }
  void value(long long v) {
    if (!afterKey_) preValue();
    out_ << v;
    postValue();
  }
  void value(unsigned long long v) {
    if (!afterKey_) preValue();
    out_ << v;
    postValue();
  }
  void value(int v) { value((long long)v); }
  void value(bool v) {
    if (!afterKey_) preValue();
    out_ << (v ? "true" : "false");
    postValue();
  }
  void nullValue() {
    if (!afterKey_) preValue();
    out_ << "null";
    postValue();
  }

private:
  enum class Scope { Object, Array };
  struct Frame {
    Scope scope{Scope::Object};
    bool first{true};
  };

  void beginScope(Scope s) {
    if (!afterKey_) preValue();
    out_ << (s == Scope::Object ? '{' : '[');
    stack_.push_back(Frame{s, true});
    if (pretty_) {
      out_ << "\n";
    }
    afterKey_ = false;
  }

  void endScope(Scope s) {
    if (stack_.empty() || stack_.back().scope != s) return;

    const bool wasEmpty = stack_.back().first;
    stack_.pop_back();

    if (pretty_) {
      if (!wasEmpty) out_ << "\n";
      indent();
    }
    out_ << (s == Scope::Object ? '}' : ']');
    postValue();
  }

  void preValue() {
    if (!stack_.empty()) {
      auto& f = stack_.back();
      // In objects, value must come after key.
      if (f.scope == Scope::Object && !afterKey_) {
        // If caller forgot to provide a key, we still separate entries correctly.
      }

      if (!f.first) {
        out_ << ',';
        if (pretty_) out_ << "\n";
      } else {
        // First element in a new scope: newline already emitted on beginScope().
      }

      if (pretty_) indent();
      f.first = false;
    }
  }

  void postValue() {
    afterKey_ = false;
  }

  void indent() {
    if (!pretty_) return;
    const int depth = (int)stack_.size();
    for (int i = 0; i < depth * indentSpaces_; ++i) out_ << ' ';
  }

  static void writeEscaped(std::ostream& o, std::string_view s) {
    for (char c : s) {
      switch (c) {
        case '"': o << "\\\""; break;
        case '\\': o << "\\\\"; break;
        case '\b': o << "\\b"; break;
        case '\f': o << "\\f"; break;
        case '\n': o << "\\n"; break;
        case '\r': o << "\\r"; break;
        case '\t': o << "\\t"; break;
        default:
          // Control chars -> \u00XX
          if ((unsigned char)c < 0x20) {
            static const char* hex = "0123456789abcdef";
            const unsigned char uc = (unsigned char)c;
            o << "\\u00" << hex[(uc >> 4) & 0xF] << hex[uc & 0xF];
          } else {
            o << c;
          }
          break;
      }
    }
  }

  void writeString(std::string_view s) {
    out_ << '"';
    writeEscaped(out_, s);
    out_ << '"';
  }

  std::ostream& out_;
  bool pretty_{true};
  int indentSpaces_{2};
  bool afterKey_{false};
  std::vector<Frame> stack_;
};

} // namespace stellar::core
