#pragma once

#include <cmath>
#include <ostream>

namespace stellar::math {

struct Vec3d {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  constexpr Vec3d() = default;
  constexpr Vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  constexpr Vec3d operator+(const Vec3d& rhs) const { return {x + rhs.x, y + rhs.y, z + rhs.z}; }
  constexpr Vec3d operator-(const Vec3d& rhs) const { return {x - rhs.x, y - rhs.y, z - rhs.z}; }
  constexpr Vec3d operator*(double s) const { return {x * s, y * s, z * s}; }
  constexpr Vec3d operator/(double s) const { return {x / s, y / s, z / s}; }

  constexpr Vec3d& operator+=(const Vec3d& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
  constexpr Vec3d& operator-=(const Vec3d& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
  constexpr Vec3d& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
  constexpr Vec3d& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

  double lengthSquared() const { return x*x + y*y + z*z; }
  double length() const { return std::sqrt(lengthSquared()); }

  Vec3d normalized() const {
    const double len = length();
    if (len <= 0.0) return {0.0, 0.0, 0.0};
    return *this / len;
  }
};

inline constexpr Vec3d operator*(double s, const Vec3d& v) { return v * s; }

inline double dot(const Vec3d& a, const Vec3d& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Vec3d cross(const Vec3d& a, const Vec3d& b) {
  return {
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
}

inline std::ostream& operator<<(std::ostream& os, const Vec3d& v) {
  os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  return os;
}

} // namespace stellar::math
