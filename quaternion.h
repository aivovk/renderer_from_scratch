#ifndef QUATERNION_H
#define QUATERNION_H

#include <cmath>

#include "geometry.h"

template <typename T = float>
struct Quaternion{
  T w;
  union{
    struct{T i, j, k;};
    Vector<T,3> v;
  };
  Quaternion(T w = 0, T i = 0, T j = 0, T k = 0):w(w), i(i), j(j), k(k) {}
  Quaternion(T w, Vector<T,3> v):w(w), v(v) {}

  Quaternion conjugate() const{
    return Quaternion{w,-v};
  }

  // this is the normalization SQUARED
  T norm() const{
    return w*w + v*v;
  }

  void normalize(){
    auto n = sqrt(norm());
    w /= n;
    v /= n;
  }

  Quaternion inverse() const{
    return conjugate() / norm();
  }

  T selection() const{
    return w;
  }

  T dot(const Quaternion& rhs) const{
    return w*w.rhs + v*v.rhs;
  }

  Quaternion product(const Quaternion& rhs) const{
    return Quaternion{w*rhs.w-v*rhs.v,
      w*rhs.v + rhs.w * v + cross(v, rhs.v)};
  }
  Quaternion operator*(const Quaternion& rhs) const{
    return product(rhs);
  }
  
  Quaternion& operator+=(const Quaternion& rhs){
    w+=rhs.w;
    i+=rhs.i;
    j+=rhs.j;
    k+=rhs.k;
  }

  bool operator==(const Quaternion& rhs) const{
    return w==rhs.w && v==rhs.v;
  }
};

template <class T, std::size_t n>
void Vector<T,n>::rotate(const Vec3f& axis, float angle){
  Quaternion<> q{cos(0.5 * angle), axis * static_cast<float>(sin(0.5 * angle))};
  auto rotated = q * Quaternion<>{0, Vec3f{x, y, z}} * q.conjugate();
  x = rotated.i;
  y = rotated.j;
  z = rotated.k;
}

#endif
