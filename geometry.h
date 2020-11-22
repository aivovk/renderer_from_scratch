/*
 * Generic Vector and Matrix data structures and operations
 *
 * Used for Points, Vertices, Transformations
 *
 * Triangle Classes
 *
 * TODO:
 * - one triangle template
 * - specialize the union for 2D vectors: x,y, instead of x,y,z..
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

template <class T, std::size_t n> struct Vector;
template<class T, std::size_t nrows, std::size_t mcols> class Matrix;

#include "geometry.tpp"

typedef Vector<int,2> Point;
typedef Vector<float,4> Vec4f;
typedef Vector<float,3> Vec3f;
typedef Vector<int,3> Vec3i;
typedef Vector<int,2> Vec2i;
typedef Vector<unsigned char, 4> Color;
typedef Matrix<float, 3, 3> Matrix3f;
typedef Matrix<float, 4, 4> Matrix4f;

struct Triangle;
struct PixelTriangle;



/**
 * Template length vector
 */
template <class T, std::size_t n>
struct Vector
{
  union{
    T data[n];
    struct{T x,y,z;};
    struct{T u,v,d;};
  };
  
  //inline Vector<T,n>():data{} {}
  
  template <typename... Ts>
  inline Vector<T,n>(Ts... data_list):data{data_list...} {}

  template <typename U, std::size_t m>
  inline Vector<T,n>(const Vector<U,m>& rhs):data{} {
    for(std::size_t i = 0 ; i < std::min(n, m) ; i++)
      data[i] = rhs.data[i];
  }
  
  inline T& operator[](std::size_t i) {
    return data[i];
  }

  inline const T& operator[](std::size_t i) const {
    return data[i];
  }

  inline Vector<T,n> operator-() const {
    Vector<T,n> a;
    for( std::size_t i = 0 ; i < n ; i++)
      a.data[i] = -data[i];
    return a;
  }
  inline Vector<T,n>& operator+=(const Vector<T,n> & r) {
    for(std::size_t i = 0 ; i < n ; i++)
      data[i] += r.data[i];
    return *this;
  }
  inline Vector<T,n>& operator-=(const Vector<T,n> & r) {
    for(std::size_t i = 0 ; i < n ; i++)
      data[i] -= r.data[i];
    return *this;
  }
  inline T magnitude() const { return sqrt((*this) * (*this)); }
  inline T magnitudeSquared() const {return (*this) * (*this); } // or this * this;
  inline Vector<T,n>& normalize() {
    for(std::size_t i = 0 ; i < n ; i++)
      data[i] /= magnitude();
    return *this;
  }
  
  bool operator<(const Vector<T,n>& other) const{
    if(this->data[1] != other.data[1])
      return this->data[1] < other.data[1];
    else
      return this->data[0] < other.data[0];
  }
  
};

// output
template <class T, std::size_t n>
std::ostream & operator<<(std::ostream & out, const Vector<T, n> & r);

/// individually multiply x, y, and z of two vectors together
/// (not dot product)
template <class T, std::size_t n, std::size_t m>
inline Vector<T,n> mult(const Vector<T,n> & r1, const Vector<T,m> & r2) {
  Vector<T,n> a = r1;
  for( std::size_t i = 0 ; i < std::min(n, m) ; i++)
    a.data[i] = r1.data[i] * r2.data[i];
  return a;
}
template <class T, std::size_t n>
inline Vector<T,n> operator+(const Vector<T,n> & r1, const Vector<T,n> & r2) {
  Vector<T,n> a;
  for( std::size_t i = 0 ; i < n ; i++)
    a.data[i] = r1.data[i] + r2.data[i];
  return a;
}
template <class T, std::size_t n>
inline Vector<T,n> operator-(const Vector<T,n> & r1, const Vector<T,n> & r2) {
  Vector<T,n> a;
  for( std::size_t i = 0 ; i < n ; i++)
    a.data[i] = r1.data[i] - r2.data[i];
  return a;
}

// dot product
template <class T, std::size_t n>
inline T operator*(const Vector<T,n> & r1, const Vector<T,n> & r2) {
  T res = 0;
  for(std::size_t i = 0 ; i < n ; i++)
    res += r1.data[i] * r2.data[i];
  return res;
}

// multiplication and division by a constant
template <class T, std::size_t n>
inline Vector<T,n> operator*(const Vector<T,n> & r, const T c) {
  Vector<T,n> a;
  for( std::size_t i = 0 ; i < n ; i++)
    a.data[i] = r.data[i] * c;
  return a;
}
template <class T, std::size_t n>
inline Vector<T,n> operator*(const T c, const Vector<T,n> & r) {
  return r * c;
}
template <class T, std::size_t n>
inline Vector<T,n> operator/(const Vector<T,n> & r, const T c) {
  Vector<T,n> a;
  for( std::size_t i = 0 ; i < n ; i++)
    a.data[i] = r.data[i] / c;
  return a;
}

//cross product
template <class T, std::size_t n>
inline Vector<T,n> cross(const Vector<T,n> & r1, const Vector<T,n> & r2) {
  Vector<T,n> result;
  for(std::size_t i = 0 ; i < n ; i++) {
    result.data[i] = r1.data[(i+1)%n] * r2.data[(i+2)%n] - r1.data[(i+2)%n] * r2.data[(i+1)%n];
  }
  return result;
}


template<class T, std::size_t nrows, std::size_t mcols>
class Matrix {
public:
  T data[nrows][mcols];
  Matrix():data{} {};
  static Matrix identity(){
    Matrix m;
    for (std::size_t i = 0 ; i < std::min(nrows, mcols) ; i++)
      m(i,i) = 1;
    return m;
  }
  T& operator() (std::size_t row, std::size_t col) { return data[row][col]; }
  const T& operator() (std::size_t row, std::size_t col) const { return data[row][col]; }

  Matrix<T, mcols, nrows> transposed() {
    Matrix<T, mcols, nrows> m;
    for(int row = 0 ; row < nrows ; row++)
      for(int col = 0 ; col < mcols ; col++)
	m.data[col][row] = data[row][col];
    return m;
  }
};

template<class T, std::size_t arows, std::size_t acols_brows, std::size_t bcols>
Matrix<T, arows, bcols> operator*(const Matrix<T, arows, acols_brows>& A,
					      const Matrix<T, acols_brows, bcols>& B) {
  Matrix<T, arows, bcols> m;
  for( std::size_t row = 0 ; row < arows ; row++)
    for( std::size_t col = 0 ; col < bcols ; col++) 
      for( std::size_t i = 0 ; i < acols_brows; i++)
	m.data[row][col] += A.data[row][i] * B.data[i][col];
  return m;
}

template<class T, std::size_t nrows, std::size_t mcols>
Vector<T, nrows> operator*(const Matrix<T, nrows, mcols>& m,
			   const Vector<T, mcols>& v){
  Vector<T, nrows> a;
  for(std::size_t row = 0 ; row < nrows ; row++)
    a.data[row] = std::inner_product(v.data, v.data + mcols,
				     m.data[row],
				     static_cast<T>(0));
  return a;
}

template <class T, std::size_t nrows, std::size_t mcols>
std::ostream & operator<<(std::ostream & out, const Matrix<T, nrows, mcols> & r);

struct Triangle{
  union{
    Vec4f data[3];
    struct{Vec4f A, B, C;} Vertices;
  };

  Triangle():data{} {}
  
  inline Vec4f& operator[](std::size_t i){
    return data[i];
  }
  inline const Vec4f& operator[](std::size_t i) const{
    return data[i];
  }
};

struct PixelTriangle{
  union{
    Vec2i data[3];
    struct{Vec2i A, B, C;} Vertices;
  };

  PixelTriangle():data{} {}
		
  PixelTriangle(const Triangle& t){
    data[0] = t[0];
    data[1] = t[1];
    data[2] = t[2];
  }
  inline Vector<int,2>& operator[](std::size_t i){
    return data[i];
  }
  inline const Vector<int,2>& operator[](std::size_t i) const{
    return data[i];
  }
};

#endif
