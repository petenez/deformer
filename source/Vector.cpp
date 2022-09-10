/* 
 * File:   Vector.cpp
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#include <cmath>
#include <assert.h>

#include "Vector.h"

// Constructor
Vector::Vector(const double p_x, const double p_y) : m_x {p_x}, m_y {p_y} {
  assert(("Error: Vector::Vector: nan detected!", m_x == m_x && m_y == m_y));
}

// Returns the x component of the vector.
double Vector::x() const { return m_x; }

// Returns the y component of the vector.
double Vector::y() const { return m_y; }

// Returns the magnitude of the vector.
double Vector::magnitude() const {
  return std::sqrt(m_x*m_x + m_y*m_y);
}

// operator overloads for vector arithmetic

Vector Vector::operator-() const {
  return -1.0*(*this);
}

Vector Vector::operator+(const Vector& p_v) const {
  return Vector(m_x + p_v.m_x, m_y + p_v.m_y);
}

Vector Vector::operator-(const Vector& p_v) const {
  return Vector(m_x - p_v.m_x, m_y - p_v.m_y);
}

Vector operator*(const Vector& p_v, const double p_d) {
  return Vector(p_d*p_v.m_x, p_d*p_v.m_y);
}

Vector operator*(const double p_d, const Vector& p_v) {
  return Vector(p_d*p_v.m_x, p_d*p_v.m_y);
}

Vector Vector::operator/(const double p_d) const {
  return Vector(m_x/p_d, m_y/p_d);
}

Vector& Vector::operator+=(const Vector& p_v) {
  return (*this = *this + p_v);
}

Vector& Vector::operator-=(const Vector& p_v) {
  return (*this = *this - p_v);
}

Vector& Vector::operator*=(const double p_d) {
  return (*this = p_d*(*this));
}

Vector& Vector::operator/=(const double p_d) {
  return (*this = *this/p_d);
}

// Returns a normalized version of the vector.
Vector Vector::normalize() const {
  double mag {magnitude()};
  return *this/mag;
}

// Returns a rotated version of the vector.
Vector Vector::rotate(double p_d) const {
  double cosd {std::cos(p_d)};
  double sind {std::sin(p_d)};
  return Vector(m_x*cosd - m_y*sind, m_x*sind + m_y*cosd);
}

// Returns the angle between the vectors.
double Vector::angle(const Vector& p_u, const Vector& p_v) {
  return std::acos(dot_product(p_u, p_v)/(p_u.magnitude()*p_v.magnitude()));
}

// Returns the dot product of the vectors
double Vector::dot_product(const Vector& p_u, const Vector& p_v) {
  return p_u.m_x*p_v.m_x + p_u.m_y*p_v.m_y;
}

// Returns the cross product of vectors (0, 0, p_d) and p_u.
Vector Vector::cross_product(const double& p_d, const Vector& p_u) {
  return Vector(-p_d*p_u.m_y, p_d*p_u.m_x);
}

// Returns the cross product of the vectors.
double Vector::cross_product(const Vector& p_u, const Vector& p_v) {
  return p_u.m_x*p_v.m_y - p_u.m_y*p_v.m_x;
}

// Returns the vector projection of p_u onto p_v.
Vector Vector::projection(const Vector& p_u, const Vector& p_v) {
  return dot_product(p_u, p_v)/dot_product(p_v, p_v)*p_v;
}