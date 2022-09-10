/* 
 * File:   Vector.h
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#ifndef VECTOR_H
#define VECTOR_H

// 2D Vectors
class Vector final {
    
private:
    
  // vector components
  double m_x;
  double m_y;
    
public:
  /* 
   * Constructor
   * Parameters:
   *   p_x - x component
   *   p_y - y component
   */
  Vector(double p_x = 0.0, double p_y = 0.0);
  
  // Returns the x component of the vector.
  double x() const;
  
  // Returns the y component of the vector.
  double y() const;
  
  // Returns the magnitude of the vector.
  double magnitude() const;
  
  // operator overloads for vector arithmetic
  Vector operator-() const;
  Vector operator+(const Vector& p_v) const;
  Vector operator-(const Vector& p_v) const;
  friend Vector operator*(const Vector& p_v, const double p_d);
  friend Vector operator*(const double p_d, const Vector& p_v);
  Vector operator/(const double) const;
  Vector& operator+=(const Vector& p_v);
  Vector& operator-=(const Vector& p_v);
  Vector& operator*=(const double p_d);
  Vector& operator/=(const double p_d);
  
  // Returns a normalized version of the vector.
  Vector normalize() const;
  
  /* 
   * Returns a rotated version of the vector.
   * Parameters:
   *   p_d - rotation angle
   */
  Vector rotate(double p_d) const;
  
  /*
   * Returns the angle between the vectors.
   * Parameters:
   *   p_u - a vector
   *   p_v - a vector
   * Returns:
   *   The angle between the vectors in radians
   */
  static double angle(const Vector& p_u, const Vector& p_v);
  
  /*
   * Returns the dot product of the vectors.
   * Parameters:
   *   p_u - a vector
   *   p_v - a vector
   * Returns:
   *   The dot product of the vectors
   */
  static double dot_product(const Vector& p_u, const Vector& p_v);
  
  /*
   * Returns the cross product of vectors (0, 0, p_d) and p_u.
   * Parameters:
   *   p_d - z component of vector (0, 0, p_d)
   *   p_u - a vector
   * Returns:
   *   The cross product of vectors (0, 0, p_d) and p_u
   */
  static Vector cross_product(const double& p_d, const Vector& p_u);
  
  /*
   * Returns the cross product (0, 0, z) of the vectors.
   * Parameters:
   *   p_u - a vector
   *   p_v - a vector
   * Returns:
   *   The cross product (0, 0, z) of the vectors
   */
  static double cross_product(const Vector& p_u, const Vector& p_v);
  
  /*
   * Returns the vector projection of p_u onto p_v.
   * Parameters:
   *   p_u - a vector
   *   p_v - a vector
   * Returns:
   *   The vector projection of p_u onto p_v
   */
  static Vector projection(const Vector& p_u, const Vector& p_v);
};

#endif /* VECTOR_H */

