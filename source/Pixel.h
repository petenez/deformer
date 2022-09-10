/* 
 * File:   Pixel.h
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#ifndef PIXEL_H
#define PIXEL_H

// Pixels with three color components.
class Pixel final {
    
private:
    
  // pixel color components
  double m_red;
  double m_green;
  double m_blue;
  // pixel validity
  bool m_valid;
  
public:
    
  // Constructor
  Pixel(double p_red = 0.0, double p_green = 0.0, double p_blue = 0.0, bool m_valid = true);
  
  // Returns the red component.
  const double get_red() const;
  
  // Sets the red component.
  void set_red(double p_d);
  
  // Returns the green component.
  const double get_green() const;
  
  // Sets the green component.
  void set_green(double p_d);
  
  // Returns the blue component.
  const double get_blue() const;
  
  // Sets the blue component.
  void set_blue(double p_d);
  
  // Returns whether the pixel is valid.
  const bool is_valid() const;
  
  // Sets the pixel's validity.
  void set_valid(bool p_b);
  
  // Addition overload.
  Pixel operator+(const Pixel& p_pixel) const;
  
  // Subtraction overload.
  Pixel operator-(const Pixel& p_pixel) const;
  
  // Multiplication overload.
  friend Pixel operator*(const Pixel& p_pixel, const double p_d);
  
  // Multiplication overload.
  friend Pixel operator*(const double p_d, const Pixel& p_pixel);
  
  // Division overload.
  Pixel operator/(const double) const;
  
  // Addition assignment overload.
  Pixel& operator+=(const Pixel& p_pixel);
};

#endif /* PIXEL_H */

