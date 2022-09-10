/* 
 * File:   Pixel.cpp
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#include "Pixel.h"

Pixel::Pixel(double p_red, double p_green, double p_blue, bool p_valid) : m_red {p_red}, m_green {p_green}, m_blue {p_blue}, m_valid {p_valid} {}

const double Pixel::get_red() const {
  return m_red;
}

void Pixel::set_red(double p_d) {
  m_red = p_d;
}

const double Pixel::get_green() const {
  return m_green;
}

void Pixel::set_green(double p_d) {
  m_green = p_d;
}

const double Pixel::get_blue() const {
  return m_blue;
}

void Pixel::set_blue(double p_d) {
  m_blue = p_d;
}
  
const bool Pixel::is_valid() const {
  return m_valid;
}
  
void Pixel::set_valid(bool p_b) {
  m_valid = p_b;
}

Pixel Pixel::operator+(const Pixel& p_pixel) const {
  return Pixel(m_red + p_pixel.m_red, m_green + p_pixel.m_green, m_blue + p_pixel.m_blue);
}

Pixel Pixel::operator-(const Pixel& p_pixel) const {
  return Pixel(m_red - p_pixel.m_red, m_green - p_pixel.m_green, m_blue - p_pixel.m_blue);
}

Pixel operator*(const Pixel& p_pixel, const double p_d) {
  return Pixel(p_d*p_pixel.m_red, p_d*p_pixel.m_green, p_d*p_pixel.m_blue);
}

Pixel operator*(const double p_d, const Pixel& p_pixel) {
  return Pixel(p_d*p_pixel.m_red, p_d*p_pixel.m_green, p_d*p_pixel.m_blue);
}

Pixel Pixel::operator/(const double p_d) const {
  return Pixel(m_red/p_d, m_green/p_d, m_blue/p_d);
}

Pixel& Pixel::operator+=(const Pixel& p_pixel) {
  return (*this = *this + p_pixel);
}