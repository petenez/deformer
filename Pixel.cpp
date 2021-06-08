/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Pixel.cpp
 * Author: pete
 * 
 * Created on May 29, 2021, 8:32 AM
 */

#include "Pixel.h"

Pixel::Pixel(int i, int j, double r, double g, double b)
: i0_(i), j0_(j), i_(i), j_(j), r_(r), g_(g), b_(b) {}

const double Pixel::r() const {
  return r_;
}

const double Pixel::g() const {
  return g_;
}

const double Pixel::b() const {
  return b_;
}