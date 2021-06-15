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

// constructs a blank pixel
Pixel::Pixel(int i, int j)
: i_(i), j_(j) {}

// constructs a pixel with given red, green and blue values
Pixel::Pixel(int i, int j, double r, double g, double b)
: i_(i), j_(j), r_(r), g_(g), b_(b) {}

// returns the red value
const double Pixel::get_r() const {
  return r_;
}

// sets the red value
void Pixel::set_r(const double d) {
  r_ = d;
}

// same for green
const double Pixel::get_g() const {
  return g_;
}

void Pixel::set_g(const double d) {
  g_ = d;
}

// same for blue
const double Pixel::get_b() const {
  return b_;
}

void Pixel::set_b(const double d) {
  b_ = d;
}

// returns the horizontal position of a pixel in a deformed image
const double Pixel::get_i() const {
  return i_;
}

// sets the horizontal position of a pixel in a deformed image
void Pixel::set_i(const double d) {
  i_ = d;
}

// same for vertical position
const double Pixel::get_j() const {
  return j_;
}

void Pixel::set_j(const double d) {
  j_ = d;
}

// counts the adjacent (nearest) neighbors of a pixel
const int Pixel::count_adjacent() const {
  return adjacent_.size();
}

// returns the ith adjacent neighbor of a pixel
std::shared_ptr<Pixel> Pixel::get_adjacent(const int i) const {
  return adjacent_[i];
}

// adds an adjacent neighbor to a pixel
void Pixel::add_adjacent(std::shared_ptr<Pixel> pixel) {
  adjacent_.push_back(pixel);
}

// same for diagonal (second-nearest) neighbors
const int Pixel::count_diagonal() const {
  return diagonal_.size();
}

std::shared_ptr<Pixel> Pixel::get_diagonal(const int i) const {
  return diagonal_[i];
}

void Pixel::add_diagonal(std::shared_ptr<Pixel> pixel) {
  diagonal_.push_back(pixel);
}