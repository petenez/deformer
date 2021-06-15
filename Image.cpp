/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Image.cpp
 * Author: pete
 * 
 * Created on June 13, 2021, 10:29 PM
 */

#include <cmath>

#include "Image.h"

// constructs an image from an input file
Image::Image(std::string input) {
  std::ifstream reader(input, std::ios::binary);
  reader.read((char*)&Nx_, sizeof(int));
  reader.read((char*)&Ny_, sizeof(int));
  pixels_.reserve(Nx_*Ny_);
  double r;
  double g;
  double b;
  for(int j = 0; j < Ny_; ++j) {
    for(int i = 0; i < Nx_; ++i) {
      reader.read((char*)&r, sizeof(double));
      reader.read((char*)&g, sizeof(double));
      reader.read((char*)&b, sizeof(double));
      pixels_.push_back(std::make_shared<Pixel>(i, j, r, g, b));
    }
  }
  reader.close();
  set_neighbors();
}

// constructs a blank image
Image::Image(const int Nx, const int Ny) : Nx_(Nx), Ny_(Ny) {
  pixels_.reserve(Nx_*Ny_);
  for(int j = 0; j < Ny; ++j) {
    for(int i = 0; i < Nx; ++i) {
      pixels_.push_back(std::make_shared<Pixel>(i, j));
    }
  }
}

// sets the neighboring pixels (for relaxing displacements in the future)
void Image::set_neighbors() {
  int n = 0;
  for(int j = 0; j < Ny_; ++j) {
    for(int i = 0; i < Nx_; ++i) {
      for(int j_ = std::max(0, j - 1); j_ < std::min(Ny_, j + 2); ++j_) {
        for(int i_ = std::max(0, i - 1); i_ < std::min(Nx_, i + 2); ++i_) {
          int n_ = Nx_*j_ + i_;
          if(abs(i_ - i) + abs(j_ - j) == 1) {
            pixels_[n]->add_adjacent(pixels_[n_]);
          }
          if(abs(i_ - i) + abs(j_ - j) == 2) {
            pixels_[n]->add_diagonal(pixels_[n_]);
          }
        }
      }
      n++;
    }
  }
}

// returns the width of an image
const int Image::get_Nx() const {
  return Nx_;
}

// returns the height of an image
const int Image::get_Ny() const {
  return Ny_;
}

// applies arbitrary deformation to an image
void Image::deform() {
  for(int j = 0; j < Ny_; ++j) {
    for(int i = 0; i < Nx_; ++i) {
      double u = 2.5*sin(2.0*M_PI*i/50.0);
      double v = 1.5*cos(2.0*M_PI*i/30.0) + 1.0*sin(2.0*M_PI*j/40.0);
      int n = Nx_*j + i;
      pixels_[n]->set_i(i + u);
      pixels_[n]->set_j(j + v);
    }
  }
}

// TODO: replace with interpolation from deformed image (no Gaussian blur)
void Image::read(Image image) {
  double sigma = 1.0;
  double a = -0.5/sigma/sigma;
  for(int j = 0; j < image.Ny_; ++j) {
    for(int i = 0; i < image.Nx_; ++i) {
      int n = image.Nx_*j + i;
      double i1 = image.pixels_[n]->get_i();
      double j1 = image.pixels_[n]->get_j();
      int i0 = std::max(0, (int)floor(i1 - 5.0*sigma));
      int i2 = std::min(Nx_ - 1, (int)ceil(i1 + 5.0*sigma));
      int j0 = std::max(0, (int)floor(j1 - 5.0*sigma));
      int j2 = std::min(Ny_ - 1, (int)ceil(j1 + 5.0*sigma));
      for(int j_ = j0; j_ <= j2; ++j_) {
        double dj = j_ - j1;
        for(int i_ = i0; i_ <= i2; ++i_) {
          double di = i_ - i1;
          int n_ = Nx_*j_ + i_;
          double w = exp(a*(di*di + dj*dj));
          pixels_[n_]->set_r(pixels_[n_]->get_r() + w*image.pixels_[n]->get_r());
          pixels_[n_]->set_g(pixels_[n_]->get_g() + w*image.pixels_[n]->get_g());
          pixels_[n_]->set_b(pixels_[n_]->get_b() + w*image.pixels_[n]->get_b());
        }
      }
    }
  }
}

// normalizes an image
void Image::normalize() {
  double max = 0.0;
  for(int n = 0; n < pixels_.size(); ++n) {
    max = std::max(max, pixels_[n]->get_r());
    max = std::max(max, pixels_[n]->get_g());
    max = std::max(max, pixels_[n]->get_b());
  }
  for(int n = 0; n < pixels_.size(); ++n) {
    pixels_[n]->set_r(std::min(1.0, pixels_[n]->get_r()/max));
    pixels_[n]->set_g(std::min(1.0, pixels_[n]->get_g()/max));
    pixels_[n]->set_b(std::min(1.0, pixels_[n]->get_b()/max));
  }
}

// writes an image into an output file
void Image::write(std::string output) {
  normalize();
  std::ofstream writer(output, std::ios::binary);
  writer.write((char*)&Nx_, sizeof(int));
  writer.write((char*)&Ny_, sizeof(int));
  for(int n = 0; n < pixels_.size(); ++n) {
    double r = pixels_[n]->get_r();
    double g = pixels_[n]->get_g();
    double b = pixels_[n]->get_b();
    writer.write((char*)&r, sizeof(double));
    writer.write((char*)&g, sizeof(double));
    writer.write((char*)&b, sizeof(double));
  }
  writer.close();
}