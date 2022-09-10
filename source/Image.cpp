/* 
 * File:   Image.cpp
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#include <cmath>
#include <iostream>
#include <fstream>
// remember to build release version with -DNDEBUG --> disables assertions
#include <assert.h>

#include "Image.h"

namespace {
  // for looping over a pixel's neighbors
  // mth neighbor of nth pixel (i, j) has indices i + di[m] and j + dj[m]
  const static int di[] { 1, 0, -1, 0 };
  const static int dj[] { 0, 1, 0, -1 };
}

// *****

// private:

// Constructor
// Note: Constructs an empty and undeformed image.
Image::Image(const int p_Nx, const int p_Ny) : m_Nx {p_Nx}, m_Ny {p_Ny} {
  if(m_Nx <= 0 || m_Ny <= 0) {
    throw std::invalid_argument("Error: Image::Image: Both image dimensions must be greater than zero!");
  }
  m_pixels = std::vector<Pixel>(m_Nx*m_Ny);
  m_positions = std::vector<Vector>(m_Nx*m_Ny);
}

// Returns products recurring in cubic Lagrangian polynomials.
// Note: Returns \prod_{j = -1..2, j != i} (x - xj)
double Image::product(const int p_i, const double p_x) {
  assert(("Error: Image::product: The index must be -1, 0, 1 or 2!", p_i >= -1 && p_i <= 2));
  double p {1.0};
  for(int j {-1}; j < 3; ++j) {
    if(j != p_i) {
      p *= p_x - j;
    }
  }
  return p;
}

// Returns a cubic Lagrange polynomial at 'p_x'.
// Note: x0 = -1, x1 = 0, x2 = 1, x3 = 2
template <typename T>
T Image::Lagrange(const double p_x, const std::vector<T>& p_y) {
  return product(-1, p_x)/product(-1, -1.0)*p_y[0]
       + product(0, p_x)/product(0, 0.0)*p_y[1]
       + product(1, p_x)/product(1, 1.0)*p_y[2]
       + product(2, p_x)/product(2, 2.0)*p_y[3];
}

// Interpolates a value to 'p_p' from 'p_array' at 'p_q'.
// Note: Uses bicubic interpolation.
template <typename T>
bool Image::interpolate(const Vector& p_q, T& p_p, std::vector<T>& p_array, const int p_Nx) {
  assert(("Error: Image::interpolate: Data width must be greater than zero!", p_Nx > 0));
  const int Ny {(int)p_array.size()/p_Nx};
  const int x_ = std::floor(p_q.x());
  const int y_ = std::floor(p_q.y());
  const double x {p_q.x() - x_};
  const double y {p_q.y() - y_};
  const int i0 {x_ - 1};
  const int j0 {y_ - 1};
  const int i3 {i0 + 3};
  const int j3 {j0 + 3};
  // if surrounding 4-by-4 grid for bicubic interpolation within data
  if(i0 >= 0 && i3 < p_Nx && j0 >= 0 && j3 < Ny) {
    std::vector<T> py;
    for(int j {j0}; j <= j3; ++j) {
      std::vector<T> px;
      for(int i {i0}; i <= i3; ++i) {
        px.push_back(p_array[p_Nx*j + i]);
      }
      py.push_back(Lagrange(x, px));
    }
    p_p = Lagrange(y, py);
    return true;
  }
  return false;
}

// Returns the minimum color component in 'p_pixels'.
double Image::get_min(const std::vector<Pixel>& p_pixels) {
  double min {0.0};
  #pragma omp parallel for reduction(min: min)
  for(int n = 0; n < p_pixels.size(); ++n) {
    if(p_pixels[n].is_valid()) {
      min = std::min(min, p_pixels[n].get_red());
      min = std::min(min, p_pixels[n].get_green());
      min = std::min(min, p_pixels[n].get_blue());
    }
  }
  return min;
}

// Returns the maximum color component in 'p_pixels'.
double Image::get_max(const std::vector<Pixel>& p_pixels) {
  double max {0.0};
  #pragma omp parallel for reduction(max: max)
  for(int n = 0; n < p_pixels.size(); ++n) {
    if(p_pixels[n].is_valid()) {
      max = std::max(max, p_pixels[n].get_red());
      max = std::max(max, p_pixels[n].get_green());
      max = std::max(max, p_pixels[n].get_blue());
    }
  }
  return max;
}
// Sets the pixels according to the convolution of the source image and a Gaussian kernel.
// Note: See the documentation for more details.
void Image::convolve(const Image& p_source, const double p_sigma) {
  assert(("Error: Image::convolve: Gaussian kernel spread must be greater than zero!", p_sigma > 0.0));
  // Gaussian kernel: exp(a*r^2)
  const double a {-0.5/p_sigma/p_sigma};
  // loop over source image pixels
  for(int j_ {0}; j_ < p_source.m_Ny; ++j_) {
    for(int i_ {0}; i_ < p_source.m_Nx; ++i_) {
      // linear index
      const int n_ {p_source.m_Nx*j_ + i_};
      if(p_source.m_pixels[n_].is_valid()) {
        const Vector& p {p_source.m_positions[n_]};
        const double i1 {p.x()};
        const double j1 {p.y()};
        // kernel is truncated to approx. 5 sigma
        const int i0 {std::max(0, (int)std::floor(i1 - 5.0*p_sigma))};
        const int i2 {std::min(m_Nx - 1, (int)std::ceil(i1 + 5.0*p_sigma))};
        const int j0 {std::max(0, (int)std::floor(j1 - 5.0*p_sigma))};
        const int j2 {std::min(m_Ny - 1, (int)std::ceil(j1 + 5.0*p_sigma))};
        // loop over target image pixels
        for(int j {j0}; j <= j2; ++j) {
          const double dj2 {(j - j1)*(j - j1)};
          for(int i {i0}; i <= i2; ++i) {
            const double di {i - i1};
            const int n {m_Nx*j + i};
            const double w {std::exp(a*(di*di + dj2))};
            m_pixels[n] += w*p_source.m_pixels[n_];
          }
        }
      }
    }
  }
}

// Normalizes an image.
void Image::normalize(const double p_gamma) {
  if(p_gamma <= 0.0) {
    std::cout << "Warning: Image::normalize: Gamma correction factor must be greater than zero! "
      "No changes were made!" << std::endl;
  }
  const double max {get_max(m_pixels)};
  if(max == 0.0) {
    std::cout << "Warning: Image::normalize: Maximum pixel value must not be zero! "
      "No changes were made!" << std::endl;
    return;
  }
  const double _max {1.0/max};
  #pragma omp parallel for
  for(int n = 0; n < m_pixels.size(); n++) {
    Pixel& pixel {m_pixels[n]};
    if(pixel.is_valid()) {
      // make sure capped to 1.0
      pixel.set_red(std::pow(std::min(1.0, pixel.get_red()*_max), p_gamma));
      pixel.set_green(std::pow(std::min(1.0, pixel.get_green()*_max), p_gamma));
      pixel.set_blue(std::pow(std::min(1.0, pixel.get_blue()*_max), p_gamma));
    }
  }
}

/* 
 * Interpolates a pixel of the undeformed target image from the deformed source image.
 * Note:
 *   See the documentation for more details.
 *   p_q is the position of a pixel in the space spanned by the source image pixel indices.
 *   p_p is the nearest target image pixel to the aforementioned pixel's physical position.
 */
bool Image::interpolate(const Vector& p_p, const Vector& p_q, Image& p_source) {
  int i1 {0};
  int j1 {0};
  double x {p_q.x()};
  double y {p_q.y()};
  // ten iterations seems to be plenty
  for(int k {0}; k < 10; ++k) {
    const int x_ {(int)std::floor(x)};
    const int y_ {(int)std::floor(y)};
    i1 += x_;
    j1 += y_;
    const int i2 {i1 + 1};
    const int j2 {j1 + 1};
    x -= x_;
    y -= y_;
    // if surrounding 4-by-4 grid for bicubic interpolation within image
    if(i1 >= 1 && i2 < p_source.m_Nx - 1 && j1 >= 1 && j2 < p_source.m_Ny - 1) {
      const int n00 = p_source.m_Nx*j1 + i1;
      const int n10 = p_source.m_Nx*j1 + i2;
      const int n01 = p_source.m_Nx*j2 + i1;
      const int n11 = p_source.m_Nx*j2 + i2;
      if(
        !p_source.m_pixels[n00].is_valid()
        || !p_source.m_pixels[n10].is_valid()
        || !p_source.m_pixels[n01].is_valid()
        || !p_source.m_pixels[n11].is_valid()
      ) {
        return false;
      }
      const Vector p00 {p_source.m_positions[n00]};
      const Vector p10 {p_source.m_positions[n10]};
      const Vector p01 {p_source.m_positions[n01]};
      const Vector p11 {p_source.m_positions[n11]};
      const Vector pi {p10 - p00};
      const Vector pj {p01 - p00};
      const Vector pppp {p11 - p01 - p10 + p00};
      const Vector f {p00 + x*pi + y*pj + x*y*pppp - p_p};
      // if converged
      if(f.magnitude() < 0.001) {
        // bicubic interpolation of pixel value
        const int n {m_Nx*(int)std::round(p_p.y()) + (int)std::round(p_p.x())};
        const Vector q {i1 + x, j1 + y};
        interpolate(q, m_pixels[n], p_source.m_pixels, p_source.m_Nx);
        m_positions[n] = q;
        return true;  // done
      }
      // Newton's method
      const Vector dfdx {pi + y*pppp};
      const Vector dfdy {pj + x*pppp};
      const double a {dfdx.x()};
      const double b {dfdy.x()};
      const double c {dfdx.y()};
      const double d {dfdy.y()};
      const double _det {1.0/(a*d - b*c)};
      x -= _det*Vector::dot_product(Vector {d, -b}, f);
      y -= _det*Vector::dot_product(Vector {-c, a}, f);
    }
    else {
      // iteration fell outside the source image
      return false;
    }
  }
  // iteration didn't converge
  return false;
}

// Sets target image pixels starting from the source image.
// Note: See the documentation for more details.
// for each source image pixel
// 1. find its position in target image coordinates
// 2. find the nearest target image pixel
// 3. if the target image pixel is unset:
// 4.   try to find the target image pixel's position in source image coordinates
// 5.   if 4. is successful:
// 6.     set the target image pixel's value
// 7.     save its position
void Image::interpolate_source(Image& p_source, std::vector<bool>& p_vec_set_pixels) {
  // source image dimensions
  const int Nx {p_source.m_Nx};
  const int Ny {p_source.m_Ny};
  // loop over source image pixels
  for(int v {0}; v < Ny; ++v) {
    for(int u {0}; u < Nx; ++u) {
      // 1D index of pixel (u, v)
      const int m {Nx*v + u};
      // position of mth pixel in target image coordinates
      const Vector& position {p_source.m_positions[m]};
      // 2D indices of the nearest pixel within target image
      const int i = std::min(m_Nx - 1, std::max(0, (int)std::round(position.x())));
      const int j = std::min(m_Ny - 1, std::max(0, (int)std::round(position.y())));
      // 1D index of nearest pixel
      const int n {m_Nx*j + i};
      // if position of nth pixel isn't nullptr i.e. it is unset
      if(!p_vec_set_pixels[n]) {
        // try to find the position p of nth pixel in source image coordinates with Newton's method starting at q = (u, v)
        // try to set nth pixel's value using bicubic interpolation at p
        // if successful, save p in map[n]
        const Vector p {(double)i, (double)j};
        const Vector q {(double)u, (double)v};
        if(interpolate(p, q, p_source)) {
          p_vec_set_pixels[n] = true;
        }
      }
    }
  }
}

// Tries to set the remaining target image pixels starting from target image pixels already set.
// Note: See the documentation for more details.
// find the positions of unset pixels starting from those of their set neighbors and interpolate their values
// 1. while the number of unset pixels is decreasing
// 2.   for each unset pixel
// 3.     loop over its neighbors
// 4.       if the neighbor is valid and set:
// 5.         try to interpolate the unset pixel's value
// 6.         if successful, erase from the set of unset pixels
void Image::interpolate_target(Image& p_source, std::vector<bool>& p_vec_set_pixels, std::set<int>& p_set_unset_pixels) {
  int prev_size {(int)p_set_unset_pixels.size() + 1};
  // while the number of unset pixels is decreasing
  while(p_set_unset_pixels.size() < prev_size) {
    prev_size = p_set_unset_pixels.size();
    std::set<int>::iterator iter {p_set_unset_pixels.begin()};
    // loop over unset pixels
    while(iter != p_set_unset_pixels.end()) {
      const int n {*iter};
      const int i {n%m_Nx};
      const int j {n/m_Nx};
      // loop over neighbors of nth pixel
      for(int m {0}; m < 4; ++m) {
        const int i_ {i + di[m]};
        const int j_ {j + dj[m]};
        // if within target image
        if(i_ >= 0 && i_ < m_Nx && j_ >= 0 && j_ < m_Ny) {
          const int n_ {m_Nx*j_ + i_};
          // if the position of neighboring pixel isn't nullptr i.e. it's set
          if(p_vec_set_pixels[n_]) {
            // initial guess from neighbor's position
            const Vector p {(double)i, (double)j};
            const Vector q {m_positions[n_]};
            // if interpolation is successful
            if(interpolate(p, q, p_source)) {
              p_vec_set_pixels[n] = true;
              // done with nth pixel --> erase
              iter = p_set_unset_pixels.erase(iter);
              break;
            }
          }
        }
      }
      iter++;
    }
  }
}

// Returns a vector of unset pixels connected to the initial unset pixel.
// Note: See the documentation for more details.
// 1. while there are unset pixels:
// 2.   add the first unset pixel into a vector (nicer to loop over) and a set (faster to check if contains elements) of connected unset pixels
// 3.   loop over the vector:
// 4.     if the current pixel is at the edge, save in a bool
// 5.     loop over the current pixel's neighbors:
// 6.       if neighbor is valid, unset and not added, add to the vector and the set
void Image::connect_unset(std::vector<int>& p_vec_connected, const std::vector<bool>& p_vec_set_pixels, bool& p_edge) const {
  std::set<int> set_connected;
  set_connected.insert(p_vec_connected[0]);
  // loop over pixels in vec
  for(int l {0}; l < p_vec_connected.size(); ++l) {
    const int n {p_vec_connected[l]};
    const int i {n%m_Nx};
    const int j {n/m_Nx};
    // if lth pixel is at an edge
    if(i == 0 || i == m_Nx - 1 || j == 0 || j == m_Ny - 1) {
      p_edge = true;
    }
    // loop over neighbors of lth pixel
    for(int m {0}; m < 4; ++m) {
      const int i_ {i + di[m]};
      const int j_ {j + dj[m]};
      // if within target image
      if(i_ >= 0 && i_ < m_Nx && j_ >= 0 && j_ < m_Ny) {
        const int n_ {m_Nx*j_ + i_};
        // if neighboring pixel is unset and it hasn't been added to the set yet
        if(!p_vec_set_pixels[n_] && !set_connected.count(n_)) {
          // add to the vector and set
          p_vec_connected.push_back(n_);
          set_connected.insert(n_);
        }
      }
    }
  }
}

// Diffuses the neighboring pixels' color into the connected pixels.
// Note: See the documentation for more details.
void Image::diffuse_connected(const std::vector<int>& p_vec_connected, std::set<int>& p_set_unset_pixels) {
  // arbitrary convergence threshold
  const double threshold {0.01};
  // arbitrary integration step size
  const double step_size {0.05};
  double max {1.0};
  std::map<int, Pixel> map_connected;
  // iterate until threshold is met
  while(max > threshold) {
    // loop over connected unset pixels
    for(const int n : p_vec_connected) {
      if(m_pixels[n].is_valid()) {
        map_connected[n] = -4.0*m_pixels[n];
        const int i {n%m_Nx};
        const int j {n/m_Nx};
        // loop over neighbors of nth pixel
        for(int m {0}; m < 4; ++m) {
          const int i_ {i + di[m]};
          const int j_ {j + dj[m]};
          const int n_ {m_Nx*j_ + i_};
          // diffuse pixel and its neighbors' values
          map_connected[n] += m_pixels[n_];
        }
      }
    }
    max = 0.0;
    // update pixel values
    for(const int n : p_vec_connected) {
      if(m_pixels[n].is_valid()) {
        const Pixel diff {step_size*map_connected[n]};
        m_pixels[n] += diff;
        max = std::max(max, diff.get_red());
        max = std::max(max, diff.get_green());
        max = std::max(max, diff.get_blue());
      }
    }
  }
  // done with these connected unset pixels --> erase
  for(const int n : p_vec_connected) {
    p_set_unset_pixels.erase(n);
  }
}

// Interpolates an undeformed target image's pixels from a deformed source image.
// Note: See the documentation for more details.
void Image::interpolate(Image& p_source) {
  // a vector to indicate set pixels
  std::vector<bool> vec_set_pixels(m_Nx*m_Ny);
  interpolate_source(p_source, vec_set_pixels);
  // a set to store 1D indices of unset pixels
  std::set<int> set_unset_pixels;
  // loop over source image pixels
  for(int n {0}; n < m_Nx*m_Ny; ++n) {
    // if position of nth pixel is nullptr i.e. it is unset
    if(!vec_set_pixels[n]) {
      set_unset_pixels.insert(n);
    }
  }
  interpolate_target(p_source, vec_set_pixels, set_unset_pixels);
  // deal with possible remaining unset pixels:
  // a. omit ones connected to the edges of the target image
  // b. diffuse the values of neighboring set pixels into unset pixels
  // while there are unset pixels
  while(set_unset_pixels.size() > 0) {
    // (1) find regions of connected unset pixels
    // a vector for connected unset pixels
    std::vector<int> vec_connected {*set_unset_pixels.begin()};
    bool edge {false};
    connect_unset(vec_connected, vec_set_pixels, edge);
    // (2a) if one or more of the connected unset pixels at an edge --> done with them --> erase
    if(edge) {
      for(const int n : vec_connected) {
        set_unset_pixels.erase(n);
      }
    }
    // (2b) diffuse connected unset pixels' values from those of their neighbors
    else {
      diffuse_connected(vec_connected, set_unset_pixels);
    }
  }
}

// *****

// public:

// Constructor
Image::Image(const std::string& p_input) {
  std::ifstream reader(p_input, std::ios::binary);
  if(!reader.is_open()) {
    throw std::invalid_argument("Error: Image::Image: Failed to open the input file!");
  }
  reader.read((char*)&m_Nx, sizeof(int));
  reader.read((char*)&m_Ny, sizeof(int));
  if(m_Nx <= 0 || m_Ny <= 0) {
    throw std::invalid_argument("Error: Image::Image: Both image dimensions must be greater than zero!");
  }
  m_pixels.reserve(m_Nx*m_Ny);
  m_positions.reserve(m_Nx*m_Ny);
  double r;
  double g;
  double b;
  for(int j {0}; j < m_Ny; ++j) {
    for(int i {0}; i < m_Nx; ++i) {
      reader.read((char*)&r, sizeof(double));
      reader.read((char*)&g, sizeof(double));
      reader.read((char*)&b, sizeof(double));
      r = std::min(1.0, std::max(0.0, r));
      g = std::min(1.0, std::max(0.0, g));
      b = std::min(1.0, std::max(0.0, b));
      m_pixels.push_back({r, g, b});
      m_positions.push_back({(double)i, (double)j});
    }
  }
  reader.close();
}

// Scales the positions of pixels.
void Image::scale_positions(const double p_factor) {
  if(p_factor <= 0.0) {
    std::cout << "Warning: Image::scale_positions: The scaling factor must be greater than zero! "
      "No changes were made!" << std::endl;
    return;
  }
  #pragma omp parallel for
  for(int n = 0; n < m_positions.size(); ++n) {
    if(m_pixels[n].is_valid()) {
      m_positions[n] *= p_factor;
    }
  }
}

// Sets the positions of pixels based on the leaf nodes' positions in the quadtree.
void Image::set_positions(const std::shared_ptr<Node>& p_root, const Vector p_origin, const Vector p_du) {
  if(!p_root) {
    std::cout << "Warning: Image::set_positions: The root node must not be a nullptr! "
      "No changes were made!" << std::endl;
    return;
  }
  if(p_du.magnitude() == 0.0) {
    std::cout << "Warning: Image::set_positions: The direction vector must not be a null vector! "
      "No changes were made!" << std::endl;
    return;
  }
  const int size {p_root->get_size()};
  const int size2 {size*size};
  // vector of positions from leaf nodes
  std::vector<Vector> positions(size2);
  #pragma omp parallel for
  for(int n = 0; n < size2; ++n) {
    const int i {n%size};
    const int j {n/size};
    positions[n] = p_root->get_position(i, j);
  }
  const Vector dv {Vector::cross_product(1.0, p_du)};
  #pragma omp parallel for
  for(int v = 0; v < m_Ny; ++v) {
    for(int u = 0; u < m_Nx; ++u) {
      const int n {m_Nx*v + u};
      const Vector q {p_origin + u*p_du + v*dv};
      if(!interpolate(q, m_positions[n], positions, size)) {
        m_pixels[n].set_valid(false);
      }
    }
  }
}

// Returns a set of linear pixel indices from thresholding an image.
std::set<int> Image::threshold(const double p_limit, const bool p_invert) const {
  std::set<int> pixels;
  const double min {get_min(m_pixels)};
  const double max {get_max(m_pixels)};
  if(max > min) {
    if(p_invert) {
      // inverted selection
      #pragma omp parallel
      {
        std::set<int> pixels_local;
        #pragma omp for nowait
        for(int n = 0; n < m_Nx*m_Ny; ++n) {
          if(m_pixels[n].is_valid() && m_pixels[n].get_red() < p_limit && m_pixels[n].get_green() < p_limit && m_pixels[n].get_blue() < p_limit) {
            pixels_local.insert(n);
          }
        }
        #pragma omp critical
        {
          pixels.insert(pixels_local.begin(), pixels_local.end());
        }
      }
    }
    else {
      // direct selection
      #pragma omp parallel
      {
        std::set<int> pixels_local;
        #pragma omp for nowait
        for(int n = 0; n < m_Nx*m_Ny; ++n) {
          if(m_pixels[n].is_valid() && (m_pixels[n].get_red() > p_limit || m_pixels[n].get_green() > p_limit || m_pixels[n].get_blue() > p_limit)) {
            pixels_local.insert(n);
          }
        }
        #pragma omp critical
        {
          pixels.insert(pixels_local.begin(), pixels_local.end());
        }
      }
    }
  }
  return pixels;
}

// Convolves an image from a deformed source image and writes it into an output file.
void Image::write_convolved(const std::string& p_output, const int p_Nx, const int p_Ny, const double p_sigma, const double p_gamma) const {
  if(p_Nx <= 0 || p_Ny <= 0) {
    std::cout << "Warning: Image::write_convolved: Image dimensions must be greater than zero! "
      "No output was generated!" << std::endl;
    return;
  }
  if(p_sigma <= 0) {
    std::cout << "Warning: Image::write_convolved: Gaussian kernel spread must be greater than zero! "
      "No output was generated!" << std::endl;
    return;
  }
  Image target {p_Nx, p_Ny};
  target.convolve(*this, p_sigma);
  target.normalize(p_gamma);
  std::ofstream writer {p_output, std::ios::binary};
  if(!writer.is_open()) {
    std::cout << "Warning: Image::write_convolved: Failed to open the output file! "
      "No output was generated!" << std::endl;
    return;
  }
  writer.write((char*)&p_Nx, sizeof(int));
  writer.write((char*)&p_Ny, sizeof(int));
  for(const Pixel& pixel : target.m_pixels) {
    const bool is_valid {pixel.is_valid()};
    const double red {is_valid ? pixel.get_red() : 0.0};
    const double green {is_valid ? pixel.get_green() : 0.0};
    const double blue {is_valid ? pixel.get_blue() : 0.0};
    writer.write((char*)&red, sizeof(double));
    writer.write((char*)&green, sizeof(double));
    writer.write((char*)&blue, sizeof(double));
  }
  writer.close();
}

// Interpolates an image from a deformed source image and writes it into an output file.
void Image::write_interpolated(const std::string& p_output, const int p_Nx, const int p_Ny, const double p_gamma) {
  if(p_Nx <= 0 || p_Ny <= 0) {
    std::cout << "Warning: Image::write_interpolated: Image dimensions must be greater than zero! "
      "No output was generated!" << std::endl;
    return;
  }
  if(p_gamma <= 0.0) {
    std::cout << "Warning: Image::write_interpolated: Gamma correction factor must be greater than zero! "
      "No output was generated!" << std::endl;
    return;
  }
  Image target(p_Nx, p_Ny);
  target.interpolate(*this);
  target.normalize(p_gamma);
  std::ofstream writer {p_output, std::ios::binary};
  if(!writer.is_open()) {
    std::cout << "Warning: Image::write_interpolated: Failed to open the output file! "
      "No output was generated!" << std::endl;
    return;
  }
  writer.write((char*)&p_Nx, sizeof(int));
  writer.write((char*)&p_Ny, sizeof(int));
  for(const Pixel& pixel : target.m_pixels) {
    const bool is_valid {pixel.is_valid()};
    const double red {is_valid ? pixel.get_red() : 0.0};
    const double green {is_valid ? pixel.get_green() : 0.0};
    const double blue {is_valid ? pixel.get_blue() : 0.0};
    writer.write((char*)&red, sizeof(double));
    writer.write((char*)&green, sizeof(double));
    writer.write((char*)&blue, sizeof(double));
  }
  writer.close();
}