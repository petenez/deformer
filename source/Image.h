/* 
 * File:   Image.h
 * Author: Petri Hirvonen
 * 8 April 2022
 */

#ifndef IMAGE_H
#define IMAGE_H

#include <vector>

#include "Pixel.h"
#include "Node.h"

// Images whose pixels have positions in addition to the color information.
// Used for rendering the deformed images.
class Image final {
    
private:
    
  // Image dimensions
  int m_Nx;
  int m_Ny;
  // Color of pixels
  std::vector<Pixel> m_pixels;
  // positions of pixels
  std::vector<Vector> m_positions;
  
  /*
   * Constructor
   * Parameters:
   *   p_Nx - image x size
   *   p_Ny - image y size
   */
  Image(int p_Nx, int p_Ny);
  
  /*
   * Returns the product \prod_{j = 0..3, j != i} (x - xj). Used for Lagrange polynomials.
   * Parameters:
   *   p_i - an integer in [0, 3]
   *   p_x - a real number
   */
  static double product(const int p_i, const double p_x);
  
  /*
   * Returns the value of a cubic Lagrange polynomial at p_x.
   * Parameters:
   *   p_x - the point where the polynomial is to be evaluated at
   *   p_y - the values of the polynomial at neighboring grid points
   */
  template <typename T> static T Lagrange(const double p_x, const std::vector<T>& p_y);
  
  /*
   * Interpolates, using bicubic interpolation, a value from data at a given point.
   * Parameters:
   *   p_q - the point to interpolate at
   *   p_p - the value to be interpolated
   *   p_array - 2D data array
   *   p_Nx - linear size of the data array
   */
  template <typename T> static bool interpolate(const Vector& p_q, T& p_p, std::vector<T>& p_array, const int p_Nx);
  
  // Returns the minimum of all pixels' color components.
  // Parameter: p_pixels - the pixels
  static double get_min(const std::vector<Pixel>& p_pixels);
  
  // Returns the maximum of all pixels' color components.
  // Parameter: p_pixels - the pixels
  static double get_max(const std::vector<Pixel>& p_pixels);
  
  /*
   * Sets the pixels according to the convolution of the source image and a Gaussian kernel.
   * Parameters:
   *   p_source - source image
   *   p_sigma - spread of the Gaussian kernel
   */
  void convolve(const Image& p_source, const double p_sigma = 1.0);
  
  // Normalizes the pixels so that the maximum color component of any pixel is one.
  // Parameters: p_gamma - gamma correction factor
  void normalize(const double p_gamma = 1.0);
  
  /*
   * Interpolates a pixel of the undeformed target image from the deformed source image.
   * Note:
   *   Bicubic interpolation is used to interpolate the colors and positions between the target image pixels.
   *   Newton's method is used to find the positions in the source image corresponding to the pixels' positions in the target image.
   *   The positions found in the source image corresponding to the pixels' positions in the target image are saved.
   *   The saved positions are used as a starting point later when the remaining pixels are interpolated.
   * Parameters:
   *   p_p - position of the target image pixel closest to the position of the source image pixel
   *   p_q - position of the source image pixel
   */
  bool interpolate(const Vector& p_p, const Vector& p_q, Image& p_source);
  
  /*
   * Tries to interpolate the target image pixels with positions closest to the source image pixels.
   * Note:
   *   Keeps track of pixels set successfully.
   * Parameters:
   *   p_source - source image
   *   p_vec_set_pixels - Boolean map indicating the set pixels
   */
  void interpolate_source(Image& p_source, std::vector<bool>& p_vec_set_pixels);
  
  /*
   * Tries to interpolate the target image pixels not already set.
   * Note:
   *   Uses the positions saved in function 'interpolate_source' as starting points.
   * Parameters:
   *   p_source - source image
   *   p_vec_set_pixels - Boolean map indicating the set pixels
   *   p_set_unset_pixels - set indicating the unset pixels
   */
  void interpolate_target(Image& p_source, std::vector<bool>& p_vec_set_pixels, std::set<int>& p_set_unset_pixels);
  
  /*
   * Adds the unset pixels that are connected to the first in the set of unset pixels to a vector.
   * Parameters:
   *   p_vec_connected - vector for the connected unset pixels
   *   p_vec_set_pixels - vector of set pixels
   *   p_edge - Boolean indicating if the connected pixels are connected to an edge
   */
  void connect_unset(std::vector<int>& p_vec_connected, const std::vector<bool>& p_vec_set_pixels, bool& p_edge) const;
  
  /*
   * Diffuses the colors of the connected unset pixels from the neighboring set pixels.
   * Parameters:
   *   p_vec_connected - vector of connected unset pixels
   *   p_set_unset_pixels - set indicating the unset pixels
   */
  void diffuse_connected(const std::vector<int>& p_vec_connected, std::set<int>& p_set_unset_pixels);
  
  /*
   * Interpolates the undeformed target image pixels from the deformed source image pixels.
   * Note:
   *   If some target image pixels cannot be interpolated, their colors are diffused from neighboring pixels.
   * Parameters:
   *   p_source - source image
   */
  void interpolate(Image& p_source);
  
  // Copy constructor
  // Note: Copying is not allowed.
  Image(const Image& p_image) = delete;
  
public:
    
  // Constructor
  // Parameter: p_input - input image file name
  Image(const std::string& p_input);
  
  // Scales the pixels' positions.
  // Parameter: p_d - scaling factor
  void scale_positions(const double p_factor);
  
  /*
   * Sets the positions of the pixels based on the nodes' positions.
   * Parameters:
   *   p_root - root node
   *   p_origin - origin
   *   p_du - direction of the grid of nodes
   */
  void set_positions(const std::shared_ptr<Node>& p_root, const Vector p_origin = Vector{}, const Vector p_du = Vector{1.0});
  
  /*
   * Returns a set of linear indices of pixels from thresholding the image.
   * Parameters:
   *   p_d - a real value to threshold with respect to
   *   p_invert - a bool indicating whether to invert the selection
   */
  std::set<int> threshold(const double p_limit = 0.5, const bool p_invert = false) const;
  
  /*
   * Convolves the image and writes it into an output file.
   * Parameters:
   *   p_output - output filename
   *   p_Nx - image width
   *   p_Ny - image height
   *   p_sigma - spread of Gaussian kernel
   *   p_gamma - gamma correction factor
   */
  void write_convolved(const std::string& p_output, const int p_Nx, const int p_Ny, const double p_sigma = 1.0, const double p_gamma = 1.0) const;
  
  /*
   * Interpolates the image and writes it into an output file.
   * Parameters:
   *   p_output - output filename
   *   p_Nx - image width
   *   p_Ny - image height
   *   p_gamma - gamma correction factor
   */
  void write_interpolated(const std::string& p_output, const int p_Nx, const int p_Ny, const double p_gamma = 1.0);
};

#endif /* IMAGE_H */

