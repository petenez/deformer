/* 
 * File:   main.cpp
 * Author: Petri Hirvonen
 * 10 September 2022
 */

#include <cmath>
#include <random>
#include <cstdlib>
#include <iostream>
#include <string>
#include <chrono>
#include <set>

#include "Node.h"
#include "Image.h"
#include "Vector.h"

namespace {
  // demo cases
  enum class Demo {
    RESAMPLE,
    UPSAMPLE,
    ZOOM,
    FAN,
    SKEW,
    STRETCH,
    TWIST,
    PINCHES,
    QUENCH,
    SPREAD,
    RESIZE,
    MOSAIC,
    WORLD,
    BOTCH,
  };
}

// input images (I'm providing some links to elsewhere due to copyright reasons)
// baboon: https://cms.uni-konstanz.de/fileadmin/archive/informatik-saupe/fileadmin/informatik/ag-saupe/Webpages/lehre/dip_w0910/pictures/baboon.png
// world map: https://southafrica-info.com/wp-content/uploads/2017/10/world_map_of_countries_by_surface_area.jpg (I cut out the legend at the bottom)
// female statue: https://commons.wikimedia.org/wiki/File:Lilith_Periodo_de_Isin_Larsa_y_Babilonia.JPG
// male statue: https://commons.wikimedia.org/wiki/File:Marble_statue_of_a_kouros_(youth)_MET_DT263.jpg

int main(int argc, char** argv) {
  
  // case to be demonstrated
  Demo demo {Demo::RESAMPLE};
  
  switch(demo) {
    case Demo::RESAMPLE:
    {
      // The input image is simply resampled and gamma correction is applied.
      // The end result is a somewhat downscaled and tilted, darker version
      // of the original.
      
      // input and output strings
      // use "./pngreader input.png output.rgb" to generate rgb files
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      // use "./pngwriter input.rgb output.png to generate png files
      const std::string output {"output.rgb"};
      // path for input and output
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/resample/"};
      
      // node grid linear size
      const int size_nodes {64};
      // output image linear size
      const int size_output {512};
      // scaling factor
      const double scale {8.0};
      // gamma factor
      const double gamma {2.0};
      // image origin
      const Vector origin {5.0, 10.0};
      // pixel offset
      const Vector offset {0.1, -0.01};
      
      // quadtree root
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      // the input image is loaded
      Image image {path + input};
      // pixel positions are interpolated from the grid of nodes
      // image origin is set to (5, 10) in grid coordinates
      // pixel offset in x direction is set to (0.1, -0.01) in grid coordinates
      image.set_positions(root, origin, offset);
      // pixel positions are scaled
      image.scale_positions(8.0);
      // the output image is written
      image.write_interpolated(path + output, size_output, size_output, gamma);
      
      break;
    }
    case Demo::UPSAMPLE:
    {
      // A small input image is upsampled both with convolution and interpolation.
      // The convolved image is composed of blurry blobs.
      // The interpolated image is continuous but also blurry due to a lot of upsampling.
      
      // input and output strings
      const std::string input {"input-mario.rgb"};
      // you can use for example the provided 16x16 Super Mario image (drew myself after
      // some examples online) here
      const std::string output_int {"output-interpolated.rgb"};
      const std::string output_con {"output-convolved.rgb"};
      // path for input and output
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/upsample/"};
      
      const int size_output {512};
      const double scale {32.0};
      // Gaussian spread ("blur radius") for convolution
      const double sigma {5.0};
      
      // the input image is loaded
      Image image {path + input};
      // pixel positions are scaled
      image.scale_positions(scale);
      // output image is interpolated (target image pixels are interpolated from the source image)
      // note that bicubic interpolation crops image edges a bit (use padded images to avoid this)
      image.write_interpolated(path + output_int, size_output, size_output);
      // output image is convolved (source image pixels are convolved with a Gaussian kernel onto the target image)
      image.write_convolved(path + output_con, size_output, size_output, sigma);
      
      break;
    }
    case Demo::ZOOM:
    {
      // Two grids of nodes are created: one with a scale (1.0 by default) and one without.
      // Two corresponding, downsampled images are created: both are stretched identically
      // at two points, but the one with a scale resist stretching whereas the one without
      // a scale relaxes by growing larger.
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      const std::string output_fix {"output-fixed-scale.rgb"};
      const std::string output_no {"output-no-scale.rgb"};
      // path for input and output
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/zoom/"};
      
      const int size_nodes {64};
      const int size_output {512};
      const int step_count {100000};
      const double scale {4.0};
      // step size for relaxation (gradient descent)
      const double step_size {0.05};
      const Vector origin {};
      const Vector offset {0.125};
      
      // a set of all nodes
      std::set<int> set_all;
      for(int n {0}; n < size_nodes*size_nodes; ++n) {
        set_all.insert(n);
      }
      std::shared_ptr<Node> root_fix {Node::create_nodes(size_nodes)};
      std::shared_ptr<Node> root_no {Node::create_nodes(size_nodes)};
      // Translations are applied to the nodes in the set to center the images.
      root_fix->apply_translation(set_all, Vector {size_nodes/2, size_nodes/2});
      root_no->apply_translation(set_all, Vector {size_nodes/2, size_nodes/2});
      // nodes are set not to have a scale
      root_no->set_scale(set_all, 0.0);
      
      // nodes are offset to apply stretching
      int i {size_nodes/4};
      int j {size_nodes/2};
      int di {-size_nodes/4};
      root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
      root_no->set_fixed(i, j);
      i = size_nodes*3/4;
      j = size_nodes/2;
      di = size_nodes/4;
      root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
      root_no->set_fixed(i, j);
      i = size_nodes/4;
      j = size_nodes/2;
      di = -size_nodes/4;
      root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
      root_fix->set_fixed(i, j);
      i = size_nodes*3/4;
      j = size_nodes/2;
      di = size_nodes/4;
      root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
      root_fix->set_fixed(i, j);
      
      // the nodes are relaxed
      root_no->relax(step_count, step_size);
      root_fix->relax(step_count, step_size);
      
      // images are created and written
      {
        Image image {path + input};
        image.set_positions(root_fix, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output_fix, size_output, size_output);
      }
      {
        Image image {path + input};
        image.set_positions(root_no, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output_no, size_output, size_output);
      }
      
      break;
    }
    case Demo::FAN:
    {
      // Two grids of nodes are created: one with a scale (1.0 by default) and one without.
      // Two corresponding, downsampled images are created: both are pinched at their bottom
      // end and stretched at their top end. The one without a scale is better able to
      // compress at the bottom and to expand at the top, but this is still limited due to
      // the angular forces between the nodes.
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      const std::string output_fix {"output-fixed-scale.rgb"};
      const std::string output_no {"output-no-scale.rgb"};
      // path for input and output
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/fan/"};
      
      const int size_nodes {64};
      const int size_output {512};
      const int step_count {100000};
      const double scale {4.0};
      const double step_size {0.05};
      const Vector origin {};
      const Vector offset {0.125};
      
      std::set<int> set_all;
      for(int n {0}; n < size_nodes*size_nodes; ++n) {
        set_all.insert(n);
      }
      std::shared_ptr<Node> root_fix {Node::create_nodes(size_nodes)};
      std::shared_ptr<Node> root_no {Node::create_nodes(size_nodes)};
      root_fix->apply_translation(set_all, Vector {size_nodes/2, size_nodes/2});
      root_no->apply_translation(set_all, Vector {size_nodes/2, size_nodes/2});
      root_no->set_scale(set_all, 0.0);
      
      // For both images, two nodes are pinched together and two are stretched apart.
      {
        int i {size_nodes/4};
        int j {size_nodes/4};
        int di {size_nodes/8};
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes/4;
        di = -size_nodes/8;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        i = size_nodes/4;
        j = size_nodes*3/4;
        di = -size_nodes/4;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes*3/4;
        di = size_nodes/4;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        
        i = size_nodes/4;
        j = size_nodes/4;
        di = size_nodes/8;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes/4;
        di = -size_nodes/8;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
        i = size_nodes/4;
        j = size_nodes*3/4;
        di = -size_nodes/4;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes*3/4;
        di = size_nodes/4;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
      }
      
      root_no->relax(step_count, step_size);
      root_fix->relax(step_count, step_size);
      
      {
        Image image {path + input};
        image.set_positions(root_fix, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output_fix, size_output, size_output);
      }
      {
        Image image {path + input};
        image.set_positions(root_no, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output_no, size_output, size_output);
      }
      
      break;
    }
    case Demo::SKEW:
    {
      // Three grids of nodes are created: one with a scale (1.0 by default) and two without.
      // The grid with a scale and one of the ones without both have two upper nodes and two
      // lower nodes offset horizontally to skew the images. The other grid without a scale
      // has its top and bottom edges offset horizontally to skew the image. The resulting
      // images show different deformations.
      
      //input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      const std::string output_fix {"output-fixed-scale.rgb"};
      const std::string output_no {"output-no-scale.rgb"};
      const std::string output_edge {"output-edge.rgb"};
      // path for input and output
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/skew/"};
      
      const int size_nodes {64};
      const int size_output {512};
      const int step_count {100000};
      const double scale {4.0};
      const double step_size {0.05};
      const Vector origin;
      const Vector offset {0.125};
      
      std::set<int> all;
      for(int n {0}; n < size_nodes*size_nodes; ++n) {
        all.insert(n);
      }
      std::shared_ptr<Node> root_no {Node::create_nodes(size_nodes)};
      std::shared_ptr<Node> root_fix {Node::create_nodes(size_nodes)};
      std::shared_ptr<Node> root_edge {Node::create_nodes(size_nodes)};
      {
        const Vector offset {size_nodes/2, size_nodes/2};
        root_no->apply_translation(all, offset);
        root_no->set_scale(all, 0.0);
        root_fix->apply_translation(all, offset);
        root_edge->set_scale(all, 0.0);
        root_edge->apply_translation(all, offset);
      }
      
      // For both images, two nodes at the center of their bottom and top quarters
      // are shifted to left and to right, respectively, to skew the images.
      {
        double i {size_nodes/4};
        double j {size_nodes/4};
        double di {-size_nodes/2};
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes/4;
        di = -size_nodes/2;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        i = size_nodes/4;
        j = size_nodes*3/4;
        di = size_nodes/2;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes*3/4;
        di = size_nodes/2;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {(double)di});
        root_no->set_fixed(i, j);
        
        i = size_nodes/4;
        j = size_nodes/4;
        di = -size_nodes/2;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes/4;
        di = -size_nodes/2;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
        i = size_nodes/4;
        j = size_nodes*3/4;
        di = size_nodes/2;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes*3/4;
        di = size_nodes/2;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {(double)di});
        root_fix->set_fixed(i, j);
        
        std::set<int> bottom;
        std::set<int> top;
        for(int n {0}; n < size_nodes; ++n) {
          bottom.insert(n);
          top.insert((size_nodes - 1)*size_nodes + n);
        }
        di = size_nodes/2;
        root_edge->apply_translation(bottom, Vector {-(double)di});
        root_edge->set_fixed(bottom);
        root_edge->apply_translation(top, Vector {(double)di});
        root_edge->set_fixed(top);
      }
      
      root_no->relax(step_count, step_size);
      root_fix->relax(step_count, step_size);
      root_edge->relax(step_count, step_size);
      
      {
        Image image {path + input};
        image.set_positions(root_no, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output_no, size_output, size_output);
      }
      {
        Image image {path + input};
        image.set_positions(root_fix, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output_fix, size_output, size_output);
      }
      {
        Image image {path + input};
        image.set_positions(root_edge, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output_edge, size_output, size_output);
      }
      
      break;
    }
    case Demo::STRETCH:
    {
      // For both images, the center points of the upper (lower) quadrants are pulled
      // up (down), and toward one another horizontally. The image without a fixed scale
      // is far better able to contract at the top and at the bottom and to expand in the
      // middle. The image with a fixed scale shows greatest deformation around the
      // manipulated points.
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      const std::string output_fix {"output-fixed-scale.rgb"};
      const std::string output_no {"output-no-scale.rgb"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/stretch/"};
      
      const int size_nodes {64};
      const int size_output {512};
      const int step_count {100000};
      const double step_size {0.05};
      const double scale {4.0};
      const Vector origin;
      const Vector offset {0.125};
      
      std::set<int> all;
      for(int n {0}; n < size_nodes*size_nodes; ++n) {
        all.insert(n);
      }
      std::shared_ptr<Node> root_fix {Node::create_nodes(size_nodes)};
      std::shared_ptr<Node> root_no {Node::create_nodes(size_nodes)};
      root_fix->apply_translation(all, Vector {size_nodes/2, size_nodes/2});
      root_no->apply_translation(all, Vector {size_nodes/2, size_nodes/2});
      root_no->set_scale(all, 0.0);
      
      {
        double i {size_nodes/4};
        double j {size_nodes/4};
        const double di {size_nodes/8};
        const double dj {size_nodes/2};
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {di, -dj});
        root_no->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes/4;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {-di, -dj});
        root_no->set_fixed(i, j);
        i = size_nodes/4;
        j = size_nodes*3/4;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {di, dj});
        root_no->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes*3/4;
        root_no->set_position(i, j, root_no->get_position(i, j) + Vector {-di, dj});
        root_no->set_fixed(i, j);
        
        i = size_nodes/4;
        j = size_nodes/4;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {di, -dj});
        root_fix->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes/4;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {-di, -dj});
        root_fix->set_fixed(i, j);
        i = size_nodes/4;
        j = size_nodes*3/4;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {di, dj});
        root_fix->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes*3/4;
        root_fix->set_position(i, j, root_fix->get_position(i, j) + Vector {-di, dj});
        root_fix->set_fixed(i, j);
      }
      
      root_no->relax(step_count, step_size);
      root_fix->relax(step_count, step_size);
      
      {
        Image image {path + input};
        image.set_positions(root_no, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + "output-no-scale.rgb", size_output, size_output);
      }
      {
        Image image {path + input};
        image.set_positions(root_fix, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + "output-fixed-scale.rgb", size_output, size_output);
      }
      
      break;
    }
    case Demo::TWIST:
    {
      // The center points of two diagonally adjacent quadrants are offset to left and
      // to right so that the image is rotated by 90 degrees. The rest of the image
      // eventually follows the offset points.
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      const std::string output {"output-baboon-"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/twist/"};
      
      const int size_nodes {64};
      const int size_output {512};
      const int step_count {100000};
      const double step_size {0.05};
      const double scale {8.0};
      const double sigma {0.75};
      const double gamma {0.5};
      const Vector origin;
      const Vector offset {0.125};
      
      std::set<int> all;
      for(int n {0}; n < size_nodes*size_nodes; ++n) {
        all.insert(n);
      }
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      root->set_scale(all, 0.0);
      
      {
        int i {size_nodes/4};
        int j {size_nodes/4};
        int di {size_nodes/2};
        root->set_position(i, j, root->get_position(i, j) + Vector {(double)di});
        root->set_fixed(i, j);
        i = size_nodes*3/4;
        j = size_nodes*3/4;
        di = -size_nodes/2;
        root->set_position(i, j, root->get_position(i, j) + Vector {(double)di});
        root->set_fixed(i, j);
      }
      
      int n {0};
      int i {0};
      while(n < step_count) {
        root->relax(n - i, step_size);
        Image image {path + input};
        image.set_positions(root, origin, offset);
        image.scale_positions(scale);
        image.write_convolved(path + output + std::to_string(n) + ".rgb", size_output, size_output, sigma, gamma);
        i = n;
        n == 0 ? n++ : n *= 10;
      }
      
      break;
    }
    case Demo::PINCHES:
    {
      // Random points are displaced randomly and fixed. The rest of the image
      // stretches elastically.
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      const std::string output {"output-baboon.rgb"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/pinches/"};
      
      const int size_nodes {128};
      const int size_output {512};
      const int pinch_count {20};
      const int step_count {100000};
      const double step_size {0.05};
      const double pinch_sigma {20.0};
      const double scale {4.0};
      const double sigma {0.75};
      const double gamma {0.5};
      const Vector origin;
      const Vector offset {0.25};
      
      std::set<int> all;
      for(int n {0}; n < size_nodes*size_nodes; ++n) {
        all.insert(n);
      }
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      root->set_scale(all, 0.0);

      std::default_random_engine rng {std::random_device {}()};
      std::uniform_int_distribution<int> discrete {0, size_nodes - 1};
      std::normal_distribution<double> normal {0.0, pinch_sigma};
      std::uniform_real_distribution<double> uniform {0.0, 2.0*M_PI};

      for(int n = 0; n < pinch_count; ++n) {
        const int i = discrete(rng);
        const int j = discrete(rng);
        const double radius = normal(rng);
        const double angle = uniform(rng);
        const Vector displacement {radius*std::cos(angle), radius*std::sin(angle)};
        root->set_position(i, j, root->get_position(i, j) + displacement);
        root->set_fixed(i, j);
      }
      
      root->relax(step_count, step_size);
      
      Image image {path + input};
      image.set_positions(root, origin, offset);
      image.scale_positions(scale);
      
      image.write_convolved(path + output, size_output, size_output, sigma, gamma);
      
      break;
    }
    case Demo::QUENCH:
    {
      // Nodes in a grid are initially offset randomly except the edges that are fixed.
      // The offset nodes are not fixed but are allowed to relax freely. The image becomes
      // smoother and smoother as the greatest displacements are eliminated fastest.
      // Note that instabilities will likely occur if much larger initial offsets are used.
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      const std::string output {"output-baboon-"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/quench/"};
      
      const int size_nodes {128};
      const int size_output {512};
      const int step_count {10000};
      const double displacement_sigma {25.0};
      const double step_size {0.0125};
      const double scale {4.0};
      const double sigma {0.75};
      const double gamma {0.5};
      const Vector origin;
      const Vector offset {0.25};
      
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      
      // edges are fixed
      for(int n {0}; n < size_nodes; ++n) {
        root->set_fixed(n, 0);
        root->set_fixed(n, size_nodes - 1);
        root->set_fixed(0, n);
        root->set_fixed(size_nodes - 1, n);
      }

      std::default_random_engine rng {std::random_device {}()};
      std::normal_distribution<double> normal {0.0, displacement_sigma};
      std::uniform_real_distribution<double> uniform {0.0, 2.0*M_PI};

      // random displacements are applied
      for(int j {1}; j < size_nodes - 1; ++j) {
        for(int i {1}; i < size_nodes - 1; ++i) {
          const double radius {normal(rng)};
          const double angle {uniform(rng)};
          const Vector displacement {radius*std::cos(angle), radius*std::sin(angle)};
          root->set_position(i, j, root->get_position(i, j) + displacement);
        }
      }
      
      int n {0};
      int i {0};
      while(i < step_count) {
        root->relax(n - i, step_size);
        Image image {path + input};
        image.set_positions(root, origin, offset);
        image.scale_positions(scale);
        image.write_convolved(path + output + std::to_string(n) + ".rgb", size_output, size_output, sigma, gamma);
        i = n;
        n == 0 ? n++ : n *= 2;
      }
      
      break;
    }
    case Demo::SPREAD:
    {
      // The nodes corresponding to the eyes of the baboon are moved to the left and to
      // the right and fixed, and the rest of the nodes are left to relax freely resulting
      // in some deformation in response. Note that only the outer edge of the right eye
      // is moved and the interior is left unmoved. The relatively small interior quickly
      // follows the outer edge.
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      // you can use for example the provided mask images here
      // note that these two images should have the same dimensions as the node grid
      const std::string input_left {"input-baboon-left.rgb"};
      const std::string input_right {"input-baboon-right.rgb"};
      const std::string output {"output-"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/spread/"};
      
      const int size_nodes {128};
      const int step_count {100000};
      const int size_output {512};
      const double step_size {0.05};
      const double scale {4.0};
      const double sigma {0.75};
      const double gamma {0.5};
      const Vector displacement {25.0};
      const Vector origin;
      const Vector offset {0.25};
      
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      
      // The left eye is displaced to the left and fixed.
      Image img_eye_left {path + input_left};
      std::set<int> sel_left {img_eye_left.threshold()};
      root->apply_translation(sel_left, -displacement);
      root->set_fixed(sel_left);
      
      // The outer edge of the right eye is displaced to the right and fixed.
      // Note how the interior of the eye follows the outer edge during relaxation.
      Image img_eye_right {path + input_right};
      std::set<int> sel_right {img_eye_right.threshold()};
      root->apply_translation(sel_right, displacement, true);
      root->set_fixed(sel_right, true);
      
      int n {0};
      int i {0};
      Image image {path + input};
      while(i < step_count) {
        root->relax(n - i, step_size, true);
        image.set_positions(root, origin, offset);
        image.scale_positions(scale);
        image.write_convolved(path + output + std::to_string(n) + ".rgb", size_output, size_output, sigma, gamma);
        i = n;
        n == 0 ? n++ : n *= 10;
      }
      
      break;
    }
    case Demo::RESIZE:
    {
      // The nose of the baboon is shrunk and its right eye is enlarged. The image
      // deforms significantly in response. The image tilts somewhat which might be a
      // consequence of conservation of angular momentum or simply some numerical drift
      // (elimination of drift and rotation is used when relaxing but it is not always
      // perfect).
      
      // input and output strings
      // you can use for example a 512x512 baboon image here
      const std::string input {"input-baboon.rgb"};
      // you can use for example the provided mask images here
      const std::string input_nose {"input-baboon-nose.rgb"};
      const std::string input_left {"input-baboon-left.rgb"};
      const std::string input_right {"input-baboon-right.rgb"};
      const std::string output {"output-"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/resize/"};
      
      const int size_nodes {128};
      const int size_output {512};
      const int step_count {100000};
      const double scale_nose {0.05};
      const double scale_left {1.0};
      const double scale_right {20.0};
      const double step_size {0.01};
      const double scale {4.0};
      const Vector origin;
      const Vector offset {0.25};
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      
      // The nose of the baboon is shrunk.
      Image img_nose {path + input_nose};
      std::set<int> sel_nose {img_nose.threshold()};
      root->set_scale(sel_nose, scale_nose);
      
      // The left eye of the baboon is unchanged.
      Image img_left {path + input_left};
      std::set<int> sel_left {img_left.threshold()};
      root->set_scale(sel_left, scale_left);
      
      // The right eye of the baboon is enlarged.
      Image img_right {path + input_right};
      std::set<int> sel_right {img_right.threshold()};
      root->set_scale(sel_right, scale_right);
      
      int n {0};
      int i {0};
      Image image {path + input};
      while(i < step_count) {
        root->relax(n - i, step_size, true);
        image.set_positions(root, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output + std::to_string(n) + ".rgb", size_output, size_output);
        i = n;
        n == 0 ? n++ : n *= 10;
      }
      
      break;
    }
    case Demo::MOSAIC:
    {
      // An 8-by-8 checkerboard is divided into 4-by-4 tiles whose scale and stiffness
      // a sampled randomly. The end result is a distorted checkerboard obviously.
      
      // input and output strings
      // you can use for example the provided checkerboard image (drew myself) here
      const std::string input {"input-checkerboard.rgb"};
      const std::string output {"output-checkerboard.rgb"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/mosaic/"};
      
      const int size_nodes {64};
      const int size_output {512};
      const int size_tiles {16};
      const int step_count {100000};
      const double mean {0.0};
      const double sigma {0.5};
      const double step_size {0.025};
      const double scale {8.0};
      const Vector origin;
      const Vector offset {0.125};
      
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      
      // The image is divided into a grid of tiles each with a random scale and stiffness.
      std::default_random_engine rng {std::random_device {}()};
      std::lognormal_distribution<double> lognormal {mean, sigma};
      for(int v {0}; v < size_nodes/size_tiles; ++v) {
        for(int u {0}; u < size_nodes/size_tiles; ++u) {
          std::set<int> tile;
          for(int j {size_tiles*v}; j < size_tiles*(v + 1); ++j) {
            for(int i {size_tiles*u}; i < size_tiles*(u + 1); ++i) {
              tile.insert(size_nodes*j + i);
            }
          }
          root->set_scale(tile, lognormal(rng));
          root->set_stiffness(tile, lognormal(rng));
        }
      }
      
      Image image {path + input};
      root->relax(step_count, step_size, true);
      image.set_positions(root, origin, offset);
      image.scale_positions(scale);
      image.write_interpolated(path + output, size_output, size_output);
      
      break;
    }
    case Demo::WORLD:
    {
      // A map of the world is deformed according to a grayscale colormap: larger
      // countries are enlarged while smaller countries are shrunk. Unfortunately the
      // relaxation is not very stable for this particular problem so small step sizes
      // must be used and the deformations are not very large.
      
      // input and output strings
      // you can use for example a 1024x1024 world map image here
      const std::string input {"input-world.rgb"};
      // you can use for example the provided mask images here (remember to make sure
      // that the actual input image and the mask images are aligned)
      const std::string input_small {"input-world-small.rgb"};
      const std::string output {"output-world-"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/world/"};
      
      const int size_nodes {512};
      const int size_output {1024};
      const int step_count {100000};
      const double step_size {0.005};
      const double scale {2.0};
      const Vector origin;
      const Vector offset {0.5};
      
      std::shared_ptr<Node> root {Node::create_nodes(size_nodes)};
      for(int n {0}; n < size_nodes; ++n) {
        root->set_fixed(n, 0);
        root->set_fixed(n, size_nodes - 1);
        root->set_fixed(0, n);
        root->set_fixed(size_nodes - 1, n);
      }
      std::set<int> all;
      for(int n {0}; n < size_nodes*size_nodes; ++n) {
        all.insert(n);
      }
      root->set_scale(all, 0.0);
      // Scale of nodes is overwritten based on the grayscale levels in the image.
      // Larger countries become enlarged while smaller countries shrink.
      // parameters for mapping grayscale levels to node scale (from trial and error)
      const double a {10.0};
      const double b {4.0};
      const double c {-0.5};
      Image image_small {path + input_small};
      for(double t {0.01}; t <= 1.01; t += 0.01) {
        std::set<int> selection {image_small.threshold(t)};
        root->set_scale(selection, std::pow(a, b*(t + c)));
      }
      
      int n {0};
      int i {0};
      Image image {path + input};
      while(i < step_count) {
        root->relax(n - i, step_size);
        image.set_positions(root, origin, offset);
        image.scale_positions(scale);
        image.write_interpolated(path + output + std::to_string(n) + ".rgb", size_output, size_output);
        i = n;
        n == 0 ? n++ : n *= 10;
      }
      
      break;
    }
    case Demo::BOTCH:
    {
      // Pictures of ancient male and female statues are enhanced to meet the modern
      // expectations on social media. Note how the poor guy's expanding muscles
      // squeeze his other parts.
      
      // input and output strings
      // you can use for example 1024x1024 male and female statue images here
      const std::string input_female {"input-female.rgb"};
      const std::string input_male {"input-male.rgb"};
      // you can use for example the provided mask images here (remember to make sure
      // that the actual input image and the mask images are aligned)
      const std::string input_female_fix {"input-female-fix.rgb"};
      const std::string input_female_shrink {"input-female-shrink.rgb"};
      const std::string input_female_enlarge {"input-female-enlarge.rgb"};
      const std::string output_female {"output-female-"};
      const std::string input_male_fix {"input-male-fix.rgb"};
      const std::string input_male_enlarge {"input-male-enlarge.rgb"};
      const std::string input_male_enhuge {"input-male-enhuge.rgb"};
      const std::string output_male {"output-male-"};
      const std::string path {"/home/pete/Dropbox/NetBeans/Deformer/demos/botch/"};
      
      int size_nodes {256};
      const int size_output {1024};
      const int step_count {100000};
      const double step_size {0.05};
      const double scale {4.0};
      const Vector origin;
      const Vector offset {0.25};
      
      std::shared_ptr<Node> root_female {Node::create_nodes(size_nodes)};
      {
        // edges are fixed
        Image image_fix {path + input_female_fix};
        std::set<int> sel_fix {image_fix.threshold()};
        root_female->set_fixed(sel_fix);
        double scale {0.0};
        std::set<int> sel_all {image_fix.threshold(-1.0)};
        root_female->set_scale(sel_all, scale);
      }
      {
        // waist is shrunk
        Image image_shrink {path + input_female_shrink};
        double scale {0.5};
        std::set<int> sel_shrink {image_shrink.threshold()};
        root_female->set_scale(sel_shrink, scale);
      }
      {
        // breasts and hips are enlarged
        Image image_enlarge {path + input_female_enlarge};
        double scale {100.0};
        std::set<int> sel_enlarge {image_enlarge.threshold()};
        root_female->set_scale(sel_enlarge, scale);
      }
      
      std::shared_ptr<Node> root_male {Node::create_nodes(size_nodes)};
      {
        // edges are fixed
        Image image_fix {path + input_male_fix};
        std::set<int> sel_fix {image_fix.threshold()};
        root_male->set_fixed(sel_fix);
        double scale {0.0};
        std::set<int> sel_all {image_fix.threshold(-1.0)};
        root_male->set_scale(sel_all, scale);
      }
      {
        // larger muscles are enlarged
        Image image_enlarge {path + input_male_enlarge};
        double scale {10.0};
        std::set<int> sel_enlarge {image_enlarge.threshold()};
        root_male->set_scale(sel_enlarge, scale);
      }
      {
        // private parts are enhuged
        Image image_enhuge {path + input_male_enhuge};
        double scale {1000.0};
        std::set<int> sel_enhuge {image_enhuge.threshold()};
        root_male->set_scale(sel_enhuge, scale);
      }
      
      int n {0};
      int i {0};
      Image image_female {path + input_female};
      Image image_male {path + input_male};
      while(i < step_count) {
        root_female->relax(n - i, step_size);
        image_female.set_positions(root_female, origin, offset);
        image_female.scale_positions(scale);
        image_female.write_interpolated(path + output_female + std::to_string(n) + ".rgb", size_output, size_output);
        
        root_male->relax(n - i, step_size);
        image_male.set_positions(root_male, origin, offset);
        image_male.scale_positions(scale);
        image_male.write_interpolated(path + output_male + std::to_string(n) + ".rgb", size_output, size_output);
        
        i = n;
        n == 0 ? n++ : n *= 10;
      }
      
      break;
    }
  }
}