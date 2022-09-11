/* 
 * File:   pngwriter.cpp
 * Author: Petri Hirvonen
 * Created on 26 May 2021
 * 
 * compile:
    1) source Spack: . ~/Software/spack/spack/share/spack/setup-env.sh
    2) load PNGwriter: spack load pngwriter
    3) compile: g++ pngwriter.cpp -o pngwriter -std=c++17 `freetype-config --cflags` -I/usr/local/include -L/usr/local/lib -lPNGwriter -lpng -lz -lfreetype
 *
 * use: do 1) and 2), and run: ./pngwriter input.rgb output.png [magnification]
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <pngwriter.h>

int main(int argc, char** argv) {

  argc--;
  if(argc < 2 || argc > 3) {
    std::cout << "Error: Invalid arguments!" << std::endl;
    return 1;
  }
  
  std::string input(argv[1]);
  std::string output(argv[2]);
  int magnification = argc == 3 ? atoi(argv[3]) : 1;
  
  std::ifstream reader(input, std::ios::binary);
  if(!reader) {
    std::cout << "Error: Could not open input!" << std::endl;
    return 1;
  }
  
  std::string ext_in = input;
  while(ext_in.find('.') != std::string::npos) {
    ext_in = ext_in.substr(ext_in.find('.') + 1);
  }
  if(ext_in.compare("rgb") != 0 && ext_in.compare("hsv") != 0) {
    std::cout << "Error: The input file extension must be either 'rgb' or 'hsv'!" << std::endl;
    return 1;
  }
  
  std::string ext_out = output;
  while(ext_out.find('.') != std::string::npos) {
    ext_out = ext_out.substr(ext_out.find('.') + 1);
  }
  if(ext_out.compare("png") != 0 && ext_out.compare("PNG") != 0) {
    std::cout << "Error: The output file extension must be 'png' (or 'PNG')!" << std::endl;
    return 1;
  }
  
  if(magnification < 1) {
    std::cout << "Error: Magnification must be greater or equal to 1!" << std::endl;
    return 1;
  }
  
  int Nx;
  int Ny;
  reader.read((char*)&Nx, sizeof(int));
  reader.read((char*)&Ny, sizeof(int));
  
  double a;
  double b;
  double c;
  pngwriter writer(magnification*Nx, magnification*Ny, 1, output.c_str());
  
  if(ext_in.compare("rgb") == 0) {
    for(int j = 0; j < Ny; ++j) {
      for(int i = 0; i < Nx; ++i) {
        reader.read((char*)&a, sizeof(double));
        reader.read((char*)&b, sizeof(double));
        reader.read((char*)&c, sizeof(double));
        for(int v = 0; v < magnification; ++v) {
          for(int u = 0; u < magnification; ++u) {
            int i_ = magnification*i + u + 1;
            int j_ = magnification*j + v + 1;
            writer.plot(i_, j_, a, b, c);
          }
        }
      }
    }
  }
  else if(ext_in.compare("hsv") == 0) {
    for(int j = 0; j < Ny; ++j) {
      for(int i = 0; i < Nx; ++i) {
        reader.read((char*)&a, sizeof(double));
        reader.read((char*)&b, sizeof(double));
        reader.read((char*)&c, sizeof(double));
        for(int v = 0; v < magnification; ++v) {
          for(int u = 0; u < magnification; ++u) {
            int i_ = magnification*i + u + 1;
            int j_ = magnification*j + v + 1;
            writer.plotHSV(i_, j_, a, b, c);
          }
        }
      }
    }
  }
  
  reader.close();
  writer.close();

  return 0;
}
