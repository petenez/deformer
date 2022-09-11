/* 
 * File:   pngreader.cpp
 * Author: Petri Hirvonen
 * Created on 26 May 2021
 * 
 * compile:
    1) source Spack: . ~/Software/spack/spack/share/spack/setup-env.sh
    2) load PNGwriter: spack load pngwriter
    3) compile: g++ pngreader.cpp -o pngreader -std=c++17 `freetype-config --cflags` -I/usr/local/include -L/usr/local/lib -lPNGwriter -lpng -lz -lfreetype
 *
 * use: do 1) and 2), and run: ./pngreader input.png output.rgb
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <pngwriter.h>

int main(int argc, char** argv) {

  argc--;
  if(argc != 2) {
    std::cout << "Error: Invalid arguments!" << std::endl;
    return 1;
  }

  std::string input(argv[1]);
  std::string output(argv[2]);
  
  pngwriter image;
  image.readfromfile(argv[1]);
  
  std::ofstream writer(output, std::ios::binary);
  if(!writer) {
    std::cout << "Error: Could not open output!" << std::endl;
    return 1;
  }
  
  std::string ext_in = input;
  while(ext_in.find('.') != std::string::npos) {
    ext_in = ext_in.substr(ext_in.find('.') + 1);
  }
  if(ext_in.compare("png") != 0 && ext_in.compare("PNG") != 0) {
    std::cout << "Error: The input file extension must be 'png' (or 'PNG')!" << std::endl;
    return 1;
  }
  
  std::string ext_out = output;
  while(ext_out.find('.') != std::string::npos) {
    ext_out = ext_out.substr(ext_out.find('.') + 1);
  }
  if(ext_out.compare("rgb") != 0 && ext_out.compare("hsv") != 0) {
    std::cout << "Error: The output file extension must be either 'rgb' or 'hsv'!" << std::endl;
    return 1;
  }
  
  int Nx = image.getwidth();
  int Ny = image.getheight();
  writer.write((char*)&Nx, sizeof(int));
  writer.write((char*)&Ny, sizeof(int));
  
  double a;
  double b;
  double c;
  
  if(ext_out.compare("rgb") == 0) {
    for(int j = 1; j <= Ny; ++j) {
      for(int i = 1; i <= Nx; ++i) {
        a = image.dread(i, j, 1);
        b = image.dread(i, j, 2);
        c = image.dread(i, j, 3);
        writer.write((char*)&a, sizeof(double));
        writer.write((char*)&b, sizeof(double));
        writer.write((char*)&c, sizeof(double));
      }
    }
  }
  else if(ext_out.compare("hsv") == 0) {
    for(int j = 1; j <= Ny; ++j) {
      for(int i = 1; i <= Nx; ++i) {
        a = image.dreadHSV(i, j, 1);
        b = image.dreadHSV(i, j, 2);
        c = image.dreadHSV(i, j, 3);
        writer.write((char*)&a, sizeof(double));
        writer.write((char*)&b, sizeof(double));
        writer.write((char*)&c, sizeof(double));
      }
    }
  }
  
  writer.close();

  return 0;
}
