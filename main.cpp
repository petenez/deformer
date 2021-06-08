/* 
 * File:   main.cpp
 * Author: Petri Hirvonen
 * Created on 26 May, 2021
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Pixel.h"

int main(int argc, char** argv) {
  
  std::string input = "/home/pete/Dropbox/NetBeans/Deformer/input.rgb";
  std::string output = "/home/pete/Dropbox/NetBeans/Deformer/output.rgb";
  
  std::ifstream reader(input, std::ios::binary);
  
  int Nx;
  int Ny;
  reader.read((char*)&Nx, sizeof(int));
  reader.read((char*)&Ny, sizeof(int));
  
  std::vector<Pixel> pixels;
  pixels.reserve(Nx*Ny);
  
  {
    double r;
    double g;
    double b;
    for(int j = 0; j < Ny; ++j) {
      for(int i = 0; i < Nx; ++i) {
        reader.read((char*)&r, sizeof(double));
        reader.read((char*)&g, sizeof(double));
        reader.read((char*)&b, sizeof(double));
        pixels.push_back(Pixel(i, j, r, g, b));
      }
    }
  }
  
  reader.close();
  
  std::ofstream writer(output, std::ios::binary);
  
  writer.write((char*)&Nx, sizeof(int));
  writer.write((char*)&Ny, sizeof(int));
  
  for(int n = 0; n < pixels.size(); ++n) {
    Pixel pixel = pixels[n];
    double r = pixel.r();
    double g = pixel.g();
    double b = pixel.b();
    writer.write((char*)&r, sizeof(double));
    writer.write((char*)&g, sizeof(double));
    writer.write((char*)&b, sizeof(double));
  }
  
  writer.close();

  return 0;
}