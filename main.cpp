/* 
 * File:   main.cpp
 * Author: Petri Hirvonen
 * Created on 26 May, 2021
 */

#include <string>

#include "Pixel.h"
#include "Image.h"

int main(int argc, char** argv) {
  
  std::string input = "/home/pete/Dropbox/NetBeans/Deformer/input.rgb";
  std::string output = "/home/pete/Dropbox/NetBeans/Deformer/output.rgb";
  
  Image deformed(input);
  Image result(deformed.get_Nx(), deformed.get_Ny());
  
  deformed.deform();
  
  result.read(deformed);
  
  result.write(output);

  return 0;
}