/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Image.h
 * Author: pete
 *
 * Created on June 13, 2021, 10:29 PM
 */

#ifndef IMAGE_H
#define IMAGE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "Pixel.h"

class Image {
private:
  int Nx_;
  int Ny_;
  std::vector<std::shared_ptr<Pixel>> pixels_;

public:
  Image(std::string input);
  Image(const int Nx, const int Ny);
  void set_neighbors();
  const int get_Nx() const;
  const int get_Ny() const;
  void deform();
  void read(Image image);
  void normalize();
  void write(std::string output);
};

#endif /* IMAGE_H */

