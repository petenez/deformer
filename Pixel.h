/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Pixel.h
 * Author: pete
 *
 * Created on May 29, 2021, 8:32 AM
 */

#ifndef PIXEL_H
#define PIXEL_H

#include <vector>

class Pixel {
private:
    bool i_fix_ = false;
    bool j_fix_ = false;
    const int i0_;
    const int j0_;
    double i_;
    double j_;
    const double r_;
    const double g_;
    const double b_;
    std::vector<Pixel> neighbors_;
    
public:
    Pixel(int i, int j, double r, double g, double b);
    const double r() const;
    const double g() const;
    const double b() const;
};

#endif /* PIXEL_H */

