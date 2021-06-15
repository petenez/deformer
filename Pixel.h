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
#include <memory>

class Pixel {
private:
    bool i_fix_ = false;
    bool j_fix_ = false;
    double i_;
    double j_;
    double i__;
    double j__;
    double r_ = 0.0;
    double g_ = 0.0;
    double b_ = 0.0;
    std::vector<std::shared_ptr<Pixel>> adjacent_;
    std::vector<std::shared_ptr<Pixel>> diagonal_;
    
public:
    Pixel(int i, int j);
    Pixel(int i, int j, double r, double g, double b);
    const double get_r() const;
    void set_r(const double d);
    const double get_g() const;
    void set_g(const double d);
    const double get_b() const;
    void set_b(const double d);
    const double get_i() const;
    void set_i(const double d);
    const double get_j() const;
    void set_j(const double d);
    const int count_adjacent() const;
    std::shared_ptr<Pixel> get_adjacent(const int i) const;
    void add_adjacent(std::shared_ptr<Pixel> pixel);
    const int count_diagonal() const;
    std::shared_ptr<Pixel> get_diagonal(const int i) const;
    void add_diagonal(std::shared_ptr<Pixel> pixel);
};

#endif /* PIXEL_H */

