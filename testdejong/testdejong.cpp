/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testintervalcutfactory.cpp
 * Author: Alexander Usov
 *
 * Created on December 6, 2016, 2:50 PM
 */

#include <iostream>
#include <cmath>
#include "dejongbndsupp.hpp"
#include "box/box.hpp"
#include "common/sgerrcheck.hpp"

#define EPSILON 0.001

int is_equal(double x, double y) 
{
    return std::fabs(x - y) < EPSILON ? 1 : 0;
}

int main() {
    const int n = 2;
    TESTNUC::DejongBoundSupplier<double> ibs(n);
    snowgoose::Box<double> bx(n);
    for(int i=0; i<n; ++i)
    {
        bx.mA[i]=-2;
        bx.mA[i]=4;
    }
    SG_ASSERT(is_equal(ibs.getBound(bx), 0.0));    
    return 0;
}
