/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testackley1.cpp
 * Author: Alexander Usov
 *
 * Created on December 10, 2016, 2:50 PM
 */

#include <iostream>
#include <cmath>
#include "ackley1bndsupp.hpp"
#include "box/box.hpp"
#include "common/sgerrcheck.hpp"

#define EPSILON 0.001


int main() {
    
    const int n = 2;
    TESTNUC::Ackley1BoundSupplier<double> ibs(n);
    snowgoose::Box<double> bx(n), bx2(n);
    bx.mA[0]=0.5;
    bx.mB[0]=1.5;
    bx.mA[1]=-0.5;
    bx.mB[1]=0.5;
    
    bx2.mA[0]=0.9;
    bx2.mB[0]=1.1;
    bx2.mA[1]=-0.1;
    bx2.mB[1]=0.1;
    
    SG_ASSERT(ibs.getBound(bx) <= ibs.getBound(bx2) <= 2.4); 

    return 0;
}
