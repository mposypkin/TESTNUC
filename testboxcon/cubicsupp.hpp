/* 
 * File:   apsupp.hpp
 * Author: medved
 *
 * Created on April 12, 2016, 8:10 PM
 */

#ifndef CUBICSUPP_HPP
#define CUBICSUPP_HPP

#include <common/interval.hpp>
#include <oneobj/contboxconstr/cubic.hpp>
#include <cutfact/lbcutfact/boundsupp.hpp>
#include <cutfact/lbcutfact/lpzbndsupp.hpp>
#include <cutfact/lipgradfact/lpzgradsupp.hpp>
#include <cutfact/liphessfact/lpzhesssupp.hpp>

class CubicBoundSupp : public NUC::BoundSupplier <double> {
public:

    double getBound(const snowgoose::Box<double>& box) {
        int n = box.mDim;
        double e = 0, f = 0;
        for (int i = 0; i < n; i++) {
            double a = box.mA[i];
            double b = box.mB[i];
            double c, d;
            snowgoose::Interval<double>::degree(3, a, b, &c, &d);
            snowgoose::Interval<double>::sum(e, f, c, d, &e, &f);
        }
        return e;
    }

};

class CubicGradSupp : public NUC::LpzGradSupp <double> {
public:

    void getLpzConst(const snowgoose::Box<double>& box, double* lpzv) {
        int n = box.mDim;
        for (int i = 0; i < n; i++) {
            lpzv[i] = 6 * SGMAX(SGABS(box.mA[i]), SGABS(box.mB[i]));
        }
    }

};

class CubicHessSupp : public NUC::LpzHessSupp <double> {
public:

    void getLpzConst(const snowgoose::Box<double>& box, double* lpzv) {
        int n = box.mDim;
        for (int i = 0; i < n; i++) {
            lpzv[i] = 6;
        }
    }

};

void initboxCubic(snowgoose::Box<double>& box) {
    const int n = box.mDim;
    for (int i = 0; i < n; i++) {
        //box.mA[i] = 1;
        box.mA[i] = -10;
        box.mB[i] = 10;
        //box.mA[i] = -0.5;
        //box.mB[i] = 0.5;
    }
}

#endif /* CUBICSUPP_HPP */

