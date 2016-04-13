/* 
 * File:   apsupp.hpp
 * Author: medved
 *
 * Created on April 12, 2016, 8:10 PM
 */

#ifndef HYPERSUPP_HPP
#define HYPERSUPP_HPP

#include <common/interval.hpp>
#include <oneobj/contboxconstr/hyper.hpp>
#include <cutfact/lbcutfact/boundsupp.hpp>
#include <cutfact/lbcutfact/lpzbndsupp.hpp>
#include <cutfact/lipgradfact/lpzgradsupp.hpp>
#include <cutfact/liphessfact/lpzhesssupp.hpp>

class HyperBoundSupp : public NUC::BoundSupplier <double> {
public:

    double getBound(const snowgoose::Box<double>& box) {
        double a = box.mA[0];
        double b = box.mB[0];
        double c, d;
        snowgoose::Interval<double>::degree(2, a, b, &c, &d);
        a = box.mA[1];
        b = box.mB[1];
        double e, f;
        snowgoose::Interval<double>::degree(2, a, b, &e, &f);
        snowgoose::Interval<double>::diff(c, d, e, f, &a, &b);
        return a;
    }
};

class HyperGradSupp : public NUC::LpzGradSupp <double> {
public:

    void getLpzConst(const snowgoose::Box<double>& box, double* lpzv) {
        lpzv[0] = 2;
        lpzv[1] = 2;
    }

};

class HyperHessSupp : public NUC::LpzHessSupp <double> {
public:

    void getLpzConst(const snowgoose::Box<double>& box, double* lpzv) {
        lpzv[0] = 0;
        lpzv[1] = 0;
        lpzv[2] = 0;
        lpzv[3] = 0;
    }

};

void initboxHyper(snowgoose::Box<double>& box) {
    const int n = box.mDim;
    for (int i = 0; i < n; i++) {
        //box.mA[i] = 1;
        box.mA[i] = -10;
        box.mB[i] = 10;
        //box.mA[i] = -0.5;
        //box.mB[i] = 0.5;
    }
}

#endif /* APSUPP_HPP */

