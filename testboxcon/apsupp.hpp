/* 
 * File:   apsupp.hpp
 * Author: medved
 *
 * Created on April 12, 2016, 8:10 PM
 */

#ifndef APSUPP_HPP
#define	APSUPP_HPP

#include <common/interval.hpp>
#include <oneobj/contboxconstr/alufpent.hpp>
#include <cutfact/lbcutfact/boundsupp.hpp>
#include <cutfact/lbcutfact/lpzbndsupp.hpp>
#include <cutfact/lipgradfact/lpzgradsupp.hpp>
#include <cutfact/liphessfact/lpzhesssupp.hpp>

class APBoundSupp : public NUC::BoundSupplier <double> {
public:
    double getBound(const snowgoose::Box<double>& box) {
        double a = box.mA[0];
        double b = box.mB[0];
        double c, d;
        snowgoose::Interval<double>::degree(4, a, b, &c, &d);
        c *= 0.25;
        d *= 0.25;
        double e,f;
        snowgoose::Interval<double>::degree(2, a, b, &e, &f);
        e *= 0.5;
        f *= 0.5;
        snowgoose::Interval<double>::diff(c, d, e, f, &c, &d);
        snowgoose::Interval<double>::sum(0.1 * a, 0.1 * b, c, d, &c, &d);
        a = box.mA[1];
        b = box.mB[1];        
        snowgoose::Interval<double>::degree(2, a, b, &e, &f);
        e *= 0.5;
        f *= 0.5;
        snowgoose::Interval<double>::sum(c, d, e, f, &c, &d);
        return c;
    }
};

class APGradSupp : public NUC::LpzGradSupp <double> {
public:

    void getLpzConst(const snowgoose::Box<double>& box, double* lpzv) {
        double a = box.mA[0];
        double b = box.mB[0];
        snowgoose::Interval<double>::degree(2, a, b, &a, &b);
        a *= 3;
        b *= 3;
        a -= 1;
        b -= 1;
        lpzv[0] = SGMAX(SGABS(a), SGABS(b));
        lpzv[1] = 1;
    }

};

class APHessSupp : public NUC::LpzHessSupp <double> {
public:

    void getLpzConst(const snowgoose::Box<double>& box, double* lpzv) {
        double a = box.mA[0];
        double b = box.mB[0];
        lpzv[0] = 6 * SGMAX(SGABS(a), SGABS(b));
        lpzv[1] = 0;
        lpzv[2] = 0;
        lpzv[3] = 0;
    }

};

void initboxAP(snowgoose::Box<double>& box) {
    const int n = box.mDim;
    for (int i = 0; i < n; i++) {
        //box.mA[i] = 1;
        box.mA[i] = -10;
        box.mB[i] = 10;
        //box.mA[i] = -0.5;
        //box.mB[i] = 0.5;
    }
}

#endif	/* APSUPP_HPP */

