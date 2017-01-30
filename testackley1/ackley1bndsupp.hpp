/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ackley1bndsupp.hpp
 * Author: Alexander Usov
 *
 * Created on December 10, 2016, 8:16 AM
 */

#ifndef ACKLEY1BNDSUPP_HPP
#define ACKLEY1BNDSUPP_HPP

#include <math.h>
#include <box/boxutils.hpp>
#include <interval/interval_air.hpp>
#include <cutfact/lbcutfact/boundsupp.hpp>

using namespace snowgoose::interval;



namespace TESTNUC {

    /**
     * Interval analisis lower bound supplier for Ackley1 function 
     */
    template <class FT> class Ackley1BoundSupplier : public NUC::BoundSupplier <FT>{
    public:

        /**
         * Constructor
         * @param n problem dimension
         */
        Ackley1BoundSupplier(int n) : mN(n) {

        }
        
        /**
         * Retrieve 
         * @param box to find the bound
         * @return bound
         */
        FT getBound(const snowgoose::Box<FT>& box) {
            Interval<FT> t1(0), t2(0);
            for(int i=0; i < mN; i++)
            {
                Interval<FT> t(box.mA[i], box.mB[i]);
                t1 = t1 + sqr(t);
                t2 = t2 + cos(2*M_PI*t);
            }
            auto interval = -20.0 * exp(-0.2 * sqrt(t1/mN)) - exp(t2/mN) + 20 + M_E;
            return interval.lb();
        }

    private:
        int mN;
    };

}

#endif /* ACKLEY1BNDSUPP_HPP */

