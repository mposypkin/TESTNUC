/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dejongbndsupp.hpp
 * Author: Alexander Usov
 *
 * Created on December 6, 2016, 12:32 PM
 */

#ifndef DEJONGBNDSUPP_HPP
#define DEJONGBNDSUPP_HPP

#include <box/boxutils.hpp>
#include <common/interval.hpp>
#include <cutfact/lbcutfact/boundsupp.hpp>


namespace TESTNUC {

    /**
     * Interval analisis lower bound supplier for Deijong1 function 
     */
    template <class FT> class DejongBoundSupplier : public NUC::BoundSupplier <FT>{
    public:

        /**
         * Constructor
         * @param n problem dimension
         */
        DejongBoundSupplier(int n) : mN(n) {

        }
        
        /**
         * Retrieve 
         * @param box to find the bound
         * @return bound
         */
        FT getBound(const snowgoose::Box<FT>& box) {
            FT lb = 0;
            FT ub = 0;
            for(int i=0; i < mN; i++)
            {
                FT min, max;
                snowgoose::Interval<FT>::sqr(box.mA[i], box.mB[i], &min, &max);
                snowgoose::Interval<FT>::sum(lb, ub, min, max, &lb, &ub);
            }          
            return lb;
        }

    private:
        int mN;
    };

}

#endif /* DEJONGBNDSUPP_HPP */

