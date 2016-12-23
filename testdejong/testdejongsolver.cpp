/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testdejongsolver.cpp
 * Author: Alexander Usov
 *
 * Created on December 7, 2016, 9:45 AM
 */

#include <limits>
#include <vector>
#include <bags/bfsbag.hpp>
#include <cutfact/lbcutfact/recordsupp.hpp>
#include <cutfact/lbcutfact/lbcutfactory.hpp>
#include <applycut/basic/serialcutapp.hpp>
#include <decomp/bisectdecomp.hpp>
#include <oneobj/contboxconstr/dejong.hpp>
#include <solver/basesolver.hpp>
#include <box/boxutils.hpp>
#include <common/interval.hpp>

#include "dejongbndsupp.hpp"

bool stopper(const NUC::BaseSolver<double>& solver) {
    static int cnt = 0;
    const int maxCnt = 100;
    if (cnt++ > maxCnt) {
        return true;
    } else {
        return false;
    }
}

int main() {
    const int n = 2;
    const double eps = 0.01;

    // Setup problem
    OPTITEST::DejongProblemFactory fact(n, -2, 4);
    COMPI::MPProblem<double> *mpp = fact.getProblem();

    //Setup bag of sub problems
    NUC::Sub<double> sub(0, std::numeric_limits<double>::max(), *(mpp->mBox));
    NUC::BFSBag<double> bag;
    bag.putSub(sub);

    //Setup Cut Factory
    NUC::RecordSupplier<double> rs(std::numeric_limits<double>::max());
    COMPI::Functor<double>* pf = mpp->mObjectives.at(0);
    TESTNUC::DejongBoundSupplier<double> ibs(n);
    NUC::LBCutFactory<double> cf(eps, rs, ibs);

    // Setup decomposer
    NUC::BisectDecomp<double> bisdec;
    // Setup cut applicator 
    NUC::SerialCutApplicator<double> cutapp;
    // Setup solver
    NUC::BaseSolver<double> solver(bag, bisdec, cf, cutapp);

    // Set stopper for solver
    solver.setStopper(stopper);
    //Setup step watchers
    auto tf = [&](const NUC::Sub<double>& s,
            const std::vector<std::shared_ptr <NUC::Cut <double> > >&,
            const std::vector< snowgoose::Box<double> >&,
            const NUC::BaseSolver<double>& slv) {
        std::cout << "Sub: \n";
        std::cout << "Score = " << s.mScore;
        std::cout << "\n Layer = " << s.mLayer;
        std::cout << "\n Record = " << rs.getBound(s.mBox);
        auto bf = [](const NUC::Sub<double>& sub){ 
            std::cout << snowgoose::BoxUtils::toString(sub.mBox) << "  ";
        };
        std::cout << "\n Bag: " << bag.size() << " "; bag.traverse(bf);
        std::cout << "\n Box = " << snowgoose::BoxUtils::toString(s.mBox) << "\n";       
    };
    solver.addStepWatcher(tf);

    double x[n];
    //Setup sub evaluators
    auto sf = [&](NUC::Sub<double>& s) {
        snowgoose::BoxUtils::getCenter(s.mBox, x);
        double v = pf->func(x);
        rs.updateRv(v);
        s.mScore = ibs.getBound(s.mBox);
    };
    solver.addSubEval(sf);

    // Run solver
    solver.solve(); 

    std::cout << "Best value found : " << rs.getBound(sub.mBox) << "\n";

    auto ttf = [] (const NUC::Sub<double>& s) {
        std::cout << "Sub: ";
        std::cout << "Score = " << s.mScore;
        std::cout << ", Layer = " << s.mLayer;
        std::cout << ", Box = " << snowgoose::BoxUtils::toString(s.mBox) << "\n";
    };

    solver.getBag().traverse(ttf);
    
    return 0;
}

