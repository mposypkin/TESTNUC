#include <limits>

#include <common/vec.hpp>
#include <cutfact/lbcutfact/recordsupp.hpp>
#include <cutfact/lbcutfact/lbcutfactory.hpp>
#include <cutfact/lipgradfact/lipgradfact.hpp>
#include <cutfact/liphessfact/liphessfact.hpp>
#include <cutfact/compfact/compcutfact.hpp>
#include <bags/bfsbag.hpp>
#include <solver/basesolver.hpp>
#include <applycut/basic/serialcutapp.hpp>
#include <decomp/bisectdecomp.hpp>
#include "apsupp.hpp"



/*
 * 
 */
int main(int argc, char** argv) {

    // Setup problem
    const int n = 2;
    OPTITEST::AluffiPentiniObjective obj;
    snowgoose::Box<double> box(n);
    initboxAP(box);

    // Setup Cut Factory 1
    const double eps = 1e-4;
    const double L = 4;
    NUC::RecordSupplier<double> rs(std::numeric_limits<double>::max());
    APBoundSupp supp;
    NUC::LBCutFactory<double> cf(eps, rs, supp);

    // Setup cut factory 2
    APGradSupp gsupp;
    NUC::LipGradCutFactory<double> lfact(obj, box, gsupp);

    // Setup cut factory 3
    APHessSupp hsupp;
    NUC::LipHessCutFactory<double> hfact(obj, box, hsupp);

    // Setup composite cut factory 
    NUC::CompositeCutFactory<double> compf;
    compf.addFactory(&cf);
    //compf.addFactory(&lfact);
    //compf.addFactory(&hfact);

    // Setup bag of sub problems
    NUC::Sub<double> sub(0, std::numeric_limits<double>::max(), box);
    NUC::BFSBag<double> bag;
    bag.putSub(sub);

    // Setup decomposer
    NUC::BisectDecomp<double> bisdec;

    // Setup cut applicator 
    NUC::SerialCutApplicator<double> cutapp;

    // Setup solver
    NUC::BaseSolver<double> solver(bag, bisdec, compf, cutapp);

    double x[n];
    double xbest[n];
    // Setup sub evaluators
    auto sf = [&](NUC::Sub<double>& s) {
        snowgoose::BoxUtils::getCenter(s.mBox, x);
        double v = obj.func(x);
        if(v < rs.getBound(sub.mBox)) {
            snowgoose::VecUtils::vecCopy(n, x, xbest);
        }
        rs.updateRv(v);
        double r = snowgoose::BoxUtils::radius(s.mBox);
        s.mScore = v - r * L;
    };

    solver.addSubEval(sf);

    int cnt = 0;
    // Setup step watchers
    auto tf = [&](const NUC::Sub<double>& s, const NUC::BaseSolver<double>& slv) {
#if 0        
        std::cout << "Sub: \n";
        std::cout << "Score = " << s.mScore;
        std::cout << "\n Layer = " << s.mLayer;
        std::cout << "\n Box = " << snowgoose::BoxUtils::toString(s.mBox) << "\n";
#endif
        cnt++;
    };
    solver.addStepWatcher(tf);

    // Preset value 
    rs.setRv(-0.352386);
    // Run solver
    solver.solve();

    std::cout << "Best value found: " << rs.getBound(sub.mBox) << "\n";
    std::cout << "In " << cnt << " steps\n";
    std::cout << "best point: " << snowgoose::VecUtils::vecPrint(n, xbest) << "\n";

    return 0;
}

