// Tests the AOP/SOP Optima-based equilibrium solver modes end-to-end via
// TNode::GEM_run() (NEED_GEM_AOP/SOP, not by calling TMultiBase methods
// directly - see NODECODECH in databr.h).
//
// Two parts:
//  1. Plain AOP solve of the untitrated tools/Cu-Pourbaix base recipe,
//     checked against GEMS3K's own native AIA reference answer
//     (pH=4.50059, Eh=0.951744).
//  2. A 5-point pH/Eh sweep driven through Set_pH_target()/Set_Eh_target()
//     + NEED_GEM_AOP, exercising the generalized EqControlCondition
//     mechanism (ipm_optima.h/.cpp).

#include "node.h"

#ifdef USE_OPTIMA_SOLVER

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    std::string lst = (argc >= 2) ? argv[1] : "Cu-Pourbaix/Cu-dat.lst";
    int nFail = 0;

    // --- Part 1: plain AOP vs. native AIA reference -------------------
    {
        TNode nodeNative;
        if (nodeNative.GEM_init(lst.c_str())) {
            std::cerr << "Error reading GEMS3K input files: " << lst << std::endl;
            return 1;
        }
        nodeNative.pCNode()->NodeStatusCH = NEED_GEM_AIA;
        long nativeCode = nodeNative.GEM_run(false);
        std::cout << "Native AIA:  status=" << nativeCode
                  << " pH=" << nodeNative.Get_pH() << " Eh=" << nodeNative.Get_Eh() << std::endl;

        TNode nodeOptima;
        if (nodeOptima.GEM_init(lst.c_str())) {
            std::cerr << "Error reading GEMS3K input files: " << lst << std::endl;
            return 1;
        }
        nodeOptima.pCNode()->NodeStatusCH = NEED_GEM_AOP;
        long optimaCode = nodeOptima.GEM_run(false);
        std::cout << "Optima AOP:  status=" << optimaCode
                  << " pH=" << nodeOptima.Get_pH() << " Eh=" << nodeOptima.Get_Eh() << std::endl;
        if (optimaCode != OK_GEM_AOP)
            std::cout << "  [" << nodeOptima.code_error_IPM() << ": " << nodeOptima.description_error_IPM() << "]" << std::endl;

        bool ok1 = (optimaCode == OK_GEM_AOP)
                && std::fabs(nodeOptima.Get_pH() - nodeNative.Get_pH()) < 1e-3
                && std::fabs(nodeOptima.Get_Eh() - nodeNative.Get_Eh()) < 1e-3;

        // Same recipe again, via SOP warm-started from the AOP result above
        // (reusing the in-memory pm.Y[] TNode still holds from the AOP call).
        nodeOptima.pCNode()->NodeStatusCH = NEED_GEM_SOP;
        long sopCode = nodeOptima.GEM_run(false);
        std::cout << "Optima SOP:  status=" << sopCode
                  << " pH=" << nodeOptima.Get_pH() << " Eh=" << nodeOptima.Get_Eh() << std::endl;
        bool ok2 = (sopCode == OK_GEM_SOP)
                && std::fabs(nodeOptima.Get_pH() - nodeNative.Get_pH()) < 1e-3
                && std::fabs(nodeOptima.Get_Eh() - nodeNative.Get_Eh()) < 1e-3;

        std::cout << "PART1 (plain AOP/SOP vs. native AIA): " << ((ok1 && ok2) ? "OK" : "FAIL") << std::endl;
        if (!(ok1 && ok2)) ++nFail;
    }

    // --- Part 2: pH/Eh sweep via generalized control conditions --------
    {
        TNode node;
        if (node.GEM_init(lst.c_str())) {
            std::cerr << "Error reading GEMS3K input files: " << lst << std::endl;
            return 1;
        }
        DATACH* dCH = node.pCSD();
        DATABR* dBR = node.pCNode();
        std::vector<double> base_bIC(dCH->nICb);
        for (long i = 0; i < dCH->nICb; ++i)
            base_bIC[i] = dBR->bIC[i];

        struct TargetPoint { double pH, Eh; };
        std::vector<TargetPoint> targets = {
            {4.0, 0.6}, {7.0, 0.3}, {9.0, -0.2}, {6.0, -0.4}, {11.0, 0.0},
        };

        int nOk = 0;
        for (const auto& t : targets) {
            for (std::size_t i = 0; i < base_bIC.size(); ++i)
                dBR->bIC[i] = base_bIC[i];

            node.Clear_ControlConditions();
            node.Set_pH_target(t.pH);
            node.Set_Eh_target(t.Eh);
            dBR->NodeStatusCH = NEED_GEM_AOP;
            long code = node.GEM_run(false);

            bool ok = (code == OK_GEM_AOP);
            std::cout << "target(pH=" << t.pH << ",Eh=" << t.Eh << "): "
                      << (ok ? "OK  " : "FAIL")
                      << " status=" << code
                      << " pH=" << node.Get_pH() << " Eh=" << node.Get_Eh()
                      << " xiH=" << node.Get_ControlCondition_titrant("pH")
                      << " xiE=" << node.Get_ControlCondition_titrant("Eh");
            if (!ok)
                std::cout << " [" << node.code_error_IPM() << ": " << node.description_error_IPM() << "]";
            std::cout << std::endl;
            if (ok) ++nOk;
        }
        std::cout << "PART2: " << nOk << "/" << targets.size()
                  << " converged (GEMS3K chemistry, generalized Optima control-condition solve)" << std::endl;
        nFail += static_cast<int>(targets.size()) - nOk;
    }

    return nFail;
}

#else

int main()
{
    std::cerr << "optima_test: GEMS3K was not built with "
                 "-DUSE_OPTIMA_SOLVER=ON; nothing to run."
              << std::endl;
    return 1;
}

#endif
