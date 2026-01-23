#ifdef OVERFLOW_EXCEPT
#ifdef __linux__
#include <cfenv>
#elif _MSC_VER
#include <float.h>
#else
#include <cfenv>
#endif
#endif

#include <iostream>
#include "ThermoFun/ThermoFun.h"
#include "v_detail.h"

bool load_thermodynamic(ThermoFun::ThermoEngine* thermo_engine, double TK, double PPa);

int main( int argc, char* argv[] )
{

#if  defined(OVERFLOW_EXCEPT)
#ifdef __linux__
    feenableexcept (FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW);
#elif _MSC_VER
    _clearfp();
    _controlfp(_controlfp(0, 0) & ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW),
               _MCW_EM);
#else

#endif
#endif
    try {

        std::string thermo_file = "UOx_jsonfun/UOx-fun.json";

        std::unique_ptr<ThermoFun::ThermoEngine> thermo_engine;
        thermo_engine.reset(new ThermoFun::ThermoEngine(thermo_file));

        std::cout << "\n 1 ";
        load_thermodynamic(thermo_engine.get(), 298.15, 100000);

        std::cout << "\n 2 ";
        load_thermodynamic(thermo_engine.get(), 298.15, 100000);

        std::cout << "\n 3 ";
        load_thermodynamic(thermo_engine.get(), 298.15, 100000);
        std::cout << std::endl;
        return 0;
    }
    catch(std::exception& e)
    {
        std::cout  << "std::exception: " << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout  << "unknown exception" << std::endl;
    }
    return -1;
}

bool load_thermodynamic(ThermoFun::ThermoEngine* thermo_engine, double TK, double PPa)
{
    if(!thermo_engine) {
        return false;
    }
    try {
        std::cout << "\nCalc ThermoEngine T: " << TK <<"  P: " << PPa;

        double G0, Vol, S0, H0, Cp0;
        double funT = TK, funP=PPa;  // T in K, P in Pa

        std::vector<std::string> subst_list = {"O(g)", "O2(g)", "O3(g)"};
        //std::vector<std::string> subst_list = thermo_engine->database().getSubstancesList();

        for(const auto& symbol: subst_list)  {
            auto propAl = thermo_engine->thermoPropertiesSubstance(funT, funP, symbol);

            G0 = propAl.gibbs_energy.val;
            Vol= propAl.volume.val*10;
            S0 = propAl.entropy.val;
            H0 = propAl.enthalpy.val;
            Cp0 = propAl.heat_capacity_cp.val;

            std::cout << "\n" << symbol << ";" << floating_point_to_string(G0)
                      << ";" << floating_point_to_string(Vol) << ";" << floating_point_to_string(S0)
                      << ";" << floating_point_to_string(H0) << ";" << floating_point_to_string(Cp0);
        }
    }
    catch(const std::runtime_error& exception)   {
        std::cout << "\nThermoEngine error: " << exception.what();
    }
    return true;
}

