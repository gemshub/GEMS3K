#ifndef SOLMODCALC_H
#define SOLMODCALC_H

#include <map>
#include "s_solmod.h"

/// Addition parameters for specific models,
/// could be add to main SolutionData
struct AddSolutionData {

    double *arZ;
    double *arM;
    double *ardenW;
    double *arepsW;
    double *arG0;
    double *arFWGT;
    double *arX;
};

class SolModCalc
{

public:

    /// The constructor
    SolModCalc(long int k, long int jb, SolutionData& sd, const AddSolutionData& addsd);

    /// Code of the mixing model
    char modCode() const {
      return mod_code;
    }

    /// Call for calculation of temperature and pressure correction
    void SolModParPT(char ModCode);
    /// Call for calculation of activity coefficients
    void SolModActCoeff(char ModCode);
    /// Call for calculation of bulk phase excess properties
    std::map<std::string, double> SolModExcessProp(char ModCode);
    /// Call for calculation of bulk phase ideal mixing properties
    std::map<std::string, double> SolModIdealProp(char ModCode);
    /// Call for retrieving bulk phase Darken quadratic terms  (!!! not implemented)
    std::map<std::string, double> SolModDarkenProp(char ModCode);
    /// Call for retrieving bulk phase standard state terms
    std::map<std::string, double> SolModStandProp(char ModCode);

    // not used in gems
    //virtual long int PureSpecies();

    /// Updates P and T in TSolMod if those have changed
    void UpdatePT( double T_k, double P_bar)
    {
        if(solmod_task) {
            solmod_task->UpdatePT(T_k, P_bar);
        }
    }

    /// Copy activity coefficients into provided array lngamma
    void Get_lnGamma(double* lngamma);
    /// Get activity coefficients into component map
    std::map<std::string, double> GetlnGamma();

    /// Set species (end member) mole fractions from provided array aWx -> dc_num
    void Set_MoleFractionsWx(double* aWx);
    /// Set species (end member) mole fractions from component map
    void SetMoleFractionsWx(const std::map<std::string, double>& awx_map, double defwx=0.);

    /// Writing input structure TSolMod to json format file
    void to_json_file(const std::string& path)
    {
        if(solmod_task) {
            solmod_task->to_json_file(path);
        }
    }
    /// Trace writing arrays TSolMod to keyvalue format file
    void to_text_file(const std::string& path, bool append=false)
    {
        if(solmod_task) {
            solmod_task->to_text_file(path, append);
        }
    }

protected:

    /// Model name (posible add to TSolMod structure for test output)
    std::string model_name;
    /// Code of the mixing model
    char mod_code;

    /// Name of phase
    std::string phase_name;
    /// Index of phase in IPM problem
    long int phase_ndx;
    /// Index first of DCs included into phase
    long int dc_ndx;
    /// Number of DCs included into phase
    long int dc_num;
    /// Names of DCs included into phase
    std::vector<std::string> dc_names;
    /// Mole fractions Wx of DC in multi-component phases -> dc_num
    double *arWx;

    /// TSolMod description
    std::shared_ptr<TSolMod> solmod_task;

    void SolMod_create(SolutionData &sd, const AddSolutionData &addsd);
    bool check_mode(char ModCode);
    std::map<std::string, double> property2map(double* dcs_size_array);
    void map2property(const std::map<std::string, double>& dsc_name_map, double* dcs_size_array, double def_value);
};

#endif // SOLMODCALC_H
