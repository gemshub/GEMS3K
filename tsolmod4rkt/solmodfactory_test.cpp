//--------------------------------------------------------------------
//
/// Demo test of usage of the of subset of TMulti class, configuration,
/// and related functions
//
// Copyright (C) S.Dmytriyeva, D.Kulik
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>
//-------------------------------------------------------------------

#include <iostream>
#include <filesystem>
#include <regex>
#include "jsonconfig.h"
#include "solmodfactory.h"
#include "v_service.h"
namespace fs = std::filesystem;


/// KeyValueFile class used to reading key-value file structure
class KeyValueFile
{
public:
    explicit KeyValueFile( const std::string& path ):
      file_path(path), key_regexp("<([^>]*)>"), comment_regexp("#([^\\\n]*)\\\n")
    {}

    /// Destructor
    virtual ~KeyValueFile()
    {}

    /// Checks if the given path corresponds to an existing file.
    bool exist() const
    {
        fs::path ps(file_path);
        return fs::exists(ps);
    }

    /// Load file to internal structure
    bool load_all()
    {
        file_data.clear();
        if( !exist() )  {
            std::cerr <<  "Error: Trying read not existing file  " <<  file_path << "\n";
            return false;
        }
        auto fdata = read_ascii_file( file_path );

        // remove comments
        if( !comment_regexp.empty() ) {
            fdata = regexp_replace( fdata, comment_regexp, " " );
        }
        auto headers = regexp_extract( fdata, key_regexp );
        auto datas = regexp_split( fdata, key_regexp );

        auto value = datas.begin()+1;
        for( const auto& key: headers ) {
            //std::cout << key << "\n" <<  *value << std::endl;
            if( value < datas.end() )  {
                file_data[key] = *value;
                value++;
            }
        }
        return true;
    }

    /// Compare with template
    template < class T >
    std::string compare_to(const std::string& key, const std::vector<T>& data)
    {
        std::string diff_str;
        auto templ_str = file_data["<"+key+">"];
        auto templ_list = regexp_split(templ_str, "\\s+");

        if(templ_list.size() != data.size() ) {
           diff_str = key + " - the number of values is different: ";
           diff_str += std::to_string(data.size()) + " - " + std::to_string(templ_list.size());
           return "\n     " + diff_str;
        }

        diff_str = diff_arrays(data, templ_list);
        if(!diff_str.empty()) {
          diff_str = "\n     " + key + diff_str;
        }
        return diff_str;
    }

    /// Compare with template
    template < class T >
    std::string compare_subset(const std::string& key, const std::vector<T>& data, size_t from, size_t to)
    {
        std::string diff_str;
        auto templ_str = file_data["<"+key+">"];
        auto templ_list = regexp_split(templ_str, "\\s+");

        if(templ_list.size() < to && !data.empty()) {
           diff_str = key + " - the number of values is different: ";
           diff_str += std::to_string(to) + " - " + std::to_string(templ_list.size());
           return "\n     " + diff_str;
        }

        std::vector<std::string> templ_list_phase{templ_list.begin()+from, templ_list.begin()+to};
        diff_str = diff_arrays(data, templ_list_phase);
        if(!diff_str.empty()) {
          diff_str = "\n     " + key + diff_str;
        }
        return diff_str;
    }

    /// Set up difference value comparing the float-type numbers
    void setEpsilon( double eps )
    {
        epsilon = eps;
    }

private:
    /// File location
    std::string file_path;
    /// Regular expression for key
    std::string key_regexp;
    /// Regular expression for comments or empty
    std::string comment_regexp;
    /// Epsilon/difference value comparing the float-type numbers
    double  epsilon = 1e-14;

    /// Internal structure of file data
    std::map<std::string,std::string> file_data;


    /// Compare with template if size eq
    template < class T >
    std::string diff_arrays(const std::vector<T>& data, const std::vector<std::string>& templ_list)
    {
        std::string diff_str;
        T val;
        for( size_t ii=0; ii<data.size(); ++ii) {
           if( is(val, templ_list[ii]) ) {
                if( !approximatelyEqual( data[ii], val, epsilon ) ) {
                    diff_str += "\n     " + std::to_string(ii) + " - different value into template: ";
                    diff_str += floating_point_to_string(data[ii]) + " - " + floating_point_to_string(val);
                }
           }
           else {
               diff_str += "\n     " + std::to_string(ii) + " - undefined value into template: " + templ_list[ii];
           }
        }
        return diff_str;
    }

    // Read whole ASCII file into string.
    std::string read_ascii_file( const std::string& path )
    {
        std::ifstream t(path);
        std::stringstream buffer;
        buffer << t.rdbuf();
        return buffer.str();
    }
    //  Function that can be used to split text using regexp.
    std::vector<std::string> regexp_split(const std::string& str, std::string rgx_str)
    {
        std::vector<std::string> lst;
        if( !str.empty() ) {
            std::regex rgx(rgx_str);
            std::sregex_token_iterator iter(str.begin(), str.end(), rgx, -1);
            std::sregex_token_iterator end;

            while (iter != end) {
                lst.push_back(*iter);
                trim(lst.back());
                ++iter;
            }
        }
        return lst;
    }
    //  Function that can be used to extract tokens using regexp.
    std::vector<std::string> regexp_extract(const std::string& str, std::string rgx_str)
    {
      std::vector<std::string> lst;
      std::regex rgx(rgx_str);
      std::sregex_token_iterator iter(str.begin(), str.end(), rgx, 0);
      std::sregex_token_iterator end;

      while (iter != end)  {
        lst.push_back(*iter);
        trim(lst.back());
        ++iter;
      }
      return lst;
    }
    //  Function that can be used to replace text using regex.
    std::string regexp_replace(const std::string& instr, const std::string& rgx_str, const std::string& replacement )
    {
       std::regex re(rgx_str);
       std::string output_str = std::regex_replace(instr, re, replacement);
       return output_str;
    }
};

int test_SolModFactory_API(SolModFactory& task, const std::string& template_file_path, const std::string& ipmfiles_lst_name );
int test_SolModEngine_API(SolModFactory& task, const std::string& template_file_path, const std::string& ipmfiles_lst_name );

//The simplest case: data exchange using disk files only
int main(int argc, char* argv[])
{
    try{
        // Analyzing command line arguments  (with defaults)
        std::string input_system_file_list_name = "Thermo-time-all/series1-dat.lst";
        if (argc >= 2 ) {
            input_system_file_list_name = argv[1];
        }

        // Load templalte data
        auto template_file = input_system_file_list_name;
        replace(template_file, "-dat.lst", "-AfterCalc.txt");

        // Initialize SolModFactory from the GEMS3K file set
        SolModFactory task(input_system_file_list_name);
        //Check the data read for SolModFactory initialization
        auto after_reading = input_system_file_list_name;
        replace(after_reading, "-dat.lst", "-AfterReading.txt");
        task.to_text_file(after_reading );

        test_SolModFactory_API(task, template_file, input_system_file_list_name);

        test_SolModEngine_API(task, template_file, input_system_file_list_name);
        return 0;
    }
    catch(TError& err)
    {
        std::cout  << err.title << err.mess << std::endl;
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


// Test restore TMulty data after read the GEMS3K fileset
int test_SolModFactory_API(SolModFactory& task, const std::string& template_file_path, const std::string& ipmfiles_lst_name )
{
    // Load templalte data
    KeyValueFile templ_data(template_file_path);
    templ_data.load_all();

    std::string diff_string;
    // element
    diff_string += templ_data.compare_to( "Awt", task.Get_ElementMolarMasses());
    diff_string += templ_data.compare_to( "B", task.Get_ElementMoleAmounts());
    // phase
    diff_string += templ_data.compare_to( "YOF", task.Get_SurfaceFreeEnergyParameter());
    diff_string += templ_data.compare_to( "XF", task.Get_PhaseMoleAmounts());
    diff_string += templ_data.compare_to( "XFA", task.Get_PhaseCarrierMoleAmounts());
    // species
    diff_string += templ_data.compare_to( "DUL", task.Get_SpeciesUpperBounds());
    diff_string += templ_data.compare_to( "DLL", task.Get_SpeciesLowerBounds());
    diff_string += templ_data.compare_to( "MM", task.Get_SpeciesMolarMasses());
    diff_string += templ_data.compare_to( "Wx", task.Get_SpeciesMoleFractions());
    diff_string += templ_data.compare_to( "X", task.Get_SpeciesMoleAmounts());
    diff_string += templ_data.compare_to( "fDQF", task.Get_IncrementsMolarG0());
    diff_string += templ_data.compare_to( "Vol", task.Get_SpeciesMolarVolumes());
    diff_string += templ_data.compare_to( "G0", task.Get_GibbsEnergyG0());
    diff_string += templ_data.compare_to( "H0", task.Get_MolarEnthalpyH0());
    diff_string += templ_data.compare_to( "S0", task.Get_MolarEnropyS0());
    diff_string += templ_data.compare_to( "Cp0", task.Get_HeatCapacityCp0());
    diff_string += templ_data.compare_to( "A0", task.Get_HelmholtzEnergyA0());
    diff_string += templ_data.compare_to( "U0", task.Get_InternalEnergyU0());

    if(diff_string.empty() ) {
        std::cout <<  "Test SolModFactory API Pass: " << ipmfiles_lst_name <<  std::endl;
    }
    else {
        std::cout << "Test SolModFactory API Fail: " << ipmfiles_lst_name <<  std::endl;
        std::cout <<  "  Difference-------------------------------------------------";
        std::cout <<  diff_string << "\n\n";
    }
    return 0;
}

// Test restore TMulty data after read the GEMS3K fileset
int test_SolModEngine_API(SolModFactory& task, const std::string& template_file_path, const std::string& ipmfiles_lst_name )
{
    // Load templalte data
    KeyValueFile templ_data(template_file_path);
    templ_data.load_all();

    size_t jb=0, je=0;
    std::string diff_string_all;

    for(size_t k=0; k<task.Get_AllPhasesNumber(); ++k) {
      std::string diff_string;
      auto phase = task.Sol_Phase(k);
      jb=je;
      je+=phase.Get_SpeciesNumber();

      // Calculate activity coefficients of endmembers (components, species)
      phase.SolModActivityCoeffs();

      // Calculate (a dict) of ideal properties of mixing in the phase2
      //auto map_ideal = phase.SolModIdealProps();
      //std::cout << "Ideal properties of mixing in phase:  \n";
      //for(const auto& item: map_ideal ) {
      //    std::cout << "   '" << item.first << "': " << item.second << "; ";
      //}

      diff_string += templ_data.compare_subset( "Wx", phase.Get_MoleFractions(), jb, je);
      diff_string += templ_data.compare_subset( "Y_m", phase.Get_Molalities(), jb, je);
      //diff_string += templ_data.compare_subset( "lnActivities", phase.Get_lnActivities(), jb, je);
      diff_string += templ_data.compare_subset( "lnGam", phase.Get_lnActivityCoeffs(), jb, je);
      diff_string += templ_data.compare_subset( "lnCnft", phase.Get_lnConfTerms(), jb, je);
      diff_string += templ_data.compare_subset( "lnRcpt", phase.Get_lnRecipTerms(), jb, je);
      diff_string += templ_data.compare_subset( "lnExet", phase.Get_lnExcessTerms(), jb, je);
      diff_string += templ_data.compare_subset( "lnDQFt", phase.Get_lnDQFTerms(), jb, je);
      diff_string += templ_data.compare_subset( "fDQF", phase.Get_G0Increments(), jb, je);
      diff_string += templ_data.compare_subset( "Vol", phase.Get_MolarVolumes(), jb, je);
      diff_string += templ_data.compare_subset( "Pparc", phase.Get_PartialPressures(), jb, je);

      if( k<task.Get_SolPhasesNumber() ) {
          std::vector<double> VPh, GPh, HPh, SPh, CPh, APh, UPh;

          auto ex_map = phase.SolModStandProps();
          GPh.push_back(ex_map["Gst"]);
          HPh.push_back(ex_map["Hst"]);
          SPh.push_back(ex_map["Sst"]);
          CPh.push_back(ex_map["CPst"]);
          VPh.push_back(ex_map["Vst"]);
          APh.push_back(ex_map["Ast"]);
          UPh.push_back(ex_map["Ust"]);

          ex_map = phase.SolModIdealProps();
          GPh.push_back(ex_map["Gid"]);
          HPh.push_back(ex_map["Hid"]);
          SPh.push_back(ex_map["Sid"]);
          CPh.push_back(ex_map["CPid"]);
          VPh.push_back(ex_map["Vid"]);
          APh.push_back(ex_map["Aid"]);
          UPh.push_back(ex_map["Uid"]);

          ex_map = phase.SolModExcessProps();
          GPh.push_back(ex_map["Gex"]);
          HPh.push_back(ex_map["Hex"]);
          SPh.push_back(ex_map["Sex"]);
          CPh.push_back(ex_map["CPex"]);
          VPh.push_back(ex_map["Vex"]);
          APh.push_back(ex_map["Aex"]);
          UPh.push_back(ex_map["Uex"]);

          ex_map = phase.SolModDarkenProps();
          GPh.push_back(ex_map["Gdq"]);
          HPh.push_back(ex_map["Hdq"]);
          SPh.push_back(ex_map["Sdq"]);
          CPh.push_back(ex_map["CPdq"]);
          VPh.push_back(ex_map["Vdq"]);
          APh.push_back(ex_map["Adq"]);
          UPh.push_back(ex_map["Udq"]);

          diff_string += templ_data.compare_subset( "VPh", VPh, k*MIXPHPROPS, k*MIXPHPROPS+4);
          diff_string += templ_data.compare_subset( "GPh", GPh, k*MIXPHPROPS, k*MIXPHPROPS+4);
          diff_string += templ_data.compare_subset( "HPh", HPh, k*MIXPHPROPS, k*MIXPHPROPS+4);
          diff_string += templ_data.compare_subset( "SPh", SPh, k*MIXPHPROPS, k*MIXPHPROPS+4);
          diff_string += templ_data.compare_subset( "APh", APh, k*MIXPHPROPS, k*MIXPHPROPS+4);
          diff_string += templ_data.compare_subset( "UPh", UPh, k*MIXPHPROPS, k*MIXPHPROPS+4);
      }

      if(!diff_string.empty() ) {
          diff_string_all += "\n phase: " + phase.Get_SolPhaseName() + diff_string;
      }

    }

    if(diff_string_all.empty() ) {
        std::cout <<  "Test SolModEngine API Pass: " << ipmfiles_lst_name <<  std::endl;
    }
    else {
        std::cout << "Test SolModEngine API Fail: " << ipmfiles_lst_name <<  std::endl;
        std::cout <<  "  Difference-------------------------------------------------";
        std::cout <<  diff_string_all << "\n\n";
    }
    return 0;
}
