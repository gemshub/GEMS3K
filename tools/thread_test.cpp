
#ifdef OVERFLOW_EXCEPT
#ifdef __linux__
#include <cfenv>
#elif _MSC_VER
#include <float.h>
#else
#include <cfenv>
#endif
#endif
#include <filesystem>
namespace fs = std::filesystem;

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include "node.h"
#include "jsonconfig.h"

std::vector<std::string> files_into_directory(const std::string& directory_path,
                                              const std::string& sample);

static void process_task(std::string path_to_lst)
{
    std::cout << "Test path:" << path_to_lst << std::endl;

    // Creates TNode structure instance accessible through the "node" pointer
    std::shared_ptr<TNode> node(new TNode());

    // (1) Initialization of GEMS3K internal data by reading  files
    //     whose names are given in the input_system_file_list_name
    if( node->GEM_init(path_to_lst.c_str()) ) {
        // error occured during reading the files
        std::cout << "error occured during reading the files: " << path_to_lst << std::endl;
        return;
    }
    // Getting direct access to work node DATABR structure which exchanges the
    // data with GEM IPM3 (already filled out by reading the DBR input file)
    DATABR* dBR = node->pCNode();
    dBR->NodeStatusCH = NEED_GEM_SIA;

    // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
    long NodeStatusCH = node->GEM_run( false );
    std::cout <<  path_to_lst << " files results: " << NodeStatusCH << std::endl;

}

// cmake -DCMAKE_BUILD_TYPE=Debug -fsanitize=thread ..
// Thermo-time-in/series1-dat.lst  Thermo-time-out/series1-dat.lst solvus-in/series1-dat.lst Kaolinite-in/pHtitr-dat.lst Neutral-fun/Neutral-dat.lst
int main(int argc, char* argv[])
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

    //Read config file
    //gemsSettings();
    gemsSettings().gems3k_update_loggers(false, "test.log", 3);

    try{
        std::vector<std::string> dat_lst_files;

        if( argc <= 2) {
            dat_lst_files = files_into_directory("gems3k", "-dat.lst");
        }
        else {
            size_t lst_count = 0;
            while( argv[++lst_count]) {
                dat_lst_files.push_back(argv[lst_count]);
            }
        }

        std::vector<std::thread> threads;
        //  Launch pool of worker threads
        for (const auto& file : dat_lst_files) {
            threads.push_back( std::thread(process_task, file) );
        }
        for (auto& th : threads)
            th.join();
        return 0;
    }
    catch( std::exception& e )
    {
        std::cout << "std::exception: " << e.what() <<  std::endl;
    }
    catch(...)
    {
        std::cout <<  "unknown exception" <<  std::endl;
    }

    return 1;
}

std::vector<std::string> files_into_directory(const std::string& directory_path,
                                              const std::string& sample)
{
    std::vector<std::string> file_names;

    if(!directory_path.empty()) {
        fs::path ps(directory_path);
        if(fs::exists(ps))  {
            for(auto& p: fs::directory_iterator(ps)) {
                if (fs::is_regular_file(p.path())) {
                    std::string file = p.path().string();
                    if ( sample.empty() || file.find(sample) != std::string::npos) {
                        //std::cout << "file =  " <<  file<<  std::endl;
                        file_names.push_back(file);
                    }
                }
                if (fs::is_directory(p.path())) {
                    std::string dir = p.path().string();
                    auto files = files_into_directory(dir, sample);
                    file_names.insert(file_names.end(), files.begin(), files.end());
                }
            }
        }
    }
    return file_names;
}
