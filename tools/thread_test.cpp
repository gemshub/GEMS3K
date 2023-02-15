
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
#include <string>
#include <vector>
#include <thread>
#include "node.h"
#include "jsonconfig.h"

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
    gemsSettings();

    try{
        size_t lst_count = 0;

        std::vector<std::thread> threads;
        //  Launch pool of worker threads
        while( argv[++lst_count]) {
            std::string lst_path = argv[lst_count];
            threads.push_back( std::thread(process_task, lst_path) );
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

