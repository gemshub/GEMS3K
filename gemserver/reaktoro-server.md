## ZeroMQ server { GEMS-Reaktoro project } ##

Now we are at the end of Stage 1 having the first working variant GEMS3K-server of GEMS-Reaktoro code.
With GEMS-Reaktoro GUI, we can set up a chemical equilibrium problem as usual and compute Systems and Processes using the GEMS3K code as a standalone server.

The accuracy of results and their consistency with the trunk GEM-Selektor version 3.7.0 has been tested on Kaolinite and Solvus projects; 
more detailed comparisons will be done in the next tasks.

Main focus at Stage 2 will be to modify Reaktoro as a standalone server using ZeroMQ messaging in a similar way as it has been already done with GEMS3K, and  to be able to run Reaktoro server alternatively to GEMS3K server and compare the results of solving the same test equilibria and processes.

This task may be facilitated by using the new JSON format of GEMS3K I/O messages, already readable/writable by GEMS3K server and possible to be used by the Reaktoro/Optima server.



### Implementation ###

Implementation of ZeroMQ GEMS3K-server you can see: https://bitbucket.org/gems4/gems3k/src/GEMS-Reactoro/gemserver/



* main.cpp

Implementation of ZeroMQ exchange. Would be the same for GEMS3K-server and Reaktoro-server


* tnodeinterface.h

NodeInterface interface, that must be implemented for GEMS3K-server and Reaktoro-server

```c++

/// Initialization of internal data from json/key-value strings
/// Parameters:
///  @param dch_json -  DATACH - the Data for CHemistry data structure as a json string
///  @param ipm_json -  Multi structure as a json string
///  @param dbr_json -  DATABR - the data bridge structure as a json string
///  @return array of error type, error title and error message in case of error or empty array if all OK

```



```c++

/// Run process of calculate equilibria into the GEMS3K/Reaktoro side
/// Parameters:
///  @param dbr_json -  DATABR - the data bridge structure as a json string
///  @return array with strings contains:
///    - NodeStatus codes with respect to GEMIPM calculations
///    - the result DATABR - the data bridge structure as a json string
///    - calculation time in seconds elapsed during the last call of GEM_run (obsolete)
///    - number of IPM loops performed ( >1 up to 6 because of PSSC() ) (obsolete)
///    - Number of completed IA EFD iterations (obsolete)
///    - Number of completed GEM IPM iterations (obsolete)
///    or if NodeStatus is error
///    - error title
///    - error message
virtual std::vector<std::string> calculateEquilibrium( const std::string& new_dbr ) = 0;

```


```c++

/// \enum NODECODECH NodeStatus codes with respect to GEMIPM calculations
/*typedef*/ enum NODECODECH {
 NO_GEM_SOLVER= 0,   ///< No GEM re-calculation needed for this node
 NEED_GEM_AIA = 1,   ///< Need GEM calculation with LPP (automatic) initial approximation (AIA)
 OK_GEM_AIA   = 2,   ///< OK after GEM calculation with LPP AIA
 BAD_GEM_AIA  = 3,   ///< Bad (not fully trustful) result after GEM calculation with LPP AIA
 ERR_GEM_AIA  = 4,   ///< Failure (no result) in GEM calculation with LPP AIA
 NEED_GEM_SIA = 5,   ///< Need GEM calculation with no-LPP (smart) IA, SIA
                     ///<   using the previous speciation (full DATABR lists only)
 OK_GEM_SIA   = 6,   ///< OK after GEM calculation with SIA
 BAD_GEM_SIA  = 7,   ///< Bad (not fully trustful) result after GEM calculation with SIA
 ERR_GEM_SIA  = 8,   ///< Failure (no result) in GEM calculation with SIA
 T_ERROR_GEM  = 9    ///< Terminal error has occurred in GEMS3K (e.g. memory corruption). Restart is required.
} /*NODECODECH*/;

```


* tnodetask.h(.cpp)

Implementation NodeInterface for GEMS3K-server








> To generate test files in json format for server you can use command "Data/Export GEMS3K files ..."

> You can use key-value format data exchange if comment ``` DEFINES  += JSON_OUT``` in client and server projects.


### Some addition applications

* jsonioDiff

  Comparing json to json; key-value to key-value or key-value to json files. The order of data fields in the document is arbitrary (but in array fields the order of values is fixed). Comparing floating point numbers choosing epsilon depends on the context, and determines how equal you want the numbers to be.

* jsonioDiff/examples/kv_to_json

  Converter key-value format file to json format file.
