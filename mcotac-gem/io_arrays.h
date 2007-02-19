
//#include <iomanip>
#include  <fstream>
#include "gstring.h"

class TPrintArrays  // print fields of structure
{
    fstream& ff;

public:

    TPrintArrays( fstream& fout ):
      ff( fout ){}

    void writeArray( char *name, char*   arr, int size, int arr_siz );
    void writeArray( char *name, short*  arr, int size, int l_size=-1  );
    void writeArray( char *name, float*  arr, int size, int l_size=-1 );
    void writeArray( char *name, double* arr, int size, int l_size=-1 );

};


struct outField
 {
   char name[20]; // name of field in structure
   short alws;    // 1 - must be read, 0 - default values can be used
   short readed;  // 0; set to 1 after reading the field from input file
};

 class TReadArrays  // read fields of structure
 {

    fstream& ff;
    short numFlds;
    outField* flds;

 public:

    TReadArrays( short aNumFlds, outField* aFlds, fstream& fin ):
      numFlds(aNumFlds), flds(aFlds), ff( fin )
    {}

    void  skipSpace();
    void reset();  // reset to 0 all flags (readed)

    short findFld( const char *Name ); // find field by name
    short findNext();  // read next name from file and find in fields list
    void  setNoAlws( short ii )
    {  flds[ii].alws = 0; }

    void  setNoAlws( const char *Name )
    {
         short ii = findFld( Name );
         if( ii >=0 )
            setNoAlws(ii);
    }

    gstring testRead();   // test for reading all always arrays;

    void readArray( char *name, short* arr, int size );
    void readArray( char *name, float* arr, int size );
    void readArray( char *name, double* arr, int size );
    void readArray( char *name, char* arr, int size, int el_size );

};

//=============================================================================
