//-------------------------------------------------------------------
// $Id: io_arrays.h 774 2006-07-26 08:45:45Z gems $
//
// Service functions for writing/reading arrays in files
//
// Copyright (C) 2006-2007 S.Dmytriyeva
// Uses  gstring class
//
// This file is part of the GEM-Vizor library and GEMIPM2K
// code package
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------

#include  <fstream>

#ifdef IPMGEMPLUGIN

#include "verror.h"

#else

#include "gstring.h"

#endif

struct outField
 {
   char name[20]; // name of field in structure
   long int alws;    // 1 - must be read, 0 - default values can be used
   long int readed;  // 0; set to 1 after reading the field from input file
};

class TRWArrays  // basic class for red/write fields of structure
 {
 protected:
    fstream& ff;
    long int numFlds;
    outField* flds;

 public:

	 TRWArrays( short aNumFlds, outField* aFlds, fstream& fin ):
    	ff( fin ), numFlds(aNumFlds), flds(aFlds)
    {}

	virtual  long int findFld( const char *Name ); // find field by name

    void  setNoAlws( long int ii )
    {  flds[ii].alws = 0; }
    void  setNoAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
            setNoAlws(ii);
    }
    void  setAlws( long int ii )
    {  flds[ii].alws = 1; }
    void  setAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
            setAlws(ii);
    }

    bool  getAlws( long int ii )
    {  return (flds[ii].alws == 1); }
    bool  getAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
           return getAlws(ii);
         else
        return false;	 
    }

};


class TPrintArrays: public  TRWArrays  // print fields of structure
{
    inline void writeValue(float val);
    inline void writeValue(double val);

public:

    TPrintArrays( short aNumFlds, outField* aFlds, fstream& fout ):
    	TRWArrays( aNumFlds, aFlds, fout)
    {}

    void writeArray( const char *name, char*   arr, long int size, long int arr_size );
    void writeArray( const char *name, float*  arr, long int size, long int l_size=-1L );
    void writeArray( const char *name, double* arr, long int size, long int l_size=-1L );
    void writeArray( const char *name, long* arr, long int size, long int l_size=-1L  );

    void writeArray( const char *name, char*   arr, int size, int arr_size );
    void writeArray( const char *name, float*  arr, int size, int l_size=-1 );
    void writeArray( const char *name, double* arr, int size, int l_size=-1 );
    void writeArray( const char *name, short* arr, int size, int l_size=-1  );

    void writeArray( const char *name, float*  arr, long int size, long int* selAr,
    		long int nColumns=1L, long int l_size=-1L );
    void writeArray( const char *name, double* arr, long int size, long int* selAr,
    		long int nColumns=1L, long int l_size=-1L );
    void writeArray( const char *name, long* arr, long int size, long int* selAr,
    		long int nColumns=1L, long int l_size=-1L );

    void writeArray( const char *name, float*  arr, int size, long int* selAr,
    		int nColumns=1, int l_size=-1 );
    void writeArray( const char *name, double* arr, int size, long int* selAr,
    		int nColumns=1, int l_size=-1 );
    void writeArray( const char *name, short* arr, int size, long int* selAr,
    		int nColumns=1, int l_size=-1 );

};


 class TReadArrays : public  TRWArrays // read fields of structure
 {
    inline void readValue(float& val);
    inline void readValue(double& val);
 
 public:

    TReadArrays( short aNumFlds, outField* aFlds, fstream& fin ):
    	TRWArrays( aNumFlds, aFlds, fin )
    {}

    void  skipSpace();
    void reset();  // reset to 0 all flags (readed)

    long int findFld( const char *Name ); // find field by name
    long int findNext();  // read next name from file and find in fields list
    void  readNext( const char* label);

    gstring testRead();   // test for reading all arrays

    void readArray( const char *name, short* arr, long int size );
    void readArray( const char *name, int* arr, long int size );
    void readArray( const char *name, long int* arr, long int size );
    void readArray( const char *name, float* arr, long int size );
    void readArray( const char *name, double* arr, long int size );
    void readArray( const char *name, char* arr, long int size, long int el_size );

};

//=============================================================================
