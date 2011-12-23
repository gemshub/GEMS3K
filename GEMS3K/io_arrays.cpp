//-------------------------------------------------------------------
// $Id: io_arrays.cpp 774 2006-07-26 08:45:45Z gems $
//
// Implementation of service functions for writing/reading arrays in files
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

#include <iomanip>
#include  <iostream>

#include "io_arrays.h"
#include "verror.h"
#include "v_user.h"


istream& f_getline(istream& is, gstring& str, char delim);

//    This constants should be 'defined' to satisfy all compilers
#define SHORT_EMPTY 	   -32768
#define LONG_EMPTY             -2147483648L
#define FLOAT_EMPTY	          1.17549435e-38F
#define DOUBLE_EMPTY         2.2250738585072014e-308
#define CHAR_EMPTY   	     '`'

inline bool IsFloatEmpty( const float v )
{
    return ( v>0. && v <= FLOAT_EMPTY);
}
inline bool IsDoubleEmpty( const double v )
{
    return ( v>0. && v <= DOUBLE_EMPTY);
}


long int TRWArrays::findFld( const char *Name )
  {
  	long int ii;

   for( ii=0; ii < numFlds; ii++ )
    if( !( strcmp( flds[ii].name, Name ) ))
      return ii;
   return -1;
  }

//---------------------------------------------------------//
// print Arrays ( fields of structure )
inline void TPrintArrays::writeValue(float val)
    {
      if( IsFloatEmpty( val ))
        ff << CHAR_EMPTY << " ";
      else
    //    ff << setprecision(10) << scientific << arr[ii] << " ";
       ff << setprecision(7) << val << " ";
    }

inline void TPrintArrays::writeValue(double val)
    {
      if( IsDoubleEmpty( val ))
        ff << CHAR_EMPTY << " ";
      else
//    ff << setprecision(18) << scientific << arr[ii] << " ";
    ff << setprecision(15) << val << " ";
    }

// If the first parameter is given as NULL then the char array
// will be printed as a comment
void TPrintArrays::writeArray( const char *name, char* arr,
		long int size, long int arr_siz )
{
 bool isComment = false;

 if( name )
     ff << endl << "<" << name << ">" << endl;
 else
 { ff << endl << "#  ";
   isComment = true;
 }
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == 40 )
    { jj=0;  ff << endl;
      if(isComment)
    	  ff << "#  ";
    }
    gstring str = gstring( arr +(ii*arr_siz), 0, arr_siz );
    strip(str);
    ff  << "\'" << str.c_str() << "\'" << " ";
 }
}

void TPrintArrays::writeArray( const char *name, char* arr,
		int size, int arr_siz )
{
 bool isComment = false;

 if( name )
     ff << endl << "<" << name << ">" << endl;
 else
 { ff << endl << "#  ";
   isComment = true;
 }
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == 40 )
    { jj=0;  ff << endl;
      if(isComment)
    	  ff << "#  ";
    }
    gstring str = gstring( arr +(ii*arr_siz), 0, arr_siz );
    strip(str);
    ff  << "\'" << str.c_str() << "\'" << " ";
 }
}

void TPrintArrays::writeArray( const char *name,  float* arr,
		long int size, long int l_size )
{
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name,  float* arr,
		int size, int l_size )
{
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		long int size, long int l_size )
{
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		int size, int l_size )
{
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    writeValue(arr[ii]);
 }
}

void TPrintArrays::writeArray( const char *name, long* arr, long int size, long int l_size  )
{
 long int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << arr[ii] << " ";
 }
}


void TPrintArrays::writeArray( const char *name, short* arr, int size, int l_size  )
{
 int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << arr[ii] << " ";
 }
}

//-------------------------------------------------------------------------

void TPrintArrays::writeArray( const char *name,  float* arr,
		long int size, long int* selArr, long int nColumns, long int l_size )
{
 if(!arr) return;
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++  )
 {
    for(long int cc=0; cc<nColumns; cc++ )
    {
    	if(jj == sz)
	    { jj=0;  ff << endl;}
         writeValue(arr[selArr[ii]*nColumns+cc]);
        jj++;
   }
 }
}

void TPrintArrays::writeArray( const char *name,  float* arr,
		int size, long int* selArr, int nColumns, int l_size )
{
 if(!arr) return;
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++  )
 {
   for( int cc=0; cc<nColumns; cc++ )
    {
    	if(jj == sz)
	    { jj=0;  ff << endl;}
       writeValue(arr[selArr[ii]*nColumns+cc]);
        jj++;
    }
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		long int size, long int* selArr, long int nColumns, long int l_size )
{
 if(!arr) return;
 long int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++  )
 {
    for(long int cc=0; cc<nColumns; cc++ )
        {
              if(jj == sz)
	        { jj=0;  ff << endl;}
              writeValue(arr[selArr[ii]*nColumns+cc]);
              jj++;
        }
 }
}

void TPrintArrays::writeArray( const char *name,  double* arr,
		int size, long int* selArr, int nColumns, int l_size )
{
 if(!arr) return;
 int sz = 40;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++  )
 {
     for( int cc=0; cc<nColumns; cc++ )
        {
           if(jj == sz)
           { jj=0;  ff << endl;}
             writeValue(arr[selArr[ii]*nColumns+cc]);
             jj++;
           }
 }
}


void TPrintArrays::writeArray( const char *name, long* arr,
		long int size, long int* selArr, long int nColumns, long int l_size )
{
  if(!arr) return;
  long int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( long int ii=0, jj=0; ii<size; ii++  )
 {
	for(long int cc=0; cc<nColumns; cc++ )
    {
		if(jj == sz)
        { jj=0;  ff << endl;}
    	ff << arr[selArr[ii]*nColumns+cc] << " ";
    	jj++;
    }
 }
}

void TPrintArrays::writeArray( const char *name, short* arr,
		int size, long int* selArr, int nColumns, int l_size )
{
  if(!arr) return;
  int sz = 40;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++  )
 {
	for( int cc=0; cc<nColumns; cc++ )
    {
		if(jj == sz)
        { jj=0;  ff << endl;}
    	ff << arr[selArr[ii]*nColumns+cc] << " ";
    	jj++;
    }
 }
}

//-------------------------------------------------------------------------
//------------------------------------------------------------------

 inline void TReadArrays::readValue(float& val)
 {
   char input;
   skipSpace();

   ff.get( input );
   if( input == CHAR_EMPTY )
      val = FLOAT_EMPTY;
   else
      {
        ff.putback(input);
        ff >> val;
      }
  }

 inline void TReadArrays::readValue(double& val)
 {
   char input;
   skipSpace();

   ff.get( input );
   if( input == CHAR_EMPTY )
      val = DOUBLE_EMPTY;
   else
      {
        ff.putback(input);
        ff >> val;
      }
 }

// skip  ' ',  '\n', '\t' and comments (from '#' to end of line)
void  TReadArrays::skipSpace()
{
  char input;
  ff.get( input );
  while( input == '#' || input == ' ' ||
        input == '\n' || input == '\t')
 {
   if( input == '#' )
    do{
         ff.get( input );
      }while( input != '\n' && input != '\0' && !ff.eof());
   if( input == '\0' || ff.eof() )
     return;
   ff.get( input );
  }
 ff.putback(input);
}

void TReadArrays::reset()
{
 for(long int ii=0; ii < numFlds; ii++ )
    flds[ii].readed = 0;
}

long int TReadArrays::findFld( const char *Name )
{
	long int ii;
 gstring str = Name;
 size_t len = str.find('>');
 str = str.substr(0, len );

 for( ii=0; ii < numFlds; ii++ )
  if( !( strcmp( flds[ii].name, str.c_str() ) ))
    return ii;
 return -1;
}

long int TReadArrays::findNext()
{
 char buf[200];
 skipSpace();

 if( ff.eof() )
   return -3;

 ff >> buf;

 if( !( memcmp( "END_DIM", buf+1, 7 )) )
  return -2;

 long int ii = findFld( buf+1 );
 if(  ii < 0 )
    Error( buf, "Format text read 01: Invalid label of data");

 flds[ii].readed = 1;
 return ii;
}


void TReadArrays::readNext( const char* label)
{
 char buf[200];
 skipSpace();

 if( ff.eof() )
   Error( label, "Format text read 02: No data where expected");

 ff >> buf;
 gstring str = buf+1;
 size_t len = str.find('>');
 str = str.substr(0, len );

 if( !( strcmp( label, str.c_str() ) ))
     return;
 Error( buf, "Format text read 03: Invalid label of data");
}

void TReadArrays::readArray( const char*, short* arr, long int size )
{
 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff >> arr[ii];
 }
}

void TReadArrays::readArray( const char*, int* arr, long int size )
{
 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff >> arr[ii];
 }
}

void TReadArrays::readArray( const char*, long int* arr, long int size )
{
 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff >> arr[ii];
 }
}
void TReadArrays::readArray( const char*, float* arr, long int size )
{
 for( long int ii=0; ii<size; ii++  )
    readValue(arr[ii]);
}

void TReadArrays::readArray( const char*, double* arr, long int size )
{
 //ff << setprecision(15);
 for( long int ii=0; ii<size; ii++  )
    readValue(arr[ii]);
}

void TReadArrays::readArray( const char*, char* arr, long int size, long int el_size )
{
 char ch;
 char buf[200];

 for( long int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff.get(ch);
//   while( ff.good() && ch != '\'' )
//       ff.get(ch);
   ff.getline( buf, el_size+1, '\'');
   copyValues( arr +(ii*el_size), buf, el_size );
 }

}

gstring TReadArrays::testRead()
{
 gstring ret = "";
 for(long int ii=0; ii < numFlds; ii++ )
  if( flds[ii].alws==1 && flds[ii].readed != 1 )
  {  if( !ret.empty() )
       ret += ", ";
     ret += flds[ii].name;
  }
 return ret;
}

//=============================================================================
// io_arrays.cpp
