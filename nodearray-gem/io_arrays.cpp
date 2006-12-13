#include <iomanip>
#include  <iostream>

#include "io_arrays.h"
#include "verror.h"

#ifdef IPMGEMPLUGIN
  istream& f_getline(istream& is, gstring& str, char delim);
#endif


//---------------------------------------------------------//
// print Arrays ( fields of structure )

void TPrintArrays::writeArray( char *name, char* arr,
                              int size, int arr_siz )
{
 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == 10)
    { jj=0;  ff << endl;}
    gstring str = gstring( arr +(ii*arr_siz), 0, arr_siz );
    str.strip();
    ff  << "\'" << str.c_str() << "\'" << " ";
 }
}

void TPrintArrays::writeArray( char *name, short* arr,
                 int size, int l_size  )
{
  int sz = 10;
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

void TPrintArrays::writeArray( char *name,  float* arr,
            int size, int l_size )
{
 int sz = 10;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << setprecision(10) << scientific << arr[ii] << " ";
 }
}

void TPrintArrays::writeArray( char *name,  double* arr,
            int size, int l_size )
{
 int sz = 10;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << setprecision(18) << scientific << arr[ii] << " ";
 }
}

//-------------------------------------------------------------------------

//------------------------------------------------------------------

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
 for(short ii=0; ii < numFlds; ii++ )
    flds[ii].readed = 0;
}

short TReadArrays::findFld( const char *Name )
{
 short ii;
 gstring str = Name;
 size_t len = str.find('>');
 str = str.substr(0, len );

 for( ii=0; ii < numFlds; ii++ )
  if( !( strcmp( flds[ii].name, str.c_str() ) ))
    return ii;
 return -1;
}

short TReadArrays::findNext()
{
 char buf[200];
 skipSpace();

 if( ff.eof() )
   return -3;

 ff >> buf;

 if( !( memcmp( "END_DIM", buf+1, 7 )) )
  return -2;

 short ii = findFld( buf+1 );
 if(  ii < 0 )
    Error( buf, "DataBR text read 01: Invalid name of array");

 flds[ii].readed = 1;
 return ii;
}


void TReadArrays::readArray( char*, short* arr, int size )
{
 for( int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff >> arr[ii];
 }
}

void TReadArrays::readArray( char*, float* arr, int size )
{
 for( int ii=0; ii<size; ii++  )
 {
     skipSpace();
     ff >> arr[ii];
 }
}

void TReadArrays::readArray( char*, double* arr, int size )
{
 for( int ii=0; ii<size; ii++  )
 {
     skipSpace();
     ff >> arr[ii];
 }
}

void TReadArrays::readArray( char*, char* arr, int size, int el_size )
{
 char ch;
 char buf[200];

 for( int ii=0; ii<size; ii++  )
 {
   skipSpace();
   ff.get(ch);
   while( ff.good() && ch != '\'' )
       ff.get(ch);
   ff.getline( buf, el_size+1, '\'');
   strncpy( arr +(ii*el_size), buf, el_size );
 }

}

gstring TReadArrays::testRead()
{
 gstring ret = "";
 for(short ii=0; ii < numFlds; ii++ )
  if( flds[ii].alws==1 && flds[ii].readed != 1 )
  {  if( !ret.empty() )
       ret += ", ";
     ret += flds[ii].name;
  }
 return ret;
}

//=============================================================================
