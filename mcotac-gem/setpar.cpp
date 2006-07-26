/*
*                            Unterprogramm: setpar
*                       Aussetzen von Partikeln
* Input:
*
*    npmax          :  Anzahl neu zu setzender Partikel
*    xmin, xmax     :  min und max -werte von x
*    ymin, ymax     :  min und max -werte von y
* Output:
*    partx[n] : Partikel coordinaten
*    party[n] : Partikel coordinaten
*
* 18.04.95
*==========================================================================*/

//#include "gwheader.h"
//#include <stdlib.h>
//#include <math.h>

extern "C" int __stdcall SETP(int npmax,double xmin,double xmax,double partx[50000],int nbox)

{


        register i ;
        double xco ;

        xco= (xmax- xmin) / (double)(npmax) ; /*                (double) */
        partx[0] = xmin + xco/2.0;
/*        printf("setpar nbox  npmax %d %d xco  %f  %f \n",*nbox,*npmax,xco,partx[0]);*/
        for (i=1; i<= npmax-1; i++)  {
                partx[i] =partx[i-1] + xco;
        }
   
return 0;
}

