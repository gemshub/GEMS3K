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

#include <stdlib.h>
#include <math.h>
#include "gwheader.h"

#ifdef __PGI
   void setpar_(int npmax,double xmin,double xmax,double partx[NUPMAX],int nbox)
#else
   void setpar(int npmax,double xmin,double xmax,double partx[NUPMAX],int nbox)
#endif
{
        register int i ;
        double xco /*, xx*/ ;


        xco= (xmax-xmin) / (double)(npmax) ; /*                (double) */
        partx[0] = xmin + xco/2.0;
        for (i=1; i<=npmax-1; i++)  {
                partx[i] =partx[i-1] + xco;
        }
}


/*************************************************************************
*                            Unterprogramm: partid
*                          Partikeln in welcher box 
* Input:
*
*    npmax          :  Anzahl neu zu setzender Partikel
*    nbox            :  Anzahl box zwischen xmin und xmax 
*    partx[n]       :     Partikel coordinaten
* Output:
*     partb(nb)   : Anz Partikel in box nb
* 28.10.92
*==========================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "gwheader.h"

#ifdef __PGI
   void partid_( int npmax, int nbox, double xmin, double xmax,
	                        int partib[NCNODEX+2], double dx[NCNODEX+2], double partx[NUPMAX])

#else
   void partid( int npmax, int nbox, double xmin, double xmax,
	                        int partib[NCNODEX+2], double dx[NCNODEX+2], double partx[NUPMAX])

#endif
{
        register int i;
        int iknx;
/*        printf(" partid1   %d %d %f %f \n",*npmax,*nbox,*xmin,*xmax); */

        for(i=0; i<= nbox ; i++)   {
  /*              printf(" partid2   %d %d %d \n",i,partib[i], *nbox); */
                partib[i]=0;

        }
        for (i=0; i<=npmax-1; i++)  {

                if(partx[i] > xmin && partx[i] < xmax)   {

                        iknx= (int ) ( (partx[i]+dx[2])  / dx[2])  ;
                        partib[iknx] += 1 ;
                    /*   printf(" partid5   %d %d %f %f \n",i,iknx,partx[i],dx[2]); */
                   }
                  else {
                     /*  printf(" partid6   %d  %f \n",i,dx[2]); */
                      printf(" partibx out of range %d %le    \n",i,partx[i]); 
                }
        }
}

/*************************************************************************
*                            Unterprogramm: concver
*                          Partikelconcentration festlegen 
* Input:
*
*    npmax          :  Anzahl neu zu setzender Partikel
*    nbox            :  Anzahl box zwischen xmin und xmax 
*    partx[n]       :     Partikel coordinaten
*    partib[ib]    :   
* Output:
*    partic[n,nspec]: partikel n mit speciesconcentrationen nspec
* 28.10.92
*==========================================================================*/

#include <stdlib.h>
#include <math.h>
#include "gwheader.h"

#ifdef __PGI
 void concver_(int npmax,int nbox,double dx[NCNODEX],double bn[NCNODEX+2][NCBASIS],
			 double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX],double partx[NUPMAX],
			 double partic[NCCOMPL+NCBASIS][NUPMAX],int ismooth,int m1,int m2)
#else
 void concver(int npmax,int nbox,double dx[NCNODEX],double bn[NCNODEX+2][NCBASIS],
			 double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX],double partx[NUPMAX],
			 double partic[NCCOMPL+NCBASIS][NUPMAX],int ismooth,int m1,int m2)
#endif
{
        double partiv[NCNODEX+2];
        register int i,j,i2;
        int iknx;

        for(i2=0;i2<=nbox; i2++)   {
           if(ismooth ==1) {
              partiv[i2]=1. ;
           } 
           else {
              if(partib[i2]!=0.)partiv[i2]=1./(double) partib[i2];
           }
         }

        for (i=0; i<=npmax-1; i++)  {
              iknx= (int )(( partx[i]+dx[1])  / dx[2] );

              for(j=0; j< m1;j++) {
                   partic[j][i]= bn[iknx][j] * partiv[iknx];
              }        
              for(j=0; j< m2;j++) { /*new 050494*/
                   partic[j+ m1][i]=cn[iknx][j] * partiv[iknx];
              }
        }

}



/*************************************************************************
*                            Unterprogramm: concneu
*                          Partikelconcentration festlegen fuer t+dt
* Input:
*
*    npmax          :  Anzahl neu zu setzender Partikel
*    nbox            :  Anzahl box zwischen xmin und xmax 
*    partx[n]       :     Partikel coordinaten
*    partib[ib]    :   
* Output:
*    partic[n,nspec]: partikel n mit speciesconcentrationen nspec
* 28.10.92
*==========================================================================*/

#include <stdlib.h>
#include <math.h>
#include "gwheader.h"

#ifdef __PGI
  void concneu_(int npmax,int nbox,int nxmax,
             double xminr,double xmaxr,double dx[NCNODEX+2],
			 double  bn[NCNODEX+2][NCCOMPL],double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],
			 double partx[NUPMAX],double partic[NCCOMPL+NCBASIS][NUPMAX],
			 double bo[NCNODEX+2][NCCOMPL], double co[NCNODEX+2][NCCOMPL], 
			 int ismooth, int m1,int m2)
#else
  void concneu(int npmax,int nbox,int nxmax,
             double xminr,double xmaxr,double dx[NCNODEX+2],
			 double  bn[NCNODEX+2][NCCOMPL],double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],
			 double partx[NUPMAX],double partic[NCCOMPL+NCBASIS][NUPMAX],
			 double bo[NCNODEX+2][NCCOMPL], double co[NCNODEX+2][NCCOMPL], 
			 int ismooth, int m1,int m2)
#endif
  {
        int i,j,n;
        int iknx;


/* set new concentrations at t +dt equal zero and sum over particles to get the new ones */
        for(i=0; i< nxmax-1; i++) {
                partib[i]=0;
        }
        for(i=1; i< nxmax-2; i++) {   /* cal_dol 1 durch 0 ersetzt */ 
/*                partib[i]=0; */
                for(j=0; j< m1;j++) {
                        bn[i][j] = 0.0;
                }
                for(j=0; j< m2;j++) { 
                       cn[i][j] =  0.0;
                }       
/*                      printf("boooooo i j %d %d %e \n",i,j,bo[i][j]); */
        }

        for (i=0; i<=npmax-1; i++)  {
/*        for (j=0; j<=*m1; j++) {
        printf("i jjbn partic partib  %d %d %e  %e %d  \n",i,j,bn[iknx][j],partic[j][i],partib[iknx]); 
        } */
/*               printf("concneu i partx [i]  %d %f %f %f \n",i, *xminr, *xmaxr, partx [i]);*/
           if(partx[i] > xminr  && partx[i] < xmaxr) {  /**/
/*               printf("i partx  %d %e \n",i,partx[i]);*/
                  iknx= (int ) (( partx[i]+dx[2]) /dx[2] ) ;
/*               printf("i iknx  %d %d \n",i,iknx);*/
                   
                  if( iknx > 0  ) {
                        partib[iknx] += 1 ;
                        for(j=0; j< m1;j++) {
                                bn[iknx][j] += partic[j][i];
                        }
                        for(j=0; j< m2;j++) {  
                              cn[iknx][j] += partic[j+ m1][i] ;
                        }
                   }
           }
         }
                        
          

/* new part of particle concentration in a box  */


        if( ismooth ==1) {
             for(n=0; n<= nxmax-1 ;n++) {  /* cal_dol 1 durch 0 ersetzt */ 
                 for(j=0; j< m1;j++) {
                    if(partib[n]!=0) bn[n][j]=bn[n][j]/(double)partib[n];
                 }
                 for(j=0; j< m2;j++) {
                    if(partib[n]!=0) cn[n][j]=cn[n][j]/(double)partib[n];
                 }
              }
        }
/*  constant concentration at x=0   */
             for(j=0; j< m1;j++) {
               bn[0][j] = bo[0][j] ;    /* inlet = const. concentration */
            }
             for(j=0; j< m2;j++) {        
                cn[0][j] = co[0][j] ;    /* inlet = const. concentration */
             }
/*	transmission boundary  'rigth side' 
		bn[*nx -2][j] = -bn[*nx-4][j]+2.*bn[*nx -3][j] ;
		cn[*nx -2][j] = -cn[*nx-4][j]+2.*cn[*nx -3][j] ;
		bn[*nx -1][j] = -bn[*nx-3][j]+2.*bn[*nx -2][j] ;
		cn[*nx -1][j] = -cn[*nx-3][j]+2.*cn[*nx -2][j] ;
	                if(bn[*nx-1][j] <0. ) bn[*nx-1][j]=0.;
	                if(cn[*nx-1][j] <0. ) cn[*nx-1][j]=0.;
*/
    /*          for(j=0; j< *m1;j++) {
                }    
                for(j=0; j< *m2;j++) {
                }
               for(j=0; j< *m1;j++) {
                }     
                for(j=0; j< *m2;j++) {
                }
*/
}



