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

void setpar_(npmax,xmin,xmax,partx,nbox)
double *xmin, *xmax;                    
double  partx[NCPMAX];
int *nbox;
int *npmax  ;

{
        register i ;
        double xco ;

        xco= (*xmax-*xmin) / (double) *npmax ;
        partx[0] = *xmin + xco/2.0;
        printf("setpar nbox  npmax %d %d xco  %f  %f \n",*nbox,*npmax,xco,partx[0]);
        for (i=1; i<=*npmax-1; i++)  {
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
#include "gwheader.h"

void partid_(npmax,nbox,xmin,xmax,partib,dx,partx)
double  partx[NCPMAX],dx[NCNODEX+2], *xmin, *xmax;
int *nbox,partib[NCNODEX];
int *npmax;

{
        register i, j;
        int iknx;

        for(i=0; i<= *nbox ; i++)   {

                partib[i]=0;
        }
        for (i=0; i<=*npmax-1; i++)  {

                if(partx[i] > *xmin && partx[i] < *xmax)   {

                        iknx= (int ) ( (partx[i]+dx[1])  / dx[2])  ;
                        partib[iknx] += 1 ;
                  }
                  else {
                      printf(" partibx %d %f  %d %d  \n",i,partx[i],iknx,partib[iknx]); 
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

void concver_(npmax,nbox,dx,bn,cn,partib,partx,partic,ismooth,m1,m2)
double  bn[NCNODEX][NCBASIS],cn[NCNODEX][NCCOMPL],partx[NCPMAX], partic[NCBASIS+NCCOMPL][NCPMAX],dx[NCNODEX+2];
int *nbox, partib[NCNODEX], *ismooth, *m1, *m2 ;
int *npmax;

{
        double dxinv, partiv[NCNODEX];
        register i,j,i1,i2,j2;
        int iknx;

        for(i2=0;i2<=*nbox; i2++)   {
           if(*ismooth ==1) {
              partiv[i2]=1. ;
           } 
           else {
              if(partib[i2]!=0.)partiv[i2]=1./(double) partib[i2];
           }
         }

        for (i=0; i<=*npmax-1; i++)  {
              iknx= (int )(( partx[i]+dx[1])  / dx[2] );

              for(j=0; j< *m1;j++) {
                   partic[j][i]= bn[iknx][j] * partiv[iknx];
              }        
              for(j=0; j< *m2;j++) { /*new 050494*/
                   partic[j+ *m1][i]=cn[iknx][j] *partiv[iknx];
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

void concneu_(npmax,nbox,nxmax,xminr,xmaxr,dx,bn,cn,partib,partx,partic,bo,
             co,ismooth,m1,m2)
double  bn[NCNODEX][NCBASIS],bo[NCNODEX][NCBASIS],cn[NCNODEX][NCCOMPL]
       ,co[NCNODEX][NCCOMPL];
double  partx[NCPMAX], partic[NCBASIS+NCCOMPL][NCPMAX];
double dx[NCNODEX+2];
double *xminr, *xmaxr;
int *nbox,*nxmax ,partib[NCNODEX], *ismooth,*m1, *m2 ;
int *npmax;

{
        register i,j,k,ky,i1,n;
        int iknx;
        double dxinv ;

/* set new concentrations at t +dt equal zero and sum over particles to get the new ones */
        for(i=0; i< *nxmax-1; i++) {
                partib[i]=0;
        }
        for(i=1; i< *nxmax-2; i++) {   /* cal_dol 1 durch 0 ersetzt */ 
/*                partib[i]=0; */
                for(j=0; j< *m1;j++) {
                        bn[i][j] = 0.0;
                }
                for(j=0; j< *m2;j++) { 
                       cn[i][j] =  0.0;
                }       
/*                      printf("boooooo i j %d %d %e \n",i,j,bo[i][j]); */
        }

        for (i=0; i<=*npmax-1; i++)  {

           if(partx[i] > *xminr && partx[i] < *xmaxr) {  /*+dx[2]*/
           
                  iknx= (int ) (( partx[i]+dx[1]) /dx[2] ) ;
                  if( iknx > 0  ) {
                        partib[iknx] += 1 ;
                        for(j=0; j< *m1;j++) {
                                bn[iknx][j] += partic[j][i];
                        }
                        for(j=0; j< *m2;j++) {  
                              cn[iknx][j] += partic[j+ *m1][i] ;
                        }
                   }
           }
         }
/*   printf("i  j bn partic partib  %d %d %e  %e %d  \n",i,j,bn[iknx][j],partic[j][i],partib[iknx]); */
                        
          

/* new part of particle concentration in a box  */


        if( *ismooth ==1) {
             for(n=0; n<= *nxmax-1 ;n++) {  /* cal_dol 1 durch 0 ersetzt */ 
                 for(j=0; j< *m1;j++) {
                    if(partib[n]!=0) bn[n][j]=bn[n][j]/(double)partib[n];
                 }
                 for(j=0; j< *m2;j++) {
                    if(partib[n]!=0) cn[n][j]=cn[n][j]/(double)partib[n];
                 }
              }
        }
/*  constant concentration at x=0   */
             for(j=0; j< *m1;j++) {
               bn[0][j] = bo[0][j] ;    /* inlet = const. concentration */
             }
             for(j=0; j< *m2;j++) {        
                cn[0][j] = co[0][j] ;    /* inlet = const. concentration */
             }
/*  const conc right side */


/*           for(j=0; j< *m1;j++) {
                 bn[*nxmax-1][j] = bo[*nxmax-1][j] ;
                }    
             for(j=0; j< *m2;j++) {
                 cn[*nxmax-1][j] = co[*nxmax-1][j] ;
                }
*/             for(j=0; j< *m1;j++) {
                   bn[*nxmax -2][j] = bn[*nxmax -3][j] ;
                   bn[*nxmax -1][j] = bn[*nxmax -2][j] ;
                }     
                for(j=0; j< *m2;j++) {
                     cn[*nxmax -2][j] = cn[*nxmax -3][j] ;
                     cn[*nxmax -1][j] = cn[*nxmax -2][j] ;
                }
 
}



