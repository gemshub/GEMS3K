/*==========================================================================*
*
*                            Unterprogramm: walk2
*                       Berechnung eines Zeitschritts
*
* Input:
*
*    nxmax        : Anzahl der Knoten in X-Richtung
*    icyc      : Zeitzyklus-Nummer
*    along, aq    : longitudinale bzw. transversale Dispersivit„t
*    de        : Zeitschrittweite =delt
*    vx[i]  : Geschwindigkeitskomponente in X-Richtung
*    dx[i]     : Knotenabst„nde in X-Richtung
*
*    ir[i][j]: ir-Maske - Transportmodell, enth„lt Randbedingungen
*
* Output:
*
*    partx[i]       : Koordinaten des Partikels i
*    npmax        : Anzahl Partikel im Modell
*
* 16.11.92
*==========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "gwheader.h"

#define SQRT12 3.464101615
#define randinv 1.0/RAND_MAX
float gasdev();

void walk2_(npmax,nxmax,ncyc,along,aquer,dm,texe,dx,vx,partx,partxo,
     xmaxr,xminr,partic,bn,cn,partib,ibpstart,x,bo,co,m1,m2)


double partx[NCPMAX],partxo[NCPMAX],dx[NCNODEX+2],vx[NCNODEX+2],*xmaxr,*xminr, *texe;
double  partic[NCBASIS+NCCOMPL][NCPMAX];
double x[NCNODEX],bn[NCNODEX][NCBASIS],cn[NCNODEX][NCCOMPL] ,bo[NCNODEX][NCBASIS],co[NCNODEX][NCCOMPL] ;
double dm[NCNODEX+2];
double *along,*aquer;
int    *nxmax,  *ncyc ,partib[NCNODEX],*ibpstart,*m1, *m2 ;
int   *npmax;

{
        double  slong,slongy,xlaenge, xxmin, xxmax;
        double dabs(),dpx, A1,A2,A3,A4 ,vpx,Z1,vabs;
        register i, ip,ipa;
/*      float gasdef(idum);*/
        int idum, iknx;
        char dummy;

        idum=1;
        ip = 0;
        ipa =0;
/*      printf("walk c %d %d %d %g %g% %g %g\n",*ibpstart,*nxmax,*ncyc,*along,*texe,*xmaxr, *xminr);
      printf("walk %d %d %d %g %g %g %g %g \n",*ibpstart,*nymax,*ncyc,*along,*dm,*texe,*ymaxr, *yminr);

        for (i=0; i<=(*nxmax)-1; i++) {
                printf("i dx vx %d %f %f \n", i, dx[i], vx[i] );
        }
*/
        xlaenge = *xmaxr - *xminr + 2 * dx[1];
        xxmin = *xminr - dx[1];
        xxmax = *xmaxr + dx[1];     /*  + dx[1]; */

        do  {  /*  teilchenloop  */
           iknx= (int ) ( (partx[ip]+dx[1])  / dx[2])  ;

           dpx=dx[1]-(partx[ip]-x[iknx]);
           vpx=vx[iknx];
           vabs=sqrt(vpx*vpx);
           if(vabs == 0.)  vabs=1.e-30; /* keine division durch 0 */
/* neu 23 05 95 */

           if (*along == 0. && *aquer == 0. && dm[iknx] == 0.)  {              /* disp + diff =0  */
             slong=0.0;
           }
           else {                                                     /* wenigstens eine Dispersivitaet > 0 */

              Z1=(double) rand()*randinv -0.5 ;
              slong = 2.*Z1* sqrt(6.* (*along*vpx  + dm[iknx]) * *texe);  /* longitu. Weg x' */

           }

           partx[ip] +=  vpx *  *texe +  slong;  /* neue position der teilchen = konv. Anteil + disp. Anteil x */

           if(partx[ip]  >=   xxmax  ){  /*keine diffusion ueber 'rechten rand' - const. conc   */
                         /*   printf("walkrand ip dx[1] partxo xmaxr xlaenge %d %f %f %f %f\n", ip, dx[1] ,partxo[ip], *xmaxr,xlaenge);*/
                         /* neue Teilchenposition als ob Teilchen "links" im Randgrid neu eingesetztwird, Konzentration dabei egal, da nicht gerechnet */
/*             if(partx[ip]- partxo[ip] > 2*dx[1])printf(">xmax ip  %d %f %f \n", ip,partxo[ip], partx[ip] );*/
             partx[ip] -=  xlaenge ;
             for (ipa=0; ipa< *m1; ipa++) {
                partic[ipa][ip]= bn[0][ipa]/(double)*ibpstart ;/* division durch 50 partikelx in randbox */

             }
             for (ipa=0; ipa< *m2; ipa++) {
               partic[ipa+ *m1][ip]= cn[0][ipa]/(double)*ibpstart;
             }
           }
           if(partx[ip] <= xxmin) {
/*             if(-1.*(partx[ip]- partxo[ip]) > dx[1])printf("< xmin ip  %d %f %f \n", ip,partxo[ip], partx[ip] );*/
             partxo[ip] =  partx[ip] ;
               partx[ip] += xlaenge   ;
               for (ipa=0; ipa< *m1; ipa++) {
                  partic[ipa][ip]= bn[*nxmax-2][ipa]/(double)*ibpstart; /* -1 durch -2 ersetzt, 090796 */
               }
               for (ipa=0; ipa< *m2; ipa++) {
                  partic[ipa+ *m1][ip]=  cn[*nxmax-2][ipa]/(double)*ibpstart; /*   -1 durch -2 ersetzt, 090796 */
               }
            }
 /*           printf("walk11 %d %d  %g %g\n",*ibpstart,ip,partx[ip], partxo[ip]);*/
        ++ip;
/*      scanf(" walk %s \n" , &dummy) ;
      printf("walk11 %d %d  %s\n",*ibpstart,ip,dummy);
*/        }   while (ip<= *npmax-1  );      /* ende teilchenloop */
}




#include <math.h>
#include <stdlib.h>
#define randinv 1.0/RAND_MAX

float gasdev(idum)
int *idum;
/* returns a normally distributed deviate with  zero mean and unit variance */

{
        int iset=0;
        float gset;
        float fac,r,v1,v2;

        if (iset ==0) {
                do {
                        v1=2.0*rand()*randinv - 1.0;
                        v2=2.0*rand()*randinv - 1.0;
                        r=v1*v1+v2*v2;
                } while( r >= 1.0 || r ==0 );
                fac = sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=1;
                return gset;
        }
}
