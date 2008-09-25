/*==========================================================================*
*
*                            Unterprogramm: walk2
*                       Berechnung eines Zeitschritts
*
* Input:
*
*    nxmax        : Anzahl der Knoten in X-Richtung
*    icyc      : Zeitzyklus-Nummer
*    along, aq    : longitudinale bzw. transversale Dispersivit�t
*    de        : Zeitschrittweite =delt
*    vx[i]  : Geschwindigkeitskomponente in X-Richtung
*    dx[i]     : Knotenabst�nde in X-Richtung
*
*    ir[i][j]: ir-Maske - Transportmodell, enth�lt Randbedingungen
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


double gasdev ( int *idum )
/*int *idum;*/
/* returns a normally distributed deviate with  zero mean and unit variance */

{
	int iset=0;
	double gset;
	double fac,r,v1,v2;
	gset=0.;
	if ( iset ==0 )
	{
		do
		{
			v1=2.0*rand ( /*idum*/ ) *randinv - 1.0;
			v2=2.0*rand ( /*idum*/ ) *randinv - 1.0;
			r=v1*v1+v2*v2;
		}
		while ( r >= 1.0 || r ==0 );
		fac = sqrt ( -2.0*log ( r ) /r );
		gset=v1*fac;
		iset=1;
		return v2*fac;
	}
	else
	{
		iset=1;
		return gset;
	}
}

/*==========================================================================*
*
*                            Unterprogramm: walk2
*                       Berechnung eines Zeitschritts
*
* Input:
*
*    nxmax        : Anzahl der Knoten in X-Richtung
*    icyc      : Zeitzyklus-Nummer
*    along, aq    : longitudinale bzw. transversale Dispersivit�t
*    de        : Zeitschrittweite =delt
*    vx[i]  : Geschwindigkeitskomponente in X-Richtung
*    dx[i]     : Knotenabst�nde in X-Richtung
*
*    ir[i][j]: ir-Maske - Transportmodell, enth�lt Randbedingungen
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
double gasdev();
#ifdef __PGI
void walk2h_ ( int npmax,int nxmax,int ncyc,double along,double aquer,double dm[NCNODEX+2]
               ,double texe,double dx[NCNODEX+2],double vx[NCNODEX+2]
               ,double partx[NCPMAX],double partxo[NCPMAX],double xmax,double xmin
               ,double partic[NCBASIS+NCCOMPL][NCPMAX],double bn[NCNODEX+2][NCBASIS]
               ,double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],int ibpstart,double x[NCNODEX+2]
               ,double bo[NCNODEX+2][NCBASIS],double co[NCNODEX+2][NCCOMPL],int m1,int m2, double por[NCNODEX+2],int ismooth )
#else
void walk2h ( int npmax,int nxmax,int ncyc,double along,double aquer,double dm[NCNODEX+2]
              ,double texe,double dx[NCNODEX+2],double vx[NCNODEX+2]
              ,double partx[NCPMAX],double partxo[NCPMAX],double xmax,double xmin
              ,double partic[NCBASIS+NCCOMPL][NCPMAX],double bn[NCNODEX+2][NCBASIS]
              ,double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],int ibpstart,double x[NCNODEX+2]
              ,double bo[NCNODEX+2][NCBASIS],double co[NCNODEX+2][NCCOMPL],int m1,int m2, double por[NCNODEX+2],int ismooth )
#endif

{
	double  slong,slongy,xlaenge, xxmin, xxmax;
	double dabs(),dpx, A1,A2,A3,A4 ,vpx,Z1,vabs,partx_dt;
	register i, ip,ipa;
	/*      float gasdef(idum);*/
	int idum, iknx,iknxx;
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
	/*        xlaenge = xmaxr - xminr + 2 * dx[1]; */
	xlaenge=xmax-xmin;  /* does only work for regular grid */
	xxmin = xmin;
	xxmax = xmax;     /*  + dx[1]; */

	do     /*  teilchenloop  */
	{
		iknx= ( int ) ( ( partx[ip] )  / dx[2] )  ;

		partxo[ip]=partx[ip];
		dpx= ( partx[ip]-x[iknx] );
		vpx=vx[iknx];
		vabs=sqrt ( vpx*vpx );
		if ( vabs == 0. )  vabs=1.e-30; /* keine division durch 0 */
		/* neu 23 05 95 */

		if ( along == 0. && aquer == 0. && dm[iknx] == 0. )               /* disp + diff =0  */
		{
			slong=0.0;
		}
		else                                                       /* wenigstens eine Dispersivitaet > 0 */
		{

			Z1= ( double ) rand() *randinv -0.5 ;
			slong = 2.*Z1* sqrt ( 6.* ( along*vpx  + dm[iknx] ) * texe );  /* longitu. Weg x' */

		}
		/*           partx[ip] +=  vpx *  *texe +  slong; */ /* neue position der teilchen = konv. Anteil + disp. Anteil x */
		partx_dt =  vpx *  texe +  slong;  /* neue position der teilchen = konv. Anteil + disp. Anteil x */
		partx[ip] +=partx_dt ;
		/* reflection of particles */

		/*
		           iknxx= (int ) ( (partx[ip]+partx_dt+dx[1])  / dx[2])  ;
		           if(iknx == iknxx || (por[iknx]==por[iknxx])){
			     partx[ip] +=partx_dt ;
			     }
			     else {

		                   if(por[iknxx] > por[iknx]) {
		                      partx[ip]+=partx_dt;
		                      }
				       Z1=((double) rand()*randinv) ;
		                       if (Z1 <= (sqrt(por[iknxx])/(sqrt(por[iknx])+sqrt(por[iknxx]))) ){
				        partx[ip]+=partx_dt;
					}
					else {

					partx[ip]= partx[ip]-(partx_dt-2*(dx[1]-dpx));
					}



			     }
		*/
		/*end reflection */
		iknx= ( int ) ( ( partx[ip] )  / dx[2] )  ;
		if ( partx[ip]>xmax ) /* - const. conc   */
		{
			partx[ip] = xmax- ( partx[ip]-xmax )  ;  /* changed to reflection boundary */
			iknx= ( int ) ( ( partx[ip] )  / dx[2] )  ;
			for ( ipa=0; ipa< m1; ipa++ )
			{
				if ( ismooth ==1 )
				{
					partic[ipa][ip]= bn[iknx][ipa] ;/* no division */
				}
				else
				{
					partic[ipa][ip]= bn[iknx][ipa]/partib[iknx] ;/* division durch partikel in randbox ---siehe setpar_*/
				}

			}
			for ( ipa=0; ipa< m2; ipa++ )
			{
				if ( ismooth ==1 )
				{
					partic[ipa+ m1][ip]= cn[iknx][ipa] ;/* no division */
				}
				else
				{
					partic[ipa+ m1][ip]= cn[iknx][ipa]/partib[iknx] ;/* division durch partikel in randbox ---siehe setpar_*/
				}
			}
		}

		if ( partx[ip]<xmin )
		{
			/*             if(-1.*(partx[ip]- partxo[ip]) > dx[1])printf("< xmin ip  %d %f %f \n", ip,partxo[ip], partx[ip] );*/
			partx[ip] = xmin+ ( xmin-partx[ip] )   ;   /* changed to reflection boundary */
			iknx= ( int ) ( ( partx[ip] )  / dx[2] )  ;
			for ( ipa=0; ipa< m1; ipa++ )
			{
				if ( ismooth ==1 )
				{
					partic[ipa][ip]= bn[iknx][ipa] ;  /* no division */
				}
				else
				{
					partic[ipa][ip]= bn[iknx][ipa]/partib[iknx] ;   /* division durch partikel in randbox ---siehe setpar_*/
				}
			}
			for ( ipa=0; ipa< m2; ipa++ )
			{
				if ( ismooth ==1 )
				{
					partic[ipa+ m1][ip]=  cn[iknx][ipa] ;/* no division */
				}
				else
				{
					partic[ipa+ m1][ip]=  cn[iknx][ipa]/partib[iknx] ;  /* division durch partikel in randbox ---siehe setpar_*/
				}
			}
		}
		/*           printf("walk11 %d %d  %g %g\n",*ibpstart,ip,partx[ip], partxo[ip]);*/
		++ip;
		/*      scanf(" walk %s \n" , &dummy) ;
		      printf("walk11 %d %d  %s\n",*ibpstart,ip,dummy);
		*/
	}
	while ( ip<= npmax-1 );      /* ende teilchenloop */
}