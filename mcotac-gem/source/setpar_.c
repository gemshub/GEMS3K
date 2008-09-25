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
void setpar_ ( int npmax,double xmin,double xmax,double partx[NUPMAX],int nbox )
#else
void setpar ( int npmax,double xmin,double xmax,double partx[NUPMAX],int nbox )
#endif
{
	register int i ;
	double xco /*, xx*/ ;


	xco= ( xmax-xmin ) / ( double ) ( npmax ) ; /*                (double) */
	partx[0] = xmin + xco/2.0;
	for ( i=1; i<npmax; i++ )
	{
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
void partid_ ( int npmax, int nbox, int, nxmax, double xmin, double xmax,
               int partib[NCNODEX+2], double dx[NCNODEX+2], double partx[NUPMAX] )

#else
void partid ( int npmax, int nbox, int nxmax, double xmin, double xmax,
              int partib[NCNODEX+2], double dx[NCNODEX+2], double partx[NUPMAX] )

#endif
{
	register int i;
	int iknx;
	/*        printf(" partid1   %d %d %f %f \n",*npmax,*nbox,*xmin,*xmax); */

	for ( i=0; i< nbox ; i++ )
	{
		/*              printf(" partid2   %d %d %d \n",i,partib[i], *nbox); */
		partib[i]=0;

	}
	for ( i=0; i<npmax; i++ )
	{

		iknx= ( int ) ( ( partx[i] )  / dx[1] )  ;
		if ( iknx < nbox && iknx >= 0 )
		{

			partib[iknx] += 1 ;
			/*   printf(" partid5   %d %d %f %f \n",i,iknx,partx[i],dx[2]); */
		}
		else
		{
			/*  printf(" partid6   %d  %f \n",i,dx[2]); */
			printf ( " partibx out of range particle %d iknx %i    \n",i,iknx );
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
void concver_ ( int npmax,int nbox,double dx[NCNODEX+2],double bn[NCNODEX+2][NCBASIS],
                double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],double partx[NUPMAX],
                double partic[NCCOMPL+NCBASIS][NUPMAX],int ismooth,int m1,int m2 )
#else
void concver ( int npmax,int nbox,double dx[NCNODEX+2],double bn[NCNODEX+2][NCBASIS],
               double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],double partx[NUPMAX],
               double partic[NCCOMPL+NCBASIS][NUPMAX],int ismooth,int m1,int m2 )
#endif
{
	double partiv[NCNODEX+2];
	register int i,j,i2;
	int iknx;

	for ( i2=0;i2<nbox; i2++ )
	{
		if ( ismooth ==1 )
		{
			partiv[i2]=1. ;
		}
		else
		{
			partiv[i2]=1. ;
			if ( partib[i2]!=0. ) partiv[i2]=1./ ( double ) partib[i2];
		}
	}

	for ( i=0; i<npmax; i++ )
	{
		iknx= ( int ) ( ( partx[i] )  / dx[2] );

		for ( j=0; j< m1;j++ )
		{
			partic[j][i]= bn[iknx][j] * partiv[iknx];
		}
		for ( j=0; j< m2;j++ )   /*new 050494*/
		{
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
void concneu_ ( int npmax,int nbox,int nxmax,
                double xmin,double xmaxr,double dx[NCNODEX+2],
                double  bn[NCNODEX+2][NCBASIS],double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],
                double partx[NUPMAX], double partxo[NUPMAX],double partic[NCCOMPL+NCBASIS][NUPMAX],
                double bo[NCNODEX+2][NCBASIS], double co[NCNODEX+2][NCCOMPL],
                double por[NCNODEX+2], int ismooth, int m1,int m2 )
#else
void concneu ( int npmax,int nbox,int nxmax,
               double xmin,double xmaxr,double dx[NCNODEX+2],
               double  bn[NCNODEX+2][NCBASIS],double cn[NCNODEX+2][NCCOMPL],int partib[NCNODEX+2],
               double partx[NUPMAX], double partxo[NUPMAX],double partic[NCCOMPL+NCBASIS][NUPMAX],
               double bo[NCNODEX+2][NCBASIS], double co[NCNODEX+2][NCCOMPL], double por[NCNODEX+2],
               int ismooth, int m1,int m2 )
#endif
{
	int i,j,n;

	int iknx,iknxo,total_part,partout;
	double dxinv , x_partib[NCNODEX+2],total_x_part,Pn,Po;

	/* set new concentrations at t +dt equal zero and sum over particles to get the new ones */
	total_part=0.0;
	total_x_part=0.0;
	partout=0;
	for ( i=0; i< nbox; i++ )
	{
		partib[i]=0;
		x_partib[i]=0.;
	}
	for ( i=0; i< nbox; i++ )
	{
		for ( j=0; j< m1;j++ )
		{
			bn[i][j] = 0.0;
		}
		for ( j=0; j< m2;j++ )
		{
			cn[i][j] =  0.0;
		}
	}

	for ( i=0; i<npmax; i++ )
	{
		/*        for (j=0; j<=*m1; j++) {
		        printf("i jjbn partic partib  %d %d %e  %e %d  \n",i,j,bn[iknx][j],partic[j][i],partib[iknx]);
		        } */
		/*               printf("concneu i partx [i]  %d %f %f %f \n",i, *xminr, *xmaxr, partx [i]);*/
		/*               printf("i partx  %d %e \n",i,partx[i]);*/
		iknx= ( int ) ( ( partx[i] ) /dx[1] ) ;
		iknxo= ( int ) ( ( partxo[i] ) /dx[1] ) ;
		partib[iknx] +=1; /* add one for particle position */
		/*               printf("i iknx  %d %d \n",i,iknx); */
		/*			if ( (iknx !=0) && (iknx != nbox-1)  && (iknxo !=0) && (iknxo != nbox-1) ){ */

		if ( por[iknxo]<=por[iknx] )
		{
			Pn=1.0;Po=0.0;
		}
		else
		{
			Pn=por[iknx] / por[iknxo] ; Po=1-Pn;
		}
		/*					if ((iknx<2) || (iknx>nbox-2) || (iknxo<2) || (iknxo>nbox-2)) {Pn=1.0;Po=0.0;}		*/
		/* Pn=por[iknx] / ( por[iknx] + por[iknxo] ); */   /*probability to enter new domain */
		/* Po=1.0-Pn;  */     /* probability for old domain */

		x_partib[iknx]+= Pn;    /* mass is distributed according to probability */
		x_partib[iknxo]+= Po;

		for ( j=0; j< m1;j++ )
		{
			bn[iknx][j] +=partic[j][i] * Pn;
			bn[iknxo][j]+=partic[j][i] * Po;
		}
		for ( j=0; j< m2;j++ )
		{
			cn[iknx][j] +=partic[j+ m1][i] * Pn;
			cn[iknxo][j]+=partic[j+ m1][i] * Po;
		}


	}



	/* new part of particle concentration in a box  */
	for ( n=0; n< nbox ;n++ )
	{
		/* partib[n]= ( int ) ( x_partib[n] ); */ /* what is this? number of particles in box/at node...then this is wrong for variable porosity?*/
		total_part+=partib[n];
		total_x_part+=x_partib[n];
	}
	printf ( " total_part: %d partout %d total_x_partx: %d \n",total_part+partout,partout, total_x_part );

	if ( ismooth ==1 )
	{
		for ( n=0; n< nbox;n++ )
		{
			for ( j=0; j< m1;j++ )
			{
				if ( partib[n]!=0 ) bn[n][j]=bn[n][j]/  partib[n];
			}
			for ( j=0; j< m2;j++ )
			{
				if ( partib[n]!=0 ) cn[n][j]=cn[n][j]/ partib[n];
			}
		}
	}



	/*  constant concentration at x=0   */

	for ( j=0; j< m1;j++ )
	{
		bn[0][j] = bo[0][j] ;
	}
	for ( j=0; j< m2;j++ )
	{
		cn[0][j] = co[0][j] ;
	}
	/*	             for(j=0; j< m1;j++) {
		               bn[1][j] = bo[0][j] ;
		            }
		             for(j=0; j< m2;j++) {
		                cn[1][j] = co[0][j] ;
		             } */

	/* constant concentration also at right boundary */
	for ( j=0; j< m1;j++ )
	{
		bn[nxmax-1][j] = bo[nxmax-1][j] ;
	}
	for ( j=0; j< m2;j++ )
	{
		cn[nxmax-1][j] = co[nxmax-1][j] ;
	}
	/*	             for(j=0; j< m1;j++) {
		               bn[nxmax-2][j] = bo[nxmax-1][j] ;
		            }
		             for(j=0; j< m2;j++) {
		                cn[nxmax-2][j] = co[nxmax-1][j] ;
		             } */


}



