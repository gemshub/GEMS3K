/****************************************************************************
*                                                                           *
*                            routine: hydro1d_.c                                 *
*          Stroemungsberechnung fuer Finite-Differenzen-Modell              *
* 21.11.02                                                                         *
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gwheader.h"

double irmaske(int ir[NCNODEX+2], int nxmax,
   double  tx[NCNODEX+2], double st[NCNODEX+2], double h0[NCNODEX+2]);
double vrech(double hb[NCNODEX+2], double tx[NCNODEX+2],
   double am[NCNODEX+2], double por[NCNODEX+2],
   double dx[NCNODEX+2], int nxmax, double vx[NCNODEX+2]);
double bilanz(char datei6[10], int nxmax,
             double dx[NCNODEX+2], int  ir[NCNODEX+2],
             double  hb[NCNODEX+2], double st[NCNODEX+2],
             double am[NCNODEX+2], double por[NCNODEX+2],
             double qw[NCNODEX+2], double vx[NCNODEX+2],
             double qbil[NCNODEX+2]);
double gsv(double hb[NCNODEX+2], double h0[NCNODEX+2],
           double qw[NCNODEX+2], double dx[NCNODEX+2],
           double tx[NCNODEX+2], double st[NCNODEX+2],
           int nxmax, int iter, float de,float re);
double info( int nxmax, int iter, double time,
             double fehler, double texe, int ij);
int wegdat1d_( int nxmax, char fname[10],
    double hb[NCNODEX+2], char text[10]);

#ifdef __unix
  extern "C"  double hydro1d_(int& nxmax,double h0[NCNODEX+2],double hb[NCNODEX+2],double tx[NCNODEX+2],
			  double am[NCNODEX+2],double st[NCNODEX+2],double por[NCNODEX+2],int	ir[NCNODEX+2]
			  ,double qw[NCNODEX+2],double qbil[NCNODEX+2],char& text,double vx[NCNODEX+2]
			  ,double dx[NCNODEX+2],int& icyc,double& texe,double& time,
              char& fname)
#else
 double hydro1d(int nxmax,double h0[NCNODEX+2],double hb[NCNODEX+2],double tx[NCNODEX+2],
			  double am[NCNODEX+2],double st[NCNODEX+2],double por[NCNODEX+2],int	ir[NCNODEX+2]
			  ,double qw[NCNODEX+2],double qbil[NCNODEX+2],char text,double vx[NCNODEX+2]
			  ,double dx[NCNODEX+2],int icyc,double texe,double time,
              char fname)
#endif
  {
	double anfkt=1.0, fehler=1. , re=1.6 ,er=1.e-6;

        double    de,xx=0.,angle=0.,
		xkomin=0.;

	int ierr = 0, iweg = 1,itermax=5000000;
	int iinst=0, ifrei=0,  ntim, ij,
		iter, itest=0, isolve=0 ,  itrans, ihb;
	int i;
	char   datei6[10];

/*        FILE *output;*/

/*	printf("hydro called %d\n", *icyc);*/

	strcpy(datei6,"bilanz.bil");
/*	strcat(datei1, ".ism");*/

	itrans = 1;
	ihb = 0;
	ij=0;
		/* } */

	ntim = 1;
	de   = 1.;

 	irmaske(ir,nxmax,tx,st,h0);


              for(i=2; i<=(nxmax)+1; i++)  {
		hb[i] = h0[i];
 /* printf("hb ", hb[i]  ) ;  */

		vx[i]   = 0;
		qbil[i] = 0;
	}

	dx[(nxmax) + 1] = dx[nxmax];


	wegdat1d_(nxmax ,"tx.dat" , tx,  "Hb");
 	/* Berechnung der Transmissivitaeten zwischen den Knoten */
/*nov2002	if( *icyc == 0) {
	for(i=1; i<=(*nxmax); i++)
		for(j=1; j<=(*nymax); j++)  {
         if( (tx[j][i] + tx[j+1][i]) != 0.)
	 	      ty[j][i] = (dy[j] + dy[j + 1]) * tx[j][i] * tx[j+1][i] /
		      (dy[j + 1] * tx[j][i] + dy[j] * tx[j+1][i]);
	          if( (ty[j][i] + ty[j+1][i]) != 0. && ty[j+1][i]!=0.)
	 	      ty[j][i] = (dy[j] + dy[j + 1]) * ty[j][i] * ty[j+1][i] /
		      (dy[j + 1] * ty[j][i] + dy[j] * ty[j+1][i]);
		  if( (tx[j][i] + tx[j][i+1]) != 0. && tx[j+1][i]!=0.)
		      tx[j][i] = anfkt * (dx[i] + dx[i + 1]) * tx[j][i] * tx[j][i+1] /
		      (dx[i + 1] * tx[j][i] + dx[i] * tx[j][i+1]);

		 tx[j][i] = 2 * tx[j][i] * dx[i] / (dy[j] + dy[j + 1]);
		 ty[j][i] = 2 * ty[j][i] * dy[j] / (dx[i] + dx[i + 1]);


        }
	}
nov 2002*/
	wegdat1d_(nxmax,"tx.dat.dat" , tx,  "Hb");

	/*  Start der Simulation  */


	for(i=1; i<=(nxmax); i++)   /*  0 ODER 1  ? */
	 {
/*	     printf("iterat	i j h0[j] %d %d %g ir= %d \n",i,j,h0[i], ir[i]); */
               h0[i] = hb[i];
	       qw[i] =  qw[i];   /*qr[i];*/
	}
	iter = 1;    /*  iterationsschleife */
	wegdat1d_(nxmax,"porb.dat" , por,  "Hb");
	if (icyc == 0) {
	     wegdat1d_(nxmax,"stbe.dat" , st,  "Hb");
	     wegdat1d_(nxmax,"ambe.dat" , am,  "Hb");
	     wegdat1d_(nxmax,"Hbe0.dat" , hb,  "Hb");
	}
	do  {

	  fehler = gsv(hb,h0,qw,dx,tx,st,nxmax,iter,de,re);

	  info(nxmax,iter,time,fehler,texe,ij);

	  for(i=1; i<=(nxmax); i++)
 	      h0 [i] = hb [i];

	  }  while( fehler>=er  && iter++<itermax);

	wegdat1d_(nxmax,"Hber.dat" , hb,  "Hb");

	vrech(hb,tx, am,por,dx,nxmax,vx);


	bilanz(datei6,nxmax,dx,ir,hb,st,am,por,qw,vx,qbil);

 	wegdat1d_(nxmax ,"qbil.dat",qbil,"Qbil" ) ;
	wegdat1d_(nxmax ,"vxx.dat" ,  vx,  "Vx" ) ;

	wegdat1d_(nxmax ,"qwbe.dat" , qw,  "Hb");
	wegdat1d_(nxmax ,"qbib.dat" , qbil,  "Hb");
/*printf("nach bilanz \n");*/


	if (iter > itermax) {
       	   wegdat1d_(nxmax,"H_err.dat" , hb,  "Hb");
           exit(1);
        }
return (iter);
}


/*==========================================================================*
*
*                            Unterprogramm: irmaske
*                       Einarbeiten der Randbedingungen
*
*  Input:
*
*     ir[i][j] : ir-Maske
*                Kennziffer 0: Knoten au�erhalb
*                           1: Knoten innerhalb
*                           2: Randknoten mit vorgegebenem Zuflu� Qr
*                              (Qr kann auch = 0 sein)
*                           4: Knoten mit fester Standrohrspiegelhoehe
*	                    5: Leakage-Knoten (-)
*                           6: Entnahme-Knoten (-)
*
*     nxmax   : Anzahl der Knoten in X- bzw. in Y-Richtung
*
*  Output: (entspr. der ir-Maske geaenderte Wertefelder)
*
*    tt[i] : Lokale Transmissivitaeten
*    st[i]  : Speicherkoeffizienten
*    h0[i] : Standrohrspiegelhoehenn
*
*  15.08.95
*===========================================================================*/
#include "gwheader.h"

double irmaske(int ir[NCNODEX+2], int nxmax,
   double  tx[NCNODEX+2], double st[NCNODEX+2], double h0[NCNODEX+2])
{

	int i;


          for(i=1; i<= (nxmax); i++)
		st[i] = 0.;       /*stationaer*/


           for(i=1; i<= (nxmax); i++)
               switch(ir[i])   {
		case 0:
/*	        printf("i j  ir %d %d %d %e\n", i,j,ir[i],tx[i]);	*/
			tx[i] = 0.;
			st[i] = 0.;
   		break;
		case 2:  /*  Nachbarknoten muss Ir==0 haben ..*/
/*jun2001*/
/*
                if( (ir[j+1][i]==0 && (ir[j+1][i-1] ==0 || ir[j+1][i+1] ==0)) ||
                    (ir[j-1][i]==0 && (ir[j-1][i-1] ==0 || ir[j-1][i+1] ==0 )) )
			ty[j][i] = 0.;
                if( (ir[j][i+1]==0 && (ir[j-1][i+1] ==0 || ir[j+1][i+1] ==0)) ||
                    (ir[j][i-1]==0 && (ir[j-1][i-1] ==0 || ir[j+1][i-1] ==0 )) )
			tx[j][i] = 0.;
                printf("case 2 %d %d %d %g %g \n",i,j,ir[j][i],tx[j][i],ty[j][i]);
*/
		break;
		case 4:
			st[i] = 1.e+20;
		break;
		case 6:  /*  Entnahmeknoten  */

		break;
	/*	default:
		break; */
 		}
	return(nxmax);
}
/*==========================================================================*
*
*                            Unterprogramm: GSV
* Loesung des Gleichungssystems nach dem Gauss-Seidel-Verfahren mit Relaxation
*
*
*   H[i]  : Standrohrspiegelhoehen
*   H0[i] : Anfangsstandrohrspiegelhoehen
*   Q[i]  : Knoten-Zufluesse
*   Dx[i]    : Zellgroessen in X-Richtung
*   Tx[i] : Lokale Transmissivitaeten in X-Richtung
*   Nx  : Anzahl der Knoten in X- bzw. in Y-Richtung
*   Iter     : Nummer des laufenden Iterationsschrittes
*   Re       : Relaxationsfaktor
*
*   return =  max. Fehler der Standrohrspiegelhoehen E
*
*===========================================================================*/
#include <stdlib.h>
#include <math.h>
#include "gwheader.h"

double gsv(double hb[NCNODEX+2], double h0[NCNODEX+2],
           double qw[NCNODEX+2], double dx[NCNODEX+2],
           double tx[NCNODEX+2], double st[NCNODEX+2],
           int nxmax, int iter, float de,float re)
/*
double hb[NCNODEX+2], h0[NCNODEX+2], qw[NCNODEX+2],
       tx[NCNODEX+2], st[NCNODEX+2],
       dx[NCNODEX+2];

float  re, de;
int (nxmax),  iter;
*/
{
	int i;
	double m1, m2, m3, m4, xne, za, hi;
	double fehler = 0.;

	m1 = m2 = m3 = m4 = 1.;	/* gespannter Aquifer */


	for (i=1; i<=(nxmax+1); i++)

		if (st[i] / dx[i]  <= 1.)  {
/*		printf("gsvgsv i hb tx %d %e %e \n",i,hb[i],tx[i]); */
			xne  = tx[i] * m2 + tx[i- 1 ] * m4
				+ st[i] / 1. /* de */ ;
			if (xne != 0.)  {
				za  = hb[i- 1] * tx[i - 1] * m4
					+ hb[i + 1] * tx[i] * m2

					 + qw[i]
					+ h0[i] * st[i] / 1. /* de */ ;
				hi = za / xne;
				fehler = max(fehler, fabs(hb[i] - hi));
				hb[i] += re * (hi - hb[i]);
			}
		}	/* if st[i] */
/*	} 	 for ii, for jj */
	return(fehler);
}

/*==========================================================================
*
*                       Unterprogramm: Info
*
*    Iter         : Anzahl Iterationen
*    E            : Max. Knotenfehler in [m]
*    TSim         : Simulierte Zeit
*    ICyc         : Zeitzyklus-Nummer
*    IJ           : Zeitschritt-Nummer
*
*==========================================================================*/
#include <stdio.h>
#include <math.h>

double info( int nxmax, int iter, double time,
             double fehler, double texe, int ij)
/*
double fehler, time, texe;
int (nxmax),  iter,   ij;
*/
{
	int iswitch;
	iswitch = ( (nxmax)>3000 ? 1 : 50);
	if (fmod((double)iter,(double)iswitch) == 0.)  {
		printf("anzahl iterationen    :  %d\n", iter);
		printf("max fehler pro knoten :  %e\n", fehler);

		printf("zeit            :  %f\n",time);
		printf("zeitintervall   :  %f\n\n", texe);
/*		printf("zeitschritt nr. :  %d\n\n\n", ij);*/

	}
	return(iter);
 }

/*==========================================================================
*
*                       Unterprogramm: Vrech
*         Berechnung des V-Feldes (Abstandsgeschwindigkeit)
*
* Input:
*
*   hb[i]  : Standrohrspiegelhoehen
*   tx[i] : Lokale Transmissivitaeten (falls Ifrei=1: Kf-Werte)
*   por[i]: Effektive Porositaet
*   am[i] : Aquifermaechtigkeit
*   dx[i]    : Zellgroessen in X-Richtung
*   Nx   : Anzahl der Knoten in X- bzw. in Y-Richtung
*
* Output:
*
*   vx[i] : Geschwindigkeitskomponente in X-Richtung
*
*   21.11.02
*==========================================================================*/
#include "gwheader.h"

double vrech(double hb[NCNODEX+2], double tx[NCNODEX+2],
             double am[NCNODEX+2], double por[NCNODEX+2],
             double dx[NCNODEX+2], int nxmax, double vx[NCNODEX+2])
/*
double hb[NCNODEX+2], tx[NCNODEX+2],
       am[NCNODEX+2], por[NCNODEX+2], vx[NCNODEX+2],
       dx[NCNODEX+2] ;

int (nxmax);
*/
{
	int i;
/*	float ma;*/

	/* gespannter aquifer: tx ist transmissivitaet */

	  for (i=1; i<=(nxmax-3); i++)  {

	    vx[i] = tx[i]  * (hb[i] - hb[i + 1]) / dx[i] / por[i];
	  }

           vx[0] = vx[1] *por[1]/ por[0];
           vx[nxmax-2] = vx[nxmax-3] *por[nxmax-3]/ por[nxmax-2];
           vx[nxmax-1] = vx[nxmax-2] *por[nxmax-2]/ por[nxmax-1];
	return(nxmax);
}          /*   ^ was *nxmax-1    Dec 2002  */


/*==========================================================================
*
*                       Unterprogramm: Bilanz
*                 Wasserbilanz fuer Finite-Differenzen-Modell
*
*   datei    : Name der Bilanz-Datei
*   nx  : Anzahl der Knoten in X-Richtung
*   icyc     : Zeitzyklus-Nummer
*   dx[i]    : Zellgroessen in X-Richtung
*   ir[i] : Ir-Maske (Kennziffernfeld fuer Randbedingungen)
*   hb[i]  : Standrohrspiegelhoehen
*   st[i]  : Speicherkoeffizienten
*   am[i] : Aquifermaechtigkeit
*   por[i]: Effektive Porositaet
*   vx[i] : (Abstands-) Geschwindigkeitskomponente in X-Richtung
*
* Output:
*
*   qbil[i] : Knotenbilanz
*
* Auf Datenfile: Einzelbilanzen, Gesamtbilanz, Bilanzmatrix
*
*  21.11.02
*==========================================================================*/
#include <stdio.h>
#include "gwheader.h"

double bilanz(char datei6[10], int nxmax,
              double dx[NCNODEX+2], int  ir[NCNODEX+2],
              double  hb[NCNODEX+2], double st[NCNODEX+2],
              double am[NCNODEX+2], double por[NCNODEX+2],
              double qw[NCNODEX+2], double vx[NCNODEX+2],
              double qbil[NCNODEX+2])
/*
char datei6[10];
double hb[NCNODEX+2], st[NCNODEX+2], am[NCNODEX+2],
       vx[NCNODEX+2],
       por[NCNODEX+2], qbil[NCNODEX+2], qw[NCNODEX+2],
       dx[NCNODEX+2];
int (nxmax),   ir[NCNODEX+2];
*/
{
	FILE *output;
	double ma, q1, q2,  qg, qknoten,  qzu, qab;
	int i;

	double qfehlzu = 0;                 /* fehler */
	double qfehlab = 0;
	double qfixzu = 0;                 /* fest vorgegebener Randzuflu� */
	double qfixab = 0;
	double qneubzu = 0;                /* q-grundwasserneubildung */
	double qfestzu = 0;                /* q-festpotential gesamtgebiet */
	double qfestab = 0;

	/* Berechnung der Bilanzen */

	for (i=1; i<=(nxmax); i++)
            {
	    ma = 1;     /*area vor */

	    q1 = -vx[i] * por[i] * ma;

	    q2 = vx[i - 1] * por[i - 1]  * ma;
	    qknoten = -q1 - q2 ;
	    qbil[i] = qknoten;
	    switch (ir[i])  {
		case 1:
/*			qneubzu += qneu[i];*/
		break;
		case 4:
			if (qknoten<0)
				qfestab += qknoten;
			else
				qfestzu += qknoten;
		break;
		case 2:
			if (qknoten<0)
				qfixab += 0.+   qw[i];  /* qknoten - qneu[i]*wf +   qr[i];*/
			else
				qfixzu += 0.+   qw[i]; /* qknoten - qneu[i]*wf +   qr[i];*/
		break;
		case 6:                /* qfixab and qfix zu including wells */
			if (qknoten<0)
				qfixab += 0.+   qw[i];  /* qknoten - qneu[i]*wf +   qr[i];*/
			else
				qfixzu += 0.+   qw[i]; /* qknoten - qneu[i]*wf +   qr[i];*/
		break;
		default:
		break;
	}
	}


	qzu = qfestzu+qfixzu;
	qab = qfestab+qfixab;
	qg  = qzu + qab;

	output = fopen(datei6,"wt");
	fprintf(output,"Berechnete Wasserbilanz (l/s)\n");
	fprintf(output,"-----------------------------\n\n");
	fprintf(output,"          Zufluesse  Abfluesse\n");

	fprintf(output," qfest :   %e     %e\n",   qfestzu*1000, qfestab*1000);
	fprintf(output," qfix  :   %e     %e\n",    qfixzu*1000,  qfixab*1000);
	fprintf(output,"--------------------------------\n");
	fprintf(output," Summe : %e   %e\n\n\n", qzu*1000, qab*1000);

	fprintf(output," Gesamt-Bilan z  : %e\n\n\n", qg*1000);
 	fclose(output);
	return(	nxmax);
}

