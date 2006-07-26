
/*==========================================================================*
*
*                            Unterprogramm: wegdat1d
*               Speichern eines Feldes in die Datei Fname.dat
*
* Input:
*
*   
*    nx   : Anzahl der Knoten in X- bzw. in Y-Richtung
*    fname    : File-Name
*    hb[i][j]  : abzuspeicherdes Feld
*    text     : Bezeichnung des Feldes
*
*   return =  Fehlernummer ierr
*
*    23.03.95 
*===========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include "gwheader.h"

int wegdat1d_(nxmax,fname,hb,text)

double hb[NCNODEX+2];
int nxmax;
char text[10], fname[10];

{
	register i;
	int ierr=0;
	FILE *output;
	int sum=0;
        float faktor=1;
        int nxx=0;

/*	input = fopen(fname,"r"); */

	output = fopen(fname,"w");

        fprintf(output, "%d %g\n", (nxmax),  faktor );

		for(i=0; i<=(nxmax)-1; i++)  {
              	    sum += fprintf(output," %10.5e",hb[i] ); 
		}
		 fprintf(output," \n"); 
	fclose(output);
/*	if (sum != (nx+2))) ierr = 3; */
	return(ierr);
}
/*==========================================================================*
*
*                            Unterprogramm: holdat1d
*               Speichern eines Feldes in die Datei Fname.dat
*
*   
*    nx   : Anzahl der Knoten in X- bzw. in Y-Richtung
*    fname    : File-Name
*    hb[i][j] : einzulesendes Feld
*    text     : Bezeichnung des Feldes
*
*   return =  Fehlernummer ierr
*
*    14.08.95 
*===========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include "gwheader.h"
#include <math.h>

int holdat1d(int nxmax,char* fname,double hb[NCNODEX+2],char*  text)
/*
double hb[NCNODEX+2];
int *nxmax;
char text[10], *fname[10];
*/
{
	register i;
	int ierr=0;
        int ihb[NCNODEX+2];
	FILE *input;
	int sum=0;
        float faktor=1.;
        int nxx=0;
        printf("a datei nx faktor %s %d %d %g\n",fname, nxmax,NCNODEX, faktor );
	input = fopen(fname,"r");
        fscanf(input, "%d %g", &nxx, &faktor );
        printf("datei nx faktor %s %d %d %g\n",fname, nxmax,NCNODEX, faktor );
           for(i=0; i<= (nxmax)-1; i++)  {
		    /* sum += fread(&h[i][j],sizeof(float),1,input);*/
              sum += fscanf(input," %d",&ihb[i] ); 	
/*        printf("1i= %d  %d %d \n",i,ihb[i],sum ); */
	      hb[i] = faktor *  (double) ihb[i];          
/*        printf("2i= %d  %g %d \n",i,hb[i],sum );  */
	   }
	fclose(input);
/*	if (sum != (NXMAX+2)*(NYMAX+2)) ierr = 3; */
	return(ierr);
}


