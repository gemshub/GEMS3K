
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
*
*   return =  Fehlernummer ierr
*
*    23.03.95
*===========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gwheader.h"

int wegdat1d_( int nxmax, char* fname, double hb[NCNODEX+2],char *text)
{
	int i;
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
*
*   return =  Fehlernummer ierr
*
*    14.08.95
*===========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include "gwheader.h"
#include <math.h>

#ifdef __PGI
int holdat1d_(int nxmax,char* fname,double hb[NCNODEX+2])
#else
int holdat1d(int nxmax,char* fname,double hb[NCNODEX+2])
#endif
{
	register int i;
	int ierr=0;
        int ihb[NCNODEX+2];
	FILE *input;
	int sum=0;
        float faktor=1.;
        int nxx=0;
        printf("a datei nx faktor %s %d %d %g\n",fname, nxmax,NCNODEX, faktor );
	input = fopen(fname,"r");
        fscanf(input, "%d %g", &nxx, &faktor );
        printf("datei nx faktor \n%s10\n %d %d %g\n",fname, nxmax,NCNODEX, faktor );
           for(i=0; i< nxmax; i++)  {
		    /* sum += fread(&h[i][j],sizeof(float),1,input);*/
              sum += fscanf(input," %d",&ihb[i] );
/*        printf("1i= %d  %d %d \n",i,ihb[i],sum ); */
	      hb[i] = faktor *  (double) ihb[i];
	      /*        printf("2i= %d  %g %d \n",i,hb[i],sum ); */ 
	   }
	fclose(input);
/*	if (sum != (NXMAX+2)*(NYMAX+2)) ierr = 3; */
	return(ierr);
}

#ifdef __PGI
int vtkout_(int number,double time,int nxmax,int m1, int m2, int m3,double dxarr[NCNODEX+2],double bn[NCNODEX+2][NCBASIS], double cn[NCNODEX+2][NCCOMPL], double pn[NCNODEX+2][NCSOLID], double por[NCNODEX+2], double eh[NCNODEX+2], double ph[NCNODEX+2], char dumb[10*NCBASIS],char dumc[10*NCCOMPL],char dump[10*NCSOLID])
#else
int vtkout(int number,double time,int nxmax,int m1, int m2, int m3,double dxarr[NCNODEX+2],double bn[NCNODEX+2][NCBASIS], double cn[NCNODEX+2][NCCOMPL], double pn[NCNODEX+2][NCSOLID], double por[NCNODEX+2], double eh[NCNODEX+2], double ph[NCNODEX+2], char dumb[10*NCBASIS],char dumc[10*NCCOMPL],char dump[10*NCSOLID])
#endif
{
	int ierr=0, i, l;
        double dout;
	char vtk_file_name[80];
	FILE *output;
        char number_char[10], file_base_name[10], cdummy[11];

        sprintf(file_base_name,".vtk");

        sprintf(vtk_file_name,"run_");
        sprintf(number_char,"%i",number);
        strcat(vtk_file_name,number_char) ;
        strcat(vtk_file_name,file_base_name);

	output = fopen(vtk_file_name,"w");

	fprintf(output, "# vtk DataFile Version 2.0\n");
	fprintf(output, "vtk file written by mcotac-gems\n");
	fprintf(output, "ASCII\n");
	fprintf(output, "DATASET RECTILINEAR_GRID\n");
	fprintf(output, "DIMENSIONS %i 3 1\n",nxmax);
	fprintf(output, "X_COORDINATES %i double\n",nxmax);
	   dout=0.0;
           for(i=0; i< nxmax; i++)  {
                dout=dout+dxarr[i];
		fprintf(output," %g\n",dout);
	   }
	fprintf(output, "Y_COORDINATES 3 double\n");
	fprintf(output," -1.0 0.0 1.0 \n");
	fprintf(output, "Z_COORDINATES 1 double\n");
	fprintf(output," 0.0 \n");

	   fprintf(output,"POINT_DATA %i\n",nxmax*3);
	for (l=0;l<m1;l++){
	    memset(cdummy,'\0',11);
	   for (i=0;i<10;i++) { 
	       if ((i>1) && !(isgraph(dumb[l*10+i])))
		{i=10;} else {cdummy[i]=dumb[l*10+i];}
	    }
	   fprintf(output,"SCALARS %s double\n",cdummy);
	   fprintf(output,"LOOKUP_TABLE default\n");
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g \n",bn[i][l]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g  \n",bn[i][l]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g \n",bn[i][l]);
	   }
	}

/*	   fprintf(output,"POINT_DATA %i\n",nxmax); */
	for (l=0;l<m2;l++){
	   memset(cdummy,'\0',11);
	   for (i=0;i<10;i++) { 
	       if ((i>1) && !(isgraph(dumc[l*10+i])))
		{i=10;} else {cdummy[i]=dumc[l*10+i];}
	    }
	   fprintf(output,"SCALARS %s double\n",cdummy);
	   fprintf(output,"LOOKUP_TABLE default\n");
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",cn[i][l]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g \n",cn[i][l]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",cn[i][l]);
	   }
	}

/*	   fprintf(output,"POINT_DATA %i\n",nxmax); */
	for (l=0;l<m3;l++){
	   memset(cdummy,'\0',11);
	   for (i=0;i<10;i++) { 
	       if ((i>1) && !(isgraph(dump[l*10+i])))
		{i=10;} else {cdummy[i]=dump[l*10+i];}
	    }
	   fprintf(output,"SCALARS %s double\n",cdummy);
	   fprintf(output,"LOOKUP_TABLE default\n");
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",pn[i][l]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",pn[i][l]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g \n",pn[i][l]);
	   }
	}
	   fprintf(output,"SCALARS Eh double\n");
	   fprintf(output,"LOOKUP_TABLE default\n");
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",eh[i]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",eh[i]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g \n",eh[i]);
	   }	   
	   fprintf(output,"SCALARS pH double\n");
	   fprintf(output,"LOOKUP_TABLE default\n");
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",ph[i]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",ph[i]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g \n",ph[i]);
	   }	   
	   fprintf(output,"SCALARS por double\n");
	   fprintf(output,"LOOKUP_TABLE default\n");
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",por[i]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g\n",por[i]);
	   }
           for(i=0; i< nxmax; i++)  {
		fprintf(output," %g \n",por[i]);
	   }	   
	fclose(output);
	ierr=1;
	return(ierr);
}
