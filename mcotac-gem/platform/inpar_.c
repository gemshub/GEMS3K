/*==========================================================================*
*
*                            Unterprogramm: inpar_
*                          Lesen der Eingangsparameter
*
* Input:
*
*    Datei     : Einlese-Datei
*
* Output:
*
*    ITest     :  wenn ITest=1 -> Protokoll der Eingabedaten
*    NCyc      :  Anz. Zeitintervalle bzw. Pumpstufen
*    NTim      :  Anz. Zeitschritte innerhalb eines Intervalls
*    De        :  Dauer eines Zeitintervalls
*    TMult     :  Zeitfaktor
*    Nx        :  Anzahl der Knoten in X-Richtung
*    Ny        :  Anzahl der Knoten in Y-Richtung

*    isteu     :  1=Transport, 2=Bahnlinien, 3=Isochronen
*    itest     :  wenn ITest=1 -> Protokoll der Eingabedaten
*    inma      :  1 = Einsetzen der Partikel Åber inmat[][]
*    ipfile    :  1 = Startpartikelpositionen werden aus einer Datei gelesen
*    ntim      :  Anz. Zeitschritte innerhalb e. Intervalls
*    de        :  Dauer eines Zeitschrittes
*    npin      :  Anzahl Partikel im ersten Zeitschritt mit Kontamination

*    nxmax        :  Anzahl der Knoten in X-Richtung
*    nx        :  Anzahl der Knoten in Y-Richtung
*    npkt      :  Anzahl der StÅtzstellen zur Definition der Eintragsfunktion
*    teil      :  legt Ausschnitt fest mit Anfangs- u. Endknoten
*    dx[]      :  Zellgrî·en in X-Richtung
*    dy[]      :  Zellgrî·en in Y-Richtung
*    backg     :  "background"-Wert der Schadstoffkonzentration
*    rd        :  Retardierungsfaktor
*    lambda    :  Zerfallskonstante
*    lambt0    :  Verzîgerung bis zum Beginn der Zerfallsprozesse
*    al        :  longitudinale DispersivitÑt
*    vxx      :   const fliessgeschwindigkeit
*    dm0       :   Diffusionskoeffizient
*    dtmax    :   max time step
*
*   return =  Fehlernummer IErr
*
*   04.01.93
*==========================================================================*/

#include <stdio.h>
#include <math.h>
#include "gwheader.h"

int inpar_(itest,ncyc,tmult,nxmax,dx,isteu,inma,ipfile,
ntim,de,npin,npkt,backg,rd,xlambda,aquer,along,vxx,dm0,dtmax,ismooth,i_sorb,j_sorb,j_decay)

double *tmult,  *backg, *rd, *xlambda, *aquer, *along, *vxx, *dm0, *dtmax;
double dx[NCNODEX+2], *de ;
int *itest, *ncyc, *ntim, *nxmax, *isteu, *inma, *ipfile, *npkt, *ismooth, *i_sorb, *j_sorb;
int *j_decay;
int *npin;

{
        FILE *protokoll,*input;
        double xb, dxx, dyy;
        int ierr=0,  nymax;
        char name[80];

        register i, j;

        input=fopen("test11.itm","r");

        ierr += fscanf(input, "%*s %*s %d", itest);
                printf( "itest   : %d\n", *itest  );
        ierr += fscanf(input, "%*s %*s %d", ncyc);
                printf( "ncyc    : %d\n", *ncyc   );
        ierr += fscanf(input, "%*s %*s %d", ntim );
                printf("ntim    : %d\n", *ntim   );
        ierr += fscanf(input, "%*s %*s %lg", de );
                printf( "de      : %lg\n", *de     );
        ierr += fscanf(input, "%*s %*s %lg", tmult);
                printf( "tmult   : %lg\n", *tmult  );
        ierr += fscanf(input, "%*s %*s %d", nxmax );
                printf( "nxmax   : %d\n", *nxmax   );
        ierr += fscanf(input, "%*s %*s %d", &nymax );
                 printf( "nxmax   : %d\n", *nxmax   );
       ierr += fscanf(input, "%*s %*s %lg", &dxx );
                printf( "dxx     : %lf\n", dxx     );
        ierr += fscanf(input, "%*s %*s %lg", &dyy );
                printf( "dyy     : %lf\n", dyy     );
        ierr += fscanf(input, "%*s %*s %d", isteu);
                printf( "isteu   : %d\n", *isteu  );
        ierr += fscanf(input, "%*s %*s %d", inma);
                printf( "inma    : %d\n", *inma   );
        ierr += fscanf(input, "%*s %*s %d", ipfile);
                printf( "ipfile  : %d\n", *ipfile );
        ierr += fscanf(input, "%*s %*s %d", npin);
                printf( "npin    : %d\n", *npin   );
        ierr += fscanf(input, "%*s %*s %d", npkt);
                printf( "npkt    : %d\n", *npkt   );
        ierr += fscanf(input, "%*s %*s %lg", backg );
                printf( "backg   : %lg\n", *backg  );
        ierr += fscanf(input, "%*s %*s %lg", rd );
                printf( "rd      : %lg\n", *rd     );
        ierr += fscanf(input, "%*s %*s %lg", xlambda);
                printf( "xlambda : %lg\n", *xlambda );
        ierr += fscanf(input, "%*s %*s %lg", aquer);
                printf( "aquer   : %lg\n", *aquer );
        ierr += fscanf(input, "%*s %*s %lg", along );
                printf( "along   : %lg\n", *along  );
        ierr += fscanf(input, "%*s %*s %lg", vxx );
                printf( "vxx     : %lg\n", *vxx    );
        ierr += fscanf(input, "%*s %*s %lg", dm0 );
                printf( "dm0     : %lg\n", *dm0);
        ierr += fscanf(input, "%*s %*s %lf", dtmax);
                printf( "dtmax   : %lf\n", *dtmax);
        ierr += fscanf(input, "%*s %*s %d" , ismooth);
                printf( "ismooth : %d\n", *ismooth);
        ierr += fscanf(input, "%*s %*s %d" ,i_sorb);
                printf( "i_sorb : %d\n", *i_sorb );
        ierr += fscanf(input, "%*s %*s %d" ,j_sorb);
                printf( "j_sorb : %d\n", *j_sorb);
        ierr += fscanf(input, "%*s %*s %d" ,j_decay);
                printf( "j_decay: %d\n", *j_decay);
	fclose(input);

        if (*itest > 0)  {
                protokoll = fopen("tm.dat","w");
                fprintf(protokoll, "itest   : %d\n", *itest  );
                fprintf(protokoll, "ncyc    : %d\n", *ncyc   );
                fprintf(protokoll, "ntim    : %d\n", *ntim   );
                fprintf(protokoll, "de      : %lg\n", *de     );
                fprintf(protokoll, "tmult   : %lg\n", *tmult  );
                fprintf(protokoll, "nxmax   : %d\n", *nxmax   );
                fprintf(protokoll, "nxmax   : %d\n", *nxmax   );
                fprintf(protokoll, "dxx     : %lf\n", dxx     );
                fprintf(protokoll, "dyy     : %lf\n", dyy     );
                fprintf(protokoll, "isteu   : %d\n", *isteu  );
                fprintf(protokoll, "inma    : %d\n", *inma   );
                fprintf(protokoll, "ipfile  : %d\n", *ipfile );
                fprintf(protokoll, "npin    : %d\n", *npin   );
                fprintf(protokoll, "npkt    : %d\n", *npkt   );
                fprintf(protokoll, "backg   : %lg\n", *backg  );
                fprintf(protokoll, "rd      : %lg\n", *rd     );
                fprintf(protokoll, "xlambda : %lg\n", *xlambda );
                fprintf(protokoll, "aquer   : %lg\n", *aquer );
                fprintf(protokoll, "along   : %lg\n", *along  );
                fprintf(protokoll, "vxx     : %lg\n", *vxx    );
                fprintf(protokoll, "dm0     : %lg\n", *dm0);
                fprintf(protokoll, "dtmax   : %lf\n", *dtmax);
                fprintf(protokoll, "ismooth : %d\n", *ismooth);
                fprintf(protokoll, "i_sorb : %d\n", *i_sorb );
                fprintf(protokoll, "j_sorb : %d\n", *j_sorb);
                fprintf(protokoll, "j_decay: %d\n", *j_decay);
        }


        /*  Zellgroessen in X- und Y-Richtung: dx[i], dy[i] */
        for(i=0; i<=*nxmax-1; i++) {
                dx[i] =  dxx;
/*                printf("iiiiiinnn %f\n",dx[i]); */
        }
        if (itest)  fclose(protokoll);
        return(ierr==65 ? 0 : 2);
}



