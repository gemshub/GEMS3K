/****************************************************************************
*
*		 		gwheader.h
*
*       Dimensionierungen
*
****************************************************************************/

/* Dimensionierung der Felder */

#define NCNODEX 510       /* Maximale Anzahl der Knoten in X-Richtung       */
#define NUPMAX 50000       /* Maximale Anzahl Partikel       */
#define NCPMAX 50000       /* Maximale Anzahl Partikel       */
#define NCBASIS 10
#define NCCOMPL 45
#define NCSOLID 20

/* min and max macros */

#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif


