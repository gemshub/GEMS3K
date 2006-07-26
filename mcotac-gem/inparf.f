c/*==========================================================================*
c*
c*                            Unterprogramm: inpar_
c*                          Lesen der Eingangsparameter
c*
c* Input:
c*
c*    Datei     : Einlese-Datei
c*
c* Output:
c*
c*    ITest     :  wenn ITest=1 -> Protokoll der Eingabedaten
c*    NCyc      :  Anz. Zeitintervalle bzw. Pumpstufen
c*    NTim      :  Anz. Zeitschritte innerhalb eines Intervalls
c*    De        :  Dauer eines Zeitintervalls
c*    TMult     :  Zeitfaktor
c*    Nx        :  Anzahl der Knoten in X-Richtung
c*    Ny        :  Anzahl der Knoten in Y-Richtung
c
c*    isteu     :  1=Transport, 2=Bahnlinien, 3=Isochronen
c*    itest     :  wenn ITest=1 -> Protokoll der Eingabedaten
c*    inma      :  1 = Einsetzen der Partikel Åber inmat[][]
c*    ipfile    :  1 = Startpartikelpositionen werden aus einer Datei gelesen
c*    ntim      :  Anz. Zeitschritte innerhalb e. Intervalls
c*    de        :  Dauer eines Zeitschrittes
c*    npin      :  Anzahl Partikel im ersten Zeitschritt mit Kontamination
c
c*    nxmax        :  Anzahl der Knoten in X-Richtung
c*    nx        :  Anzahl der Knoten in Y-Richtung
c*    npkt      :  Anzahl der StÅtzstellen zur Definition der Eintragsfunktion
c*    teil      :  legt Ausschnitt fest mit Anfangs- u. Endknoten
c*    dx[]      :  Zellgrî·en in X-Richtung
c*    dy[]      :  Zellgrî·en in Y-Richtung
c*    backg     :  "background"-Wert der Schadstoffkonzentration
c*    rd        :  Retardierungsfaktor
c*    lambda    :  Zerfallskonstante
c*    lambt0    :  Verzîgerung bis zum Beginn der Zerfallsprozesse
c*    al        :  longitudinale DispersivitÑt
c*    vxx      :   const fliessgeschwindigkeit
c*    dm0       :   Diffusionskoeffizient
c*    dtmax    :   max time step
c*
c*   return =  Fehlernummer IErr
c*
*==========================================================================*/
	subroutine inparf(itest,ncyc,nxmax,isteu,inma,
     *ipfile,ntim,npin,npkt,ismooth,i_sorb,j_sorb,j_decay
     *,backg,rd,xlambda,aquer,along,vxx,dm0,dtmax,tmult,dx,de)
	
	
      implicit double precision (a-h,o-z)

	include 'gwheader.inc'

      dimension dx(NNODEX+2)
	character*9 dum9


      open(10, file="test11.itm")
	read(10,'(a8,1x,i10)')dum9,itest
      write(*,*)dum9,itest
	read(10,'(a8,1x,i10)')dum9,ncyc
      write(*,*)dum9,ncyc
	read(10,'(a8,1x,i10)')dum9,ntim
      write(*,*)dum9,ntim
 	read(10,'(a8,1x,g12.4)')dum9,de
      write(*,*)dum9,de
	read(10,'(a8,1x,g12.4)')dum9,tmult
      write(*,*)dum9,tmult
	read(10,'(a8,1x,i10)')dum9,nxmax
      write(*,*)dum9,nxmax
	read(10,'(a8,1x,i10)')dum9,nymax
      write(*,*)dum9,nymax
	read(10,'(a8,1x,g12.4)')dum9,dxx
 	read(10,'(a8,1x,g12.4)')dum9,dyy
      write(*,*)dum9,dyy
	read(10,'(a8,1x,i10)')dum9,isteu
      write(*,*)dum9,isteu
	read(10,'(a8,1x,i10)')dum9,inma
      write(*,*)dum9,inma
	read(10,'(a8,1x,i10)')dum9,ipfile
      write(*,*)dum9,ipfile
	read(10,'(a8,1x,i10)')dum9,npin
      write(*,*)dum9,npin
	read(10,'(a8,1x,i10)')dum9,npkt
      write(*,*)dum9,npkt
 	read(10,'(a8,1x,g12.4)')dum9,backg
      write(*,*)dum9,backg
	read(10,'(a8,1x,g12.4)')dum9,rd
      write(*,*)dum9,rd
	read(10,'(a8,1x,g12.4)')dum9,xlambda
      write(*,*)dum9,xlambda
	read(10,'(a8,1x,g12.4)')dum9,aquer
      write(*,*)dum9,aquer
	read(10,'(a8,1x,g12.4)')dum9,along
      write(*,*)dum9,along
	read(10,'(a8,1x,g12.4)')dum9,vxx
      write(*,*)dum9,vxx
	read(10,'(a8,1x,g12.4)')dum9,dm0
      write(*,*)dum9,dm0
	read(10,'(a8,1x,g12.4)')dum9,dtmax
      write(*,*)dum9,dtmax
	read(10,'(a8,1x,i10)')dum9,ismooth
      write(*,*)dum9,ismooth
	read(10,'(a8,1x,i10)')dum9,i_sorb
      write(*,*)dum9,i_sorb
	read(10,'(a8,1x,i10)')dum9,j_sorb
      write(*,*)dum9,j_sorb
	read(10,'(a8,1x,i10)')dum9,j_decay
      write(*,*)dum9,j_decay

	close (10)
c        /*  Zellgroessen in X- und Y-Richtung: dx[i] */
		do 5 i=1,nxmax+1
 5		   dx(i) =  dxx;
c/*                printf("iiiiiinnn %f\n",dx[i]); */
      return 
	end
  


