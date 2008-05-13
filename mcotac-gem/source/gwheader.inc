      real*8 rlarge
      parameter (nupmax = 50000)
      parameter (nnodex = 101)
      parameter (nbasis = 10)
      parameter (ncompl = 25)
      parameter (nsolid = 10)

      parameter (zero = 0.)
      parameter (rlarge = 1.0e+12)
c parameter for GEMS
	integer gsize1, gsize2, gsize3
c     Number of Dependent Components (DCs) in chemical system
      parameter (gsize1 = 17)
C Number of Independent Components (ICs) in chemical system
      parameter (gsize2 = 7)
C Number of Phases (PHs) in chemical system
      parameter (gsize3 = 20)

