      real*8 rlarge
      parameter (nupmax = 500000)
      parameter (nnodex = 510)
      parameter (nbasis = 10)
      parameter (ncompl = 45)
      parameter (nsolid = 20)

      parameter (zero = 0.)
      parameter (rlarge = 1.0e+12)

c parameter for GEMS
	integer gsize1, gsize2, gsize3
c     Number of Dependent Components (DCs) in chemical system
      parameter (gsize1 = 19)
C Number of Independent Components (ICs) in chemical system
      parameter (gsize2 =43)
C Number of Phases (PHs) in chemical system
      parameter (gsize3 = 20)

