c	kinetics.inc

	common /rkin1/amin,km,eab,eac,kmp,kmd,eabp,eabd,eacp,eacd,
     . 	       fg,omega,catinh,ratem,a,b,ap,ad,bp,bd,otherj,otheri,
     .	       dratemdc,dqdc,dxdc,rng,volmol,pi,por,texe,gespvfi,
     .	       gespvfb,phi

	common /rkin2/vout,rout,wout

	common /ikin1/ikin,ifgp,ifgd,nmineq

	real*8 amin(nsolid,nnodex,nnodey),km(nsolid)
	real*8 eab(nsolid,nbasis),eac(nsolid,ncompl)
	real*8 kmp(nsolid)
	real*8 eabp(nsolid,nbasis),eacp(nsolid,ncompl)
	real*8 kmd(nsolid)
	real*8 eabd(nsolid,nbasis),eacd(nsolid,ncompl)
	real*8 fg(nsolid),omega(nsolid),catinh(nsolid)
	real*8 ratem(nsolid),a(nsolid),b(nsolid)
	real*8 ap(nsolid),ad(nsolid),bp(nsolid),bd(nsolid)
	real*8 otherj(nbasis),otheri(ncompl)
	real*8 dratemdc(nsolid,nbasis)
	real*8 dqdc(nsolid,nbasis),dxdc(ncompl,nbasis)
	real*8 rng(nsolid),volmol(nsolid),por(nnodex+2,nnodey+2)
	real*8 gespvfi(nsolid),gespvfb(nsolid)

	real*8 vout(nsolid,nnodex),rout(nsolid,nnodex,nnodey)
	real*8 wout(nsolid,nnodex,nnodey)

	integer ikin(nsolid),ifg(nsolid)
	integer ifgp(nsolid),ifgd(nsolid)

