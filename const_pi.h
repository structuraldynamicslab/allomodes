c	Const_pi.h		Version 1 8/29/2000		Patrice Koehl
c
c	This file contains the definition of some constants
c
	real*8	pi,twopi
	integer	ibuild_pi
c
	save	ibuild_pi,pi
c
	data ibuild_pi /0/
c
	if(ibuild_pi.eq.0) then
c
		pi = acos(-1.d0)
c
		ibuild_pi = 1
c
	endif
c
