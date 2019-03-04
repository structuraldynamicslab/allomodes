c==========================================================================================
c==========================================================================================
c	Portable Random number generator.
c
c	From Numerical Recipes, W.H.Press, B.P.Flannery
c				S.A.Teukolsky, W.T.Vetterling
c		Cambridge Univ. Press
c
c	Function that returns a uniform random deviate between 0.0 and 1.0.
c	Set idum to any negative value to initialize the sequence.
c==========================================================================================
c==========================================================================================
c
	function ran2(idum)
C	-------------------
C
	integer	m,ia,ic
	real	rm,ran2
c
	parameter (m=714025,ia=1366,ic=150889,rm=1./m)
C
	integer	ir(97), iy,iff,idum,j
	data iff /0/
c
	save iy,ir
C
C	BEGIN.
C
	if(idum.lt.0.or.iff.eq.0) then
		iff=1
		idum=mod(ic-idum,m)
		do 11 j=1,97		! init the shuffle table
			idum=mod(ia*idum+ic,m)
			ir(j)=idum
11		continue
		idum=mod(ia*idum+ic,m)
		iy=idum
	endif
	j=1+(97*iy)/m
	if(j.gt.97.or.j.lt.1) then
		write(6,*) 'j=',j
		write(6,*) 'iy,m :',iy,m
		pause
	endif
	iy=ir(j)
C
C	RETURNED VALUE
C
	ran2=iy*rm
C
	idum=mod(ia*idum+ic,m)
	ir(j)=idum
C
	return
	end
