c===================================================================================
c===================================================================================
c	Compute_bfact.f
c
c	This subroutine computes the atomic b-factors based on the eigenvalues and
c	eigenvectors of the Hessian matrix
c
c===================================================================================
c===================================================================================
c
	subroutine compute_bfact(nat,nmode1,nmode2,eigval,
     1		eigvect,bfact_calc,facb)
c
c===================================================================================
c===================================================================================
c
	integer	i,j
	integer	nat,nmode1,nmode2
c
	real*8	facb
	real*8	x,y,z,val
	real*8 	eigval(*),eigvect(3*nat,*)
	real*8	bfact_calc(*)
c
c===================================================================================
c	Initialize Bfact
c===================================================================================
c
	do 100 i = 1,nat
		bfact_calc(i) = 0.d0
100	continue
c
c===================================================================================
c	Compute first order contribution
c===================================================================================
c
	do 300 j = nmode1,nmode2
		do 200 i = 1,nat
			x = eigvect(3*i-2,j)
			y = eigvect(3*i-1,j)
			z = eigvect(3*i,j)
			val = (x*x+y*y+z*z)/eigval(j)
			bfact_calc(i) = bfact_calc(i) + val
200		continue
300	continue
c
c===================================================================================
c	Scale B-factors
c===================================================================================
c
	do 600 i = 1,nat
		bfact_calc(i) = facb*bfact_calc(i)
600	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	Compare_bfact.f
c
c	This subroutine computes the "distance" between computed and experimental
c	B-factors, where distance is either RMS or the Pearson's correlation coefficient
c
c===================================================================================
c===================================================================================
c
	subroutine compare_bfact(nat,corresp,bfact,bfact_calc,rms,
     1		correl,Knew)
c
c===================================================================================
c===================================================================================
c
	integer	i,iat
	integer	nat
	integer	corresp(*)
c
	real*8	rms,correl,Knew
	real*8	Se,See,Sec,Sc,Scc
	real*8	bfact(*),bfact_calc(*)
c
	Se  = 0.d0
	See = 0.d0
	Sc  = 0.d0
	Scc = 0.d0
	Sec = 0.d0
	do 100 i = 1,nat
		iat = corresp(i)
		Se  = Se + bfact(iat)
		See = See + bfact(iat)*bfact(iat)
		Sc  = Sc + bfact_calc(i)
		Scc = Scc + bfact_calc(i)*bfact_calc(i)
		Sec = Sec + bfact(iat)*bfact_calc(i)
100	continue
c
	Knew = Scc/Sec
c
	if(See.ne.0.d0.and.Scc.ne.0.d0) then
		rms = sqrt ( (See+Scc-2*Sec)/nat)
		correl = (nat*Sec-Se*Sc)/sqrt((nat*See-Se*Se)*
     1				(nat*Scc-Sc*Sc))
	else
		rms = 0.d0
		correl = 0.d0
	endif
c
	return
	end
c
c===================================================================================
c===================================================================================
c
	subroutine write_bfact(lun,nat,bfact,bfact_calc,label,corresp)
c
c===================================================================================
c===================================================================================
c
	integer	lun,nat,i,iat,npos
	integer	corresp(*)
c
	real*8	bfact(*),bfact_calc(*)
c
	character label(*)*30
c
1	format(i6,4x,f8.3,4x,f10.5,2x,a30)
c
	npos = 0
	do 100 i = 1,nat
		iat = corresp(i)
		if(iat.ne.0) then
			npos = npos + 1
			write(lun,1) npos,bfact(i),bfact_calc(iat),
     1				label(i)
		endif
100	continue
c
	return
	end
