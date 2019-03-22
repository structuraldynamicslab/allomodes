c ==========================================================================================
c ==========================================================================================
c       
c	Score_allostery.f
c
c	This file contains a series of functions for scoring the deformations that
c	are related to allostery
c
c ==========================================================================================
c ==========================================================================================
c
c	1. Get_diff.f
c
c	Computes differences between the two conformations considered
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine get_diff(natom,corresp,coord,coord_target,diff)
c
	integer	natom,i,iat,j
	integer	corresp(*)
c
	real*8	coord(*),coord_target(*)
	real*8	diff(*)
c
	do 200 i = 1,natom
c
		iat = corresp(i)
		do 100 j = 1,3
			diff(3*(i-1)+j) = coord_target(3*(iat-1)+j)-
     1			coord(3*(i-1)+j)
100		continue
c
200	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
c	2. score_allostery.f
c
c	Computes the score: correlation between induced deformation on active site
c	atoms, and expected deformations
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine score_allostery(nat,nactive,nligand,diff,
     1		nmode1,nmode2,eigenval,eigenvect,work,score)
c
	integer	k
	integer	nat,nactive,nligand,nmode1,nmode2
c
	real*8	val,score,ddot
	real*8	work(*),diff(*)
	real*8	eigenval(*),eigenvect(3*nat,*)
c
c ==========================================================================================
c	1. Compute force from ligand side to active side, and correlate
c	   with expected displacement
c ==========================================================================================
c
	call set_zeros(3*nactive,work)
c
	do 100 k = nmode1, nmode2
c
		val = ddot(3*nligand,diff(3*nactive+1),1,
     1			eigenvect(3*nactive+1,k),1)
		val = val/eigenval(k)
c
		call daxpy(3*nactive,val,eigenvect(1,k),1,work,1)
c
100	continue
c
	call normalize(nactive,work)
c
	score = ddot(3*nactive,diff,1,work,1)
c
	score = score/nactive
	score = 1.d0 - score
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
c	2. dscore_allostery.f
c
c	Computes the score: correlation between induced deformation on active site
c	atoms, and expected deformations
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine dscore_allostery(iat,jat,nat,coord,nactive,nligand,
     1		disp,nmodes,nmode1,nmode2,eigenval,eigenvect,
     2		work1,work2,coef,deigvect,c,score,dscore)
c
	integer	i,j,k
	integer	iat,jat
	integer	nat,nactive,nligand,nsum,nmodes,nmode1,nmode2
c
	real*8	val,val1,score,ddot,dnrm2
	real*8	dval,dscore
	real*8	sum1
	real*8	dist,deig
	real*8	diff(3),vect1(3),vect2(3)
	real*8	coord(*),disp(*),work1(*),work2(*)
	real*8	eigenval(*),eigenvect(3*nat,*)
	real*8	coef(*),deigvect(*),c(*)
c
c ==========================================================================================
c	1. Compute force and dforce from ligand side, for each normal mode
c ==========================================================================================
c
	do 100 k = 1,3
		diff(k) = coord(3*(iat-1)+k)-coord(3*(jat-1)+k)
100	continue
c
	dist = ddot(3,diff,1,diff,1)
	dist = 1.d0/dist
c
	do 300 k = 1,nmodes
		sum1 = 0.d0
		do 200 i = 1,3
			sum1 = sum1 + diff(i)*(eigenvect(3*(iat-1)+i,k)
     1				 - eigenvect(3*(jat-1)+i,k))
200		continue
		coef(k) = sum1
300     continue
c
	call set_zeros(3*nactive,work1)
	call set_zeros(3*nactive,work2)
c
	nsum = nactive + nligand
c
	do 500 k = nmode1, nmode2
c
		deig = coef(k)*coef(k)*dist
c
		call set_zeros(3*nsum,deigvect)
c
		sum1 = coef(k)*dist
c
		do 400 i = 1,nmodes
			c(i) = coef(i)/(eigenval(k)-eigenval(i))
400		continue
		c(k) = 0.d0
c
		call dgemv('N',3*nsum,nmodes,sum1,eigenvect,nmodes,c,1,
     1		0.d0,deigvect,1)
c
		val = ddot(3*nligand,disp(3*nactive+1),1,
     1			eigenvect(3*nactive+1,k),1)
		dval = ddot(3*nligand,disp(3*nactive+1),1,
     1			deigvect(3*nactive+1),1)
c
		val = val/eigenval(k)
		dval = (dval - deig*val)/eigenval(k)
c
		call daxpy(3*nactive,val,eigenvect(1,k),1,work1,1)
		call daxpy(3*nactive,dval,eigenvect(1,k),1,work2,1)
		call daxpy(3*nactive,val,deigvect,1,work2,1)
c
500	continue
c
	do 600 j = 1,nactive
		call dcopy(3,work1(3*(j-1)+1),1,vect1,1)
		call dcopy(3,work2(3*(j-1)+1),1,vect2,1)
		val = dnrm2(3,vect1,1)
		val = 1.d0/val
		call dscal(3,val,vect1,1)
		call dscal(3,val,vect2,1)
		val1 = -ddot(3,vect1,1,vect2,1)
		call daxpy(3,val1,vect1,1,vect2,1)
		call dcopy(3,vect1,1,work1(3*(j-1)+1),1)
		call dcopy(3,vect2,1,work2(3*(j-1)+1),1)
600	continue
c
	score = ddot(3*nactive,disp,1,work1,1)
	dscore = ddot(3*nactive,disp,1,work2,1)
c
	score = score/nactive
	score = 1.d0 - score
	dscore = -dscore/nactive
c
	return
	end
