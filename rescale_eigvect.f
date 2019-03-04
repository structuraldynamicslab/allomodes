c===================================================================================
c===================================================================================
c	Rescale_eigvect.f
c
c	For mass-weighted normal modes, rescale eigenvectors by mass
c
c===================================================================================
c===================================================================================
c
	subroutine rescale_eigvect(nat,corresp,sqrtmass,eigvect,nmodes)
c
	integer	i,j
	integer	nat,nmodes
	integer	corresp(*)
c
	real*8	fact
	real*8	sqrtmass(*),eigvect(3*nat,*)
c
	do 200 i = 1,nmodes
c
		do 100 j = 1,nat
			fact = 1.d0/sqrtmass(corresp(j))
			eigvect(3*(j-1)+1,i) = eigvect(3*(j-1)+1,i)*fact
			eigvect(3*(j-1)+2,i) = eigvect(3*(j-1)+2,i)*fact
			eigvect(3*(j-1)+3,i) = eigvect(3*(j-1)+3,i)*fact
100		continue
c
200	continue
c
	return
	end
