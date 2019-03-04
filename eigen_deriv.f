c ==========================================================================================
c ==========================================================================================
c       
c	deigen.f
c
c	This subroutine computes the derivatives of the eigenvalues and eigenvectors
c	of the Hessian with respect to a elastic constant
c
c ==========================================================================================
c ==========================================================================================
c       
	subroutine deigen(nkval,nmode_tot,nmode1,nmode2,nat,coord,
     1		eigenval,eigenvect,deigenval,deigenvect,coef)
c       
c ==========================================================================================
c ==========================================================================================
c
	integer	nkval, nmode_tot, nmode1, nmode2, nat
	integer	nsize, ncont
	integer	i,j,k
	integer	iatom, jatom
	integer ncontact(*),listcontact(2,*)
c
	real*8	dist, val, sum1, sum2
	real*8	diff(3)
	real*8	coord(*),coef(*)
	real*8	eigenval(*),eigenvect(3*nat,*)
	real*8	deigenval(*),deigenvect(3*nat,*)
c
	pointer (ptr_ncontact,ncontact)
	pointer (ptr_listcontact,listcontact)
c
	common /contacts/       ptr_ncontact,ptr_listcontact,nsize,ncont
c
	iatom = listcontact(1,nkval)
	jatom = listcontact(2,nkval)
c
	do 100 k = 1,3
		diff(k) = coord(3*(iatom-1)+k)-coord(3*(jatom-1)+k)
100     continue
c
	call dotvect(diff,diff,dist)
	dist = 1.d0/dist
c
	do 400 k = 1,nmode_tot
		sum1 = 0.d0
		sum2 = 0.d0
		do 300 i = 1,3
			sum1 = sum1 + diff(i)*eigenvect(3*(iatom-1)+i,k)
			sum2 = sum2 - diff(i)*eigenvect(3*(jatom-1)+i,k)
300		continue
		coef(k) = sum1 + sum2
400	continue
c
	do 500 k = nmode1,nmode2
		deigenval(k) = coef(k)*coef(k)*dist
500	continue
c
	do 900 k = nmode1,nmode2
c
		do 600 j = 1,3*nat
			deigenvect(j,k) = 0.d0
600		continue
c
		sum1 = coef(k)*dist
c
		do 800 i = 1,nmode_tot
c
			if(i.eq.k) goto 800
c
			val = -sum1*coef(i)/(eigenval(i)-eigenval(k))
c
			do 700 j = 1,3*nat
				deigenvect(j,k) = deigenvect(j,k) + 
     1					val*eigenvect(j,i)
700			continue
c
800		continue
c
900	continue
c
	return
	end
