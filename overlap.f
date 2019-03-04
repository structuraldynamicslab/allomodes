c ==========================================================================================
c ==========================================================================================
c       
c	Overlap.f
c
c	This subroutine computes the overlap between normal modes (described
c	by eigenvectors) and atomic displacements between two structures
c
c ==========================================================================================
c ==========================================================================================
c       
	subroutine compute_overlap(nmode1,nmode2,nmodes,nat,corresp,
     1		eigvect,coord_target,mass,overlaps,besto,besti)
c       
c ==========================================================================================
c ==========================================================================================
c       
	integer	i,j,iat,ierror
	integer	nmode1,nmode2,nmodes,natom,nat
	integer	besti
	integer	check(*),listatom(*)
	integer corresp(*)
c
	real*8	rms
	real*8	norm1,norm2,dotprod,dnrm2,ddot,besto
	real*8	coord_target(*),coord(*)
	real*8	mass(*)
	real*8	vect1(*),vect2(*)
	real*8	eigvect(3*nat,*)
	real*8	overlaps(*)
c
	pointer	(ptr_listatom,listatom)
	pointer	(ptr_check,check)
	pointer	(ptr_coord,coord)
	pointer	(ptr_vect1,vect1)
	pointer	(ptr_vect2,vect2)
c
c ==========================================================================================
	common /xyz/     ptr_coord,ptr_check,ptr_listatom,natom
c ==========================================================================================
c
	ptr_vect1 = malloc(8*3*nat)
	ptr_vect2 = malloc(8*3*nat)
c
	call bestfitm(coord,natom,coord_target,natom,nat,
     1		corresp,corresp,rms,ierror,mass,mass)
c
	do 200 i = 1,nat
		iat = corresp(i)
		do 100 j = 1,3
			vect1(3*(i-1)+j)=coord(3*(iat-1)+j)
     1				-coord_target(3*(iat-1)+j)
100		continue
200	continue
c
	norm1 = dnrm2(3*nat,vect1,1)
c
	do 300 i = 1,nmodes
		overlaps(i) = 0.d0
300	continue
c
	do 400 i = nmode1,nmode2
c
		call dcopy(3*nat,eigvect(1,i),1,vect2,1)
		norm2 = dnrm2(3*nat,vect2,1)
		dotprod = ddot(3*nat,vect1,1,vect2,1)
		overlaps(i) = abs(dotprod)/(norm1*norm2)
c
		if(i.eq.nmode1) then
			besti = i
			besto = overlaps(i)
		else
			if(overlaps(i).gt.besto) then
				besto = overlaps(i)
				besti = i
			endif
		endif
c
400	continue
c
	call free(ptr_vect1)
	call free(ptr_vect2)
c
	return
	end
c
c===================================================================================
c===================================================================================
c
	subroutine write_overlap(lun,nmode1,nmode2,overlap)
c
c===================================================================================
c===================================================================================
c
	integer	i,idx,lun,nmode1,nmode2
c
	real*8	over_tot,max_over,dnrm2
	real*8	cumul2,cumul
	real*8	overlap(*)
c
1	format(i6,4x,f8.3,4x,f8.3)
c
	cumul2 = 0.d0
	do 100 i = nmode1,nmode2
		cumul2 = cumul2 + overlap(i)*overlap(i)
		cumul = sqrt(cumul2)
		write(lun,1) i,overlap(i),cumul
		if(i.eq.nmode1) then
			idx = nmode1
			max_over = overlap(i)
		else
			if(max_over.lt.overlap(i)) then
				idx = i
				max_over = overlap(i)
			endif
		endif
100	continue
c
	over_tot = dnrm2(nmode2,overlap,1)
	write(lun,*) ' '
c
	write(lun,1) idx,max_over,over_tot
c
	return
	end
