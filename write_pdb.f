c===========================================================================================
c===========================================================================================
c       Write_pdb.f
c
c	THis subroutine writes a PDB file with the bfactor corresponding to the strength
c	of the links coming to the corresponding atom
c
c===========================================================================================
c===========================================================================================
c
	subroutine write_pdb(lun,nat,natom,coord,label,corresp,
     1                  bfact,nlink)
c
	integer	i,j,k,nat,natom
	integer	iatom,jatom,lun
	integer	nsize,ncont
	integer corresp(*),nlink(*)
	integer	ncontact(*),listcontact(2,*)
c
	real*8	scale,bfact_max,occ
	real*8	coord(*),bfact(*)
	real*8	kval(*)
c
	character label(*)*30
c
	pointer	(ptr_kval,kval)
	pointer (ptr_ncontact,ncontact)
	pointer (ptr_listcontact,listcontact)
c
c===========================================================================================
	common /contacts/       ptr_ncontact,ptr_listcontact,nsize,ncont
	common /elastic/        ptr_kval
c===========================================================================================
c
1	format(a30,3f8.3,2f6.2)
c
c===========================================================================================
c	Compute sume of kval for each vertex
c===========================================================================================
c
	do 100 i=1,natom
		nlink(i)=0
		bfact(i)=0.d0
100	continue
c
	do 200 i=1,ncont
		iatom=listcontact(1,i)
		jatom=listcontact(2,i)
		bfact(iatom) = bfact(iatom) + kval(i)
		bfact(jatom) = bfact(jatom) + kval(i)
		nlink(iatom) = nlink(iatom) + 1
		nlink(jatom) = nlink(jatom) + 1
200     continue
c
	bfact_max = 0.d0
	do 300 i=1,natom
		if(nlink(i).ne.0) then
			bfact(i)=bfact(i)/nlink(i)
		else
			bfact(i)=0.d0
		endif
		bfact_max = max(bfact_max,bfact(i))
300	continue
c
	scale=100.d0/bfact_max
c
	do 400 i=1,natom
		bfact(i)=bfact(i)*scale
400	continue
c
	occ = 1.d0
	do 500 i = 1,nat
c
		k = corresp(i)
		if(k.ne.0) then
			write(lun,1) label(i),(coord(3*(i-1)+j),j=1,3),
     1			occ,bfact(k)
		endif
c
500	continue
c
	return
	end
