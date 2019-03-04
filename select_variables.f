c===================================================================================
c===================================================================================
c	get_ca.f
c
c	This subroutine extracts the CA from the whole molecule
c===================================================================================
c===================================================================================
c
	subroutine get_ca(nca, corresp)
c
c===================================================================================
c===================================================================================
c
	integer	i
	integer	nca,nat
	integer	natom,nseq,nchain
	integer	itype(*),chain(*),listatom(*),check(*)
	integer	corresp(*)
c
	real*8	coord(*)
c
	character	seq(*)*4
c
	pointer (ptr_coord,coord)
	pointer (ptr_check,check)
	pointer (ptr_listatom,listatom)
	pointer	(ptr_itype,itype)
	pointer	(ptr_chain,chain)
	pointer	(ptr_seq,seq)
c
c===================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /xyz/   ptr_coord,ptr_check,ptr_listatom,natom
c===================================================================================
c
	nca = 0
	nat = 0
	do 200 i = 1,nseq
		if(chain(i).ne.0) then
			nat = nat + listatom(i)
			goto 200
		endif
		nca = nca + 1
		corresp(nca) = nat + 3
		nat = nat + listatom(i)
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	get_backbone.f
c
c	This subroutine extracts the heavy atoms from the backbone of the whole 
c	molecule
c===================================================================================
c===================================================================================
c
	subroutine get_backbone(nback, corresp)
c
c===================================================================================
c===================================================================================
c
	integer	i,j
	integer	nback,nat
	integer	natom,nseq,nchain
	integer	itype(*),chain(*),listatom(*),check(*)
	integer	corresp(*)
c
	real*8	coord(*)
c
	character	seq(*)*4
	character	atom(*)*4
c
	pointer (ptr_coord,coord)
	pointer (ptr_check,check)
	pointer (ptr_listatom,listatom)
	pointer	(ptr_atom,atom)
	pointer	(ptr_itype,itype)
	pointer	(ptr_chain,chain)
	pointer	(ptr_seq,seq)
c
c===================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /xyz/   ptr_coord,ptr_check,ptr_listatom,natom
	common /names/ ptr_atom
c===================================================================================
c
	nback = 0
	nat = listatom(1)
	do 200 i = 2,nseq-1
		if(chain(i).ne.0) then
			nat = nat + listatom(i)
			goto 200
		endif
		do 100 j = 1,min(listatom(i),6)
			if(atom(nat+j)(1:1).ne.'H'.and.
     1			atom(nat+j)(1:1).ne.'X') then
				nback = nback + 1
				corresp(nback) = nat + j
			endif
100		continue
		nat = nat + listatom(i)
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	get_heavy.f
c
c	This subroutine extracts the heavy atoms from the whole molecule
c===================================================================================
c===================================================================================
c
	subroutine get_heavy(nheavy, corresp)
c
c===================================================================================
c===================================================================================
c
	integer	i,j
	integer	nheavy,nat
	integer	natom,nseq,nchain
	integer	itype(*),chain(*),listatom(*),check(*)
	integer	corresp(*)
c
	real*8	coord(*)
c
	character	seq(*)*4
	character	atom(*)*4
c
	pointer (ptr_coord,coord)
	pointer (ptr_check,check)
	pointer (ptr_listatom,listatom)
	pointer	(ptr_atom,atom)
	pointer	(ptr_itype,itype)
	pointer	(ptr_chain,chain)
	pointer	(ptr_seq,seq)
c
c===================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /xyz/   ptr_coord,ptr_check,ptr_listatom,natom
	common /names/ ptr_atom
c===================================================================================
c
	nheavy = 0
	nat = listatom(1)
	do 200 i = 2,nseq
		if(chain(i).ne.0.and.chain(i).ne.2) then
			nat = nat + listatom(i)
			goto 200
		endif
		do 100 j = 1,listatom(i)
			if(atom(nat+j)(1:1).ne.'H'.and.
     1			atom(nat+j)(1:1).ne.'X') then
				nheavy = nheavy + 1
				corresp(nheavy) = nat + j
			endif
100		continue
		nat = nat + listatom(i)
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	center_select.f
c
c	This subroutine centers the molecule on the center of gravity of selected atoms
c===================================================================================
c===================================================================================
c
	subroutine center_select(mass,nat,corresp)
c
c===================================================================================
c===================================================================================
c
	integer	i,j
	integer	nat,iat
	integer	natom
	integer	listatom(*),check(*)
	integer	corresp(*)
c
	real*8	mtot,center(3)
	real*8	coord(*),mass(*)
c
	pointer (ptr_coord,coord)
	pointer (ptr_check,check)
	pointer (ptr_listatom,listatom)
c
c===================================================================================
	common /xyz/   ptr_coord,ptr_check,ptr_listatom,natom
	common /cg/    center
c===================================================================================
c
	mtot = 0.d0
	do 100 i = 1,3
		center(i) = 0.d0
100	continue
c
	do 300 iat = 1,nat
c
		i = corresp(iat)
		do 200 j = 1,3
			center(j) = center(j) + mass(i)*coord(3*(i-1)+j)
200		continue
		mtot = mtot + mass(i)
c
300	continue
c
	do 400 i = 1,3
		center(i) = center(i)/mtot
400	continue
c
	do 600 i = 1,natom
		do 500 j = 1,3
			coord(3*(i-1)+j) = coord(3*(i-1)+j)-center(j)
500		continue
600	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	get_select.f
c
c	This subroutine extracts the selected atoms from the whole molecule
c===================================================================================
c===================================================================================
c
	subroutine get_select(nat, corresp, coord_at)
c
c===================================================================================
c===================================================================================
c
	integer	i,j
	integer	iat,nat
	integer	natom
	integer	listatom(*),check(*)
	integer	corresp(*)
c
	real*8	coord(*),coord_at(*)
c
	pointer (ptr_coord,coord)
	pointer (ptr_check,check)
	pointer (ptr_listatom,listatom)
c
c===================================================================================
	common /xyz/   ptr_coord,ptr_check,ptr_listatom,natom
c===================================================================================
c
	do 200 i = 1,nat
		iat = corresp(i)
		do 100 j = 1,3
			coord_at(3*(i-1)+j) = coord(3*(iat-1)+j)
100		continue
200	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	select_dihed.f		Version 1 4/29/2017	Patrice Koehl
c
c	This subroutine select dihedral angles so that each has a movable part and a non
c	movable part
c ==========================================================================================
c ==========================================================================================
c
	subroutine select_dihed(nang,map,nat,corresp_at,ndihed_free,
     1				corresp)
c
c ==========================================================================================
c ==========================================================================================
c
	integer	i,nang,nat,iang,nref,nv,natom
	integer n_move,n_not,iat,iv
	integer	ndihed_free
c
	integer	corresp(*),corresp_at(*)
	integer	listatom(*)
	integer listva(*)
	integer	check(*),typevar(*)
	integer	bond(4,*),p(*)
c
	integer*1	map(nang,*)
c
	real*8	coord(*)
c
	character dihed_names(*)*6
c
	pointer	(ptr_listva, listva)
	pointer	(ptr_coord,coord)
	pointer	(ptr_check,check)
	pointer	(ptr_listatom,listatom)
	pointer	(ptr_bond,bond)
	pointer	(ptr_p,p)
	pointer	(ptr_typevar,typevar)
	pointer	(ptr_dihed_names,dihed_names)
c
c ==========================================================================================
	common /defva/ ptr_listva
	common /xyz/   ptr_coord,ptr_check,ptr_listatom,natom
	common /defda/ ptr_typevar,ptr_bond,ptr_p,nref,nv
	common /angle_names/ ptr_dihed_names
c ==========================================================================================
c
	ndihed_free = 0
c
	do 200 iang = 2,nang-1
c
		if(dihed_names(iang).eq.'OMEGA') goto 200
		if(typevar(iang+1).eq.2) goto 200
		if(typevar(iang-1).eq.7) goto 200
c
		n_move = 0
		n_not  = 0
c
		do 100 iat = 1,nat
c
			i = corresp_at(iat)
			iv = listva(i)
c
			if(map(iang,iv).eq.1) then
				n_move = n_move + 1
			else
				n_not = n_not + 1
			endif
c
100		continue
c
		if(n_move.ne.0.and.n_not.gt.1) then
			ndihed_free = ndihed_free + 1
			corresp(ndihed_free) = iang
		endif
c
200	continue
c
	return
	end
