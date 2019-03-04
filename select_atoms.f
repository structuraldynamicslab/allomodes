c ==========================================================================================
c ==========================================================================================
c...    Originally written by MD for allostery applications:
c       Defines "binding sites" (pocket)
c
c	Input:
c		nat:		total number of atoms in the molecule
c		coordmol:	coordinates of the atoms
c		cutoff:		selection of atoms are based on a cutoff...
c				this is the value of the cutoff
c		ncenter:	number of centers for defining binding site
c		list_center	list of all "centers"
c	Output:
c		npocket:	number of atoms in pockets
c		corresp_pocket	for each atom in pocket, the corresponding number
c			  	in the list of all atoms
c
c ==========================================================================================
c ==========================================================================================
c
        subroutine extract_pocket(coordmol,natom_select,corresp_select,
     1		corresp_pocket,npocket,cutoff,ncenter,list_center)
c
c ==========================================================================================
c ==========================================================================================
c
        integer i,j,k,iat,jat,jcenter,natom_select
	integer	npocket,ncenter
	integer	corresp_pocket(*),corresp_select(*),list_center(*)
c
	real*8	distsq, diff, cutoff, cutoff2
	real*8	x(3),y(3)
        real*8	coordmol(*)
c
c ==========================================================================================
c 	Brute force method (this subroutine is only called once!): for each atom,
c	check all centers if they are close enough....
c ==========================================================================================
c
	cutoff2 = cutoff*cutoff
c
	npocket = 0
c
	do 300 i=1,natom_select
c
		iat = corresp_select(i)
		x(1) = coordmol(3*(iat-1) + 1)
		x(2) = coordmol(3*(iat-1) + 2)
		x(3) = coordmol(3*iat)
c
		do 200 j=1,ncenter
c
			jcenter = list_center(j)
			jat = corresp_select(jcenter)
			y(1) = coordmol(3*(jat-1) + 1)
			y(2) = coordmol(3*(jat-1) + 2)
			y(3) = coordmol(3*jat)
c
			distsq = 0.d0
			do 100 k = 1,3
				diff = x(k)-y(k)
				distsq = distsq + diff*diff
100			continue
c
			if(distsq.lt.cutoff2) then
				npocket = npocket + 1
				corresp_pocket(npocket) = iat
			endif
c
200		continue
c
300	continue
c
        return
        end
c
c ==========================================================================================
c ==========================================================================================
c...    Originally written by MD for allostery applications:
c       Defines "active sites" (channel)
c
c
c	Input:
c		nat:		total number of atoms in the molecule
c		coordmol:	coordinates of the atoms
c		nchannel:	number of active sites considered
c		list_center	list of residue range for the channels
c	Output:
c		nchannel:	number of atoms in channels
c		corresp_channel	for each atom in channels, the corresponding number
c			  	in the list of all atoms
c
c ==========================================================================================
c ==========================================================================================
c
        subroutine extract_channel(nat,corresp_channel,
     1		nchannel, ncenter, list_center)
c
c ==========================================================================================
c ==========================================================================================
c
        integer i,j,j1,j2,nat
	integer	nchannel, ncenter
	integer	corresp_channel(*),list_center(2,*)
c
	nchannel = 0
c
	do 200 i = 1,ncenter
c
		j1 = list_center(1,i)
		j2 = list_center(2,i)
c
		do 100 j = j1+1, j2-1
			nchannel = nchannel + 1
			corresp_channel(nchannel) = j
100		continue
c
200	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
        subroutine identify_region(coord,nat,corresp_select,region,
     1		nregion,noffset,corresp,nfound)
c
c ==========================================================================================
c ==========================================================================================
c
	integer	i,j,jat,k
	integer	nat,nregion,noffset,nfound,npos
c
	integer	corresp(*),corresp_select(*)
c
	real*8	tol,distsq,diff
	real*8	x(3),y(3)
	real*8	coord(*),region(*)
c
	tol = 1.d-2
c
	do 100 i = 1,noffset
		corresp(i) = corresp_select(i)
100	continue
c
	nfound = 0
c
	do 400 i = 1,nregion
c
		x(1) = region(3*(i-1) + 1)
		x(2) = region(3*(i-1) + 2)
		x(3) = region(3*i)
c
		do 300 j = noffset+1,nat
c
			jat = corresp_select(j)
			if(jat.eq.0) goto 300
c
			y(1) = coord(3*(jat-1) + 1)
			y(2) = coord(3*(jat-1) + 2)
			y(3) = coord(3*jat)
c
			distsq = 0.d0
			do 200 k = 1,3
				diff = x(k)-y(k)
				distsq = distsq + diff*diff
200			continue
c
			distsq = sqrt(distsq)
c
			if(distsq.lt.tol) then
				nfound = nfound + 1
				corresp(noffset+nfound) = jat
				corresp_select(j) = 0
				goto 400
			endif
c
300		continue
c
400	continue
c
	npos = noffset + nfound
	do 500 j = noffset+1,nat
		jat = corresp_select(j)
		if(jat.ne.0) then
			npos = npos + 1
			corresp(npos) = jat
		endif
500	continue
c
	do 600 i = 1,nat
		corresp_select(i) = corresp(i)
600	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
	subroutine read_atom(lun,nat,coord)
c
c ==========================================================================================
c ==========================================================================================
c
c	General variables
c ==========================================================================================
c
	integer	lun,j
	integer	nat
c
	real*8	x(3),coord(*)
c
	character record*80
c
1	format(a)
2 	format(3f8.3)
c
c ==========================================================================================
c
	nat = 0
c
100	read(lun,1,end=200) record
c
	if(record(1:4).eq.'ATOM') then
		nat = nat + 1
		read(record(31:53),2) (x(j),j=1,3)
		coord(3*(nat-1)+1) = x(1)
		coord(3*(nat-1)+2) = x(2)
		coord(3*nat)       = x(3)
	endif
	goto 100
c
200	continue
c		
	return
	end
c ==========================================================================================
c ==========================================================================================
c
	subroutine get_coord_region(coord,nat,corresp,nregion,region,
     1			flag)
c
c ==========================================================================================
c ==========================================================================================
c
	integer	i,iat
	integer	nat,nregion,flag
c
	integer	corresp(*)
c
	real*8	dist, dnrm2
	real*8	x(3)
	real*8	coord(*),region(*)
c
	call set_zeros(3*nat,region)
c
	do 200 i = 1,nregion
c
		iat = corresp(i)
c
		call dcopy(3,coord(3*(iat-1)+1),1,x,1)
		if(flag.eq.1) then
			dist = dnrm2(3,x,1)
			dist = 1.d0/dist
			call dscal(3,dist,x,1)
		endif
		call dcopy(3,x,1,region(3*(iat-1)+1),1)
c
200	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
	subroutine normalize(nat,coord)
c
c ==========================================================================================
c ==========================================================================================
c
	integer	i,nat
c
	real*8	dist, dnrm2
	real*8	coord(*)
c
	do 100 i = 1,nat
c
		dist = dnrm2(3,coord(3*(i-1)+1),1)
		dist = 1.d0/dist
		call dscal(3,dist,coord(3*(i-1)+1),1)
c
100	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
	subroutine def_corresp(nat,corresp,nat_select,corresp_select)
c
c ==========================================================================================
c ==========================================================================================
c
	integer	i,j
	integer	nat,nat_select
	integer	corresp(*),corresp_select(*)
c
	do 100 i = 1,nat
		corresp(i) = 0
100	continue
c
	do 200 i = 1,nat_select
		j = corresp_select(i)
		corresp(j) = i
200	continue
c
	return
	end
