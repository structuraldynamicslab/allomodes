c===================================================================================
c===================================================================================
c	listcontact_cutoff.f
c
c       This set of subroutines computes the elastic network (contact map)
c	of a molecule, using a cutoff to find which atom pairs form an edge
c	in the network
c
c	Copyright (C) 2008 Patrice Koehl
c
c	This library is free software; you can redistribute it and/or
c	modify it under the terms of the GNU Lesser General Public
c	License as published by the Free Software Foundation; either
c	version 2.1 of the License, or (at your option) any later version.
c
c	This library is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c	Lesser General Public License for more details.
c
c	You should have received a copy of the GNU Lesser General Public
c	License along with this library; if not, write to the Free Software
c	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c
c===================================================================================
c===================================================================================
c
C	Contact.f		Version 1 14/5/1992	Patrice Koehl
c
c	This subroutine defines the list of "contacts" in the protein,
c	i.e. the list of pairs of atoms that are below a given cutoff
c	distance C
c
c===================================================================================
c===================================================================================
c	
	subroutine contact(coord,nat,C)
c
c===================================================================================
c	Variables used by defgrid
c===================================================================================
c
c	coord     contains the coordinates of the atoms considered
c	nat       is the number of atoms considered
c	C	  cutoff distance for a contact
c	ncontact  number of contacts for each atom
c	listcontact list of contact atoms for each atom
c
c===================================================================================
c	Declarations
c===================================================================================
c
	integer	natotc
	parameter	(natotc=150)
c
	integer nat,ncellx,ncelly,ncellz,ncelltot
	integer i,j,k,l,iat,jat,i1,j1,nb
	integer	ncont,nsize,npos,ncount
	integer	ndist
c
	integer natincell(*),listatincell(natotc,*)
c
	integer listb(100)
	integer ncontact(*),listcontact(2,*)
c
	real*8	distsq,C,Csq
	real*8	coord(3*nat),cube,xyzmin(3)
	real*8	dist_ref(*)
c
	pointer	(ptr_natincell,natincell)
	pointer	(ptr_listatincell,listatincell)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
	pointer	(ptr_dist_ref,dist_ref)
c
	common /grid/	  cube,xyzmin,ncellx,ncelly,ncellz
	common /grid2/	  ptr_natincell,ptr_listatincell
	common /contacts/ ptr_ncontact,ptr_listcontact,nsize,ncont
	common /distref/  ptr_dist_ref
c
	Csq = C*C
c
c===================================================================================
c	Define Cartesian grid for fast search of contacts
c===================================================================================
c
	cube = C
	call defgrid(coord,nat)
c
	ncelltot = ncellx*ncelly*ncellz
c
c===================================================================================
c	First count number of contacts, to initialize arrays
c===================================================================================
c
	ncont = 0
c
	do 600 i = 1,ncelltot
c
		if(natincell(i).eq.0) goto 600
c
c		First check inside cell
c
		do 200 j = 1,natincell(i)-1
			iat = listatincell(j,i)
			do 100 k = j+1,natincell(i)
				jat = listatincell(k,i)
				i1 = min(iat,jat)
				j1 = max(iat,jat)
				call distancesq(coord,i1,j1,distsq,
     1				3*nat)
				if(distsq.lt.Csq) then
					ncont = ncont + 1
				endif
100			continue
c
200		continue
c
c		Now check between cell i, and all neighbouring cells
c
		call neighbour_right(i,ncellx,ncelly,ncellz,nb,listb)
c
		do 500 j = 1,nb
			if(listb(j).eq.i) goto 500
c
			do 400 k = 1,natincell(listb(j))
				iat = listatincell(k,listb(j))
				do 300 l = 1,natincell(i)
					jat = listatincell(l,i)
					i1= min(iat,jat)
					j1= max(iat,jat)
					call distancesq(coord,i1,j1,
     1					distsq,3*nat)
					if(distsq.lt.Csq) then
					   ncont = ncont + 1
					endif
300				continue
400			continue
500		continue
c
600	continue
c
	nsize = ncont
	ptr_listcontact = malloc(8*nsize)
c
c===================================================================================
c	Find contacts within each cell, and between neighbouring cells
c===================================================================================
c
	ncont = 0
	ndist = 0
c
	do 1200 i = 1,ncelltot
c
		if(natincell(i).eq.0) goto 1200
c
c		First check inside cell
c
		do 800 j = 1,natincell(i)-1
			iat = listatincell(j,i)
			do 700 k = j+1,natincell(i)
				jat = listatincell(k,i)
				i1 = min(iat,jat)
				j1 = max(iat,jat)
				ndist = ndist + 1
				call distancesq(coord,i1,j1,distsq,
     1				3*nat)
				if(distsq.lt.Csq) then
					ncont = ncont + 1
					listcontact(1,ncont) = i1
					listcontact(2,ncont) = j1
				endif
700			continue
800		continue
c
c		Now check between cell i, and all neighbouring cells
c
		call neighbour_right(i,ncellx,ncelly,ncellz,nb,listb)
c
		do 1100 j = 1,nb
			if(listb(j).eq.i) goto 1100
			do 1000 k = 1,natincell(listb(j))
				iat = listatincell(k,listb(j))
				do 900 l = 1,natincell(i)
					jat = listatincell(l,i)
					i1= min(iat,jat)
					j1= max(iat,jat)
					ndist = ndist + 1
					call distancesq(coord,i1,j1,
     1					distsq,3*nat)
					if(distsq.lt.Csq) then
					   ncont = ncont + 1
					   listcontact(1,ncont) = i1
					   listcontact(2,ncont) = j1
					endif
900				continue
1000			continue
1100		continue
c
1200	continue
c
c===================================================================================
c	Now order the list of contacts
c===================================================================================
c
	call hpsort_two_int(listcontact,ncont)
c
	ptr_ncontact = malloc(4*nat)
	do 1300 i = 1,nat
		ncontact(i) = 0
1300	continue
c
	npos = 0
	do 1400 j = 1,ncont
		if(listcontact(1,j).ne.npos) then
			if(npos.eq.0) then
				npos = listcontact(1,j)
				ncount = 1
			else
				ncontact(npos) = ncount
				npos = listcontact(1,j)
				ncount = 1
			endif
		else
			ncount = ncount + 1
		endif
1400	continue
	ncontact(npos) = ncount
c
	npos = 0
	do 1500 i = 1,nat
		npos = npos + ncontact(i)
1500	continue
c
	ptr_dist_ref = malloc(8*ncont)
c
	do 1600 i = 1,ncont
		i1 = listcontact(1,i)
		j1 = listcontact(2,i)
		call distancesq(coord,i1,j1,distsq,3*nat)
		dist_ref(i) = sqrt(distsq)
1600	continue
c
	call free(ptr_natincell)
	call free(ptr_listatincell)
c
	return
	end
c
c===================================================================================
c===================================================================================
c
C	Defgrid.f		Version 1 14/5/1992	Patrice Koehl
c
c	This subroutine divides the volume which contains the molecule
c	in cubic cells to form a lattice.
c
c===================================================================================
c===================================================================================
c	
	subroutine defgrid(coord,nat)
c
c===================================================================================
c	Variables used by defgrid
c===================================================================================
c
c	coord     contains the coordinates of the atoms considered
c	nat       is the number of atoms considered
c	cube	  is the size of each small cube
c	natincell contains the number of atoms for each cell
c	listatincell contains the list of atoms for each cell
c	cell	  gives the cell number for each atom
c	xyzmax    is the upper point of the lattice
c	xyzmin    is the lower point of the lattice
c	ncellx    is the number of cells in the x dimension
c	ncelly    is the number of cells in the y dimension
c	ncellz    is the number of cells in the z dimension
c
c===================================================================================
c
c===================================================================================
c	Declarations
c===================================================================================
c
	integer	natotc
	parameter	(natotc=150)
c
	integer nat,ncellx,ncelly,ncellz,nplane,ncelltot
	integer natincell(*),listatincell(natotc,*)
	integer cell(*)
	integer i,j,ia,ja,ka,icell
c
	real*8	coord(3*nat),cube,xyzmin(3),xyzmax(3)
c
	pointer (ptr_cell,cell)
	pointer	(ptr_natincell,natincell)
	pointer	(ptr_listatincell,listatincell)
c
	common /grid/	cube,xyzmin,ncellx,ncelly,ncellz
	common /grid2/	ptr_natincell,ptr_listatincell
c
	ptr_cell = malloc(4*nat)
c
c===================================================================================
c	1. first position lattice on the molecule : find position of the
c	upper and lower point on the lattice
c===================================================================================
c
	do 100 i = 1,3
		xyzmin(i) = coord(i)
		xyzmax(i) = coord(i)
100	continue
c
	do 300 i = 2,nat
c
		do 200 j = 1,3
			if(coord(3*(i-1)+j).gt.xyzmax(j)) 
     1			xyzmax(j) = coord(3*(i-1)+j)
			if(coord(3*(i-1)+j).lt.xyzmin(j)) 
     1			xyzmin(j) = coord(3*(i-1)+j)
200		continue
c
300	continue
c
c	Expand size of the lattice on each side by 10 angstrom
c
	do 400 i = 1,3
		xyzmin(i) = xyzmin(i) - 10.d0
		xyzmax(i) = xyzmax(i) + 10.d0
400	continue
c
c	Define number of cells in each dimension, and total number
c	of cells
c
	ncellx = nint((xyzmax(1)-xyzmin(1))/cube) + 1
	ncelly = nint((xyzmax(2)-xyzmin(2))/cube) + 1
	ncellz = nint((xyzmax(3)-xyzmin(3))/cube) + 1
c
	nplane = ncellx*ncelly
	ncelltot = nplane*ncellz
c
	ptr_natincell = malloc(4*ncelltot)
	ptr_listatincell = malloc(4*natotc*ncelltot)
c
	if(ptr_natincell.eq.0.or.ptr_listatincell.eq.0) then
		write(6,*) ' '
		write(6,*) 'Problem in defgrid: could not allocate',
     1		' memory to arrays for grid. Please check'
		write(6,*) ' '
		stop
	endif
c
c	Give a cell number to each atom
c
	do 500 i = 1,nat
		ia = nint((coord(3*(i-1)+1) -xyzmin(1))/cube)+1
		ja = nint((coord(3*(i-1)+2) -xyzmin(2))/cube)+1
		ka = nint((coord(3*(i-1)+3) -xyzmin(3))/cube)+1
		icell = (ka-1)*nplane + (ja-1)*ncellx + ia
		cell(i) = icell
500	continue
c
c	Now, for each cell, give the list of atoms included
c
	do 600 i = 1,ncelltot
		natincell(i) = 0
		do 550 j = 1,natotc
			listatincell(j,i) = 0
550		continue
600	continue
c
	do 700 i = 1,nat
		natincell(cell(i)) = natincell(cell(i)) + 1
		if(natincell(cell(i)).eq.natotc) then
			write(6,*) 'too many atoms in a cell : ',
     1			natincell(cell(i))
			stop
		endif
		listatincell(natincell(cell(i)),cell(i)) = i
700	continue
c
	call free(ptr_cell)
c
	return
	end
c
c===================================================================================
c===================================================================================
c
c	Neighbour_right.f		Version 1 2/25/2001	Patrice Koehl
c
c	This subroutine finds the neighbours of a cell in a cube lattice
c	(two layers are considered, only to the right:  hence 42 neighbours 
c	are found)
c	Input	:
c			- ilayer	: number of layers
c			- icell		: cell number
c			- ncellx	: number of cells along x
c			- ncelly	: number of cells along y
c			- ncellz	: number of cells along z
c	Output 	:
c			- n		: number of neighbours 
c					  (including icell)
c			- list		: list of neighbours
c
c===================================================================================
c===================================================================================
c
	subroutine neighbour_right(icell,ncellx,ncelly,ncellz,n,list)
c
	integer	icell,ncellx,ncelly,ncellz,n,list(100)
	integer	nplane,ip,ix,iy,iz
	integer	ilayer,ik,ij,ii
	integer	i,j,k,istart,jstart
c
	nplane = ncellx*ncelly
	iz = icell/nplane + 1
	if(mod(icell,nplane).eq.0) iz = iz -1
c
	ip = icell - (iz-1)*nplane
	iy = ip/ncellx + 1
	if(mod(ip,ncellx).eq.0) iy = iy -1
c
	ix = ip - (iy-1)*ncellx
c
	ilayer = 2
c
	n = 0
	do 400 k = 0,ilayer
c
		ik = iz + k
		if(ik.gt.ncellz) goto 400
		if(k.eq.0) then
			jstart = 0
		else
			jstart = -ilayer
		endif
		do 300 j = jstart,ilayer
			ij = iy + j
			if(ij.le.0) goto 300
			if(ij.gt.ncelly) goto 300
			if(k.eq.0.and.j.eq.0) then
				istart = 0
			else
				istart = -ilayer
			endif
			do 200 i = istart,ilayer
				ii = ix + i
				if(ii.le.0) goto 200
				if(ii.gt.ncellx) goto 200
				n = n + 1
				list(n) = (ik-1)*nplane
     1				+ (ij-1)*ncellx + ii
200			continue
300		continue
c
400	continue
c
	return
	end
c===================================================================================
c===================================================================================
c
C	Contact_all.f		Version 1 14/5/1992	Patrice Koehl
c
c	This subroutine defines the list of "contacts" in the protein,
c	i.e. the list of all pairs of atoms
c
c===================================================================================
c===================================================================================
c	
	subroutine contact_all(coord,nat)
c
c===================================================================================
c	Variables used by defgrid
c===================================================================================
c
c	coord     contains the coordinates of the atoms considered
c	nat       is the number of atoms considered
c	ncontact  number of contacts for each atom
c	listcontact list of contact atoms for each atom
c
c===================================================================================
c	Declarations
c===================================================================================
c
	integer nat
	integer i,j,i1,j1
	integer	ncont,nsize,npos,ncount
c
	integer ncontact(*),listcontact(2,*)
c
	real*8	distsq
	real*8	coord(3*nat)
	real*8	dist_ref(*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
	pointer	(ptr_dist_ref,dist_ref)
c
	common /contacts/ ptr_ncontact,ptr_listcontact,nsize,ncont
	common /distref/  ptr_dist_ref
c
c===================================================================================
c	Define Cartesian grid for fast search of contacts
c===================================================================================
c
	nsize = (nat*(nat-1))/2
	ptr_listcontact = malloc(8*nsize)
c
	ncont = 0
	do 200 i = 1,nat - 1
c
		do 100 j = i+1,nat
c
			ncont = ncont + 1
			listcontact(1,ncont) = i
			listcontact(2,ncont) = j
c
100		continue
c
200	continue
c
	call hpsort_two_int(listcontact,ncont)
c
	ptr_ncontact = malloc(4*nat)
	do 300 i = 1,nat
		ncontact(i) = 0
300	continue
c
	npos = 0
	do 400 j = 1,ncont
		if(listcontact(1,j).ne.npos) then
			if(npos.eq.0) then
				npos = listcontact(1,j)
				ncount = 1
			else
				ncontact(npos) = ncount
				npos = listcontact(1,j)
				ncount = 1
			endif
		else
			ncount = ncount + 1
		endif
400	continue
	ncontact(npos) = ncount
c
	npos = 0
	do 500 i = 1,nat
		npos = npos + ncontact(i)
500	continue
c
	ptr_dist_ref = malloc(8*ncont)
c
	do 600 i = 1,ncont
		i1 = listcontact(1,i)
		j1 = listcontact(2,i)
		call distancesq(coord,i1,j1,distsq,3*nat)
		dist_ref(i) = sqrt(distsq)
600	continue
c
	return
	end
c
c===================================================================================
c===================================================================================
c
c	Distancesq.f		Version 1 14/5/1992	Patrice Koehl
c
c	This program calculates the distance**2 between two atoms n1 and n2
c	knowing the full vector of coordinates coord, of size nsize
c
c===================================================================================
c===================================================================================
c
	subroutine distancesq(coord,n1,n2,dist,nsize)
c
	integer	n1,n2,nsize,i
c
	real*8	coord(nsize),dist,p1(3),p2(3),p3(3)
c
	do 10 i = 1,3
		p1(i) = coord(3*(n1-1) + i)
		p2(i) = coord(3*(n2-1) + i)
		p3(i) = p2(i) - p1(i)
10	continue
c
	dist = 0
	do 20 i = 1,3
		dist = dist + p3(i)**2
20	continue
c
	return
	end
