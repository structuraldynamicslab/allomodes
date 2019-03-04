c============================================================================================
c============================================================================================
c	Listcontact_del.f		Version 1 6/6/2014	Patrice Koehl
c============================================================================================
c============================================================================================
c
c	This subroutine finds the list of contacts from the Delaunay triangulation
c
c	Copyright (C) 2014 Patrice Koehl
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
c============================================================================================
c
	subroutine listcontact_del(natom,coordmol)
c
c============================================================================================
c
	integer	i,j
	integer	ncont,nsize,natom
c
	integer ncontact(*),listcontact(2,*)
c
	real*8	coordmol(*),radius(*)
c
	pointer	(ptr_radius,radius)
	pointer (ptr_listcontact,listcontact)
	pointer	(ptr_ncontact,ncontact)
c
	common /contacts/       ptr_ncontact,ptr_listcontact,nsize,ncont
c
c=============================================================================
c	Set radius of each atom (point) to 1 arbitrarily
c=============================================================================
c
	ptr_radius = malloc(8*natom)
c
	do 100 i = 1,natom
		radius(i) = 1.d0
100	continue
c
c=============================================================================
c	Set up Delaunay calculation
c=============================================================================
c
	call setup(natom,coordmol,radius)
c
c=============================================================================
c	Compute weighted Delaunay triangulation
c =============================================================================
c
	call delaunay3D
c
c
c============================================================================================
c	Initialize array listcontact with "natom" edges...will be resized in subroutine
c	Delaunay_edges if necessary
c============================================================================================
c
	ncont = 0
	nsize = natom
	ptr_listcontact = malloc(2*natom*4)
c
	call Delaunay_edges
c
c============================================================================================
c	Find list of edges in Delaunay triangulation
c============================================================================================
c
	ptr_ncontact = malloc(4*natom)
c
	do 200 i = 1,natom
		ncontact(i) = 0
200	continue
c
	do 300 i = 1,ncont
c
c		Need to shift atom indices as Delaunay calculation introduced 4 dummy
c		points at positions 1-4
c
		listcontact(1,i) = listcontact(1,i) - 4
		listcontact(2,i) = listcontact(2,i) - 4
		j = listcontact(1,i)
		ncontact(j) = ncontact(j) + 1
300	continue
c
c=============================================================================
c	"Clean up" free space and remove arrays that are no more needed
c=============================================================================
c
	call cleanup
	call free(ptr_radius)
c
	return
	end
