c=================================================================================
c=================================================================================
c	Delaunay.f	
c
c	This file contains a library of fortran routines used to
c	compute the regular triangulation of a set of points in
c	3D space
c
c	Copyright (C) 2002 Patrice Koehl
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
c=================================================================================
c=================================================================================
C	Truncate_real.f		Version 1 10/15/2002	Patrice Koehl
c
c	This subroutine truncates a real number to a given accuracy
c	level (i.e. number of digits after the decimal point)
c
c=================================================================================
c=================================================================================
c
	subroutine truncate_real(x_in,x_out,ndigit)
c
c	Input:
c		x_in	: real number before truncation
c		ndigit	: number of digits to be kept
c	Output:
c		x_out	: real number after truncation
c
	real*8	x_in,x_out,y
	real*8	fact
c
	integer	i,mantissa
	integer	ndigit
	integer	digit(16)
c
	mantissa = int(x_in)
c
	y = x_in - mantissa
c
	x_out = mantissa
c
	fact = 1
	do 100 i = 1,ndigit
		fact = fact*10.d0
		digit(i) = nint(y*10)
		y = 10*(y-digit(i)/10.d0)
		x_out = x_out + digit(i)/fact
100	continue
c
	return
	end
c
c=====================================================================================
c=====================================================================================
c	adjust_nsphere.f		Version 1 6/17/2005	Patrice Koehl
c
c	This routine adjusts the number of
c	spheres included in the calculation to be always at least 4,
c	such that the regular triangulation can be computed, and the
c	alpha shape derived
c=====================================================================================
c
	subroutine adjust_nsphere
c
c=====================================================================================
c	Declaration
c=====================================================================================
c
	integer i,j,k
	integer	nvertex,npoint
	integer	sign(3,3)
c
	integer*1	redinfo(*)
c
c	Information on the vertices
c
	real*8	Rmax
	real*8	Dmax(3)
	real*8	radius(*)
	real*8	coord(*)
	real*8	coord4(*)
c
	data ((sign(i,j),j=1,3),i=1,3) /1,1,1,1,-1,-1,-1,1,1/
c
	pointer	(ptr_coord,coord)
	pointer (ptr_radius,radius)
	pointer (ptr_coord4,coord4)
	pointer (ptr_redinfo,redinfo)
c
	common /xyz_vertex/ ptr_coord,ptr_radius,ptr_coord4
	common /vertex_zone/ ptr_redinfo,npoint,nvertex
c
	save
c
c=====================================================================================
c	Do nothing if we already have at least 4 balls
c=====================================================================================
c
	if(npoint.ge.4) return
c
c=====================================================================================
c	Get bounding box of the current set of spheres, as well as max
c	radius
c=====================================================================================
c
	do 100 i = 1,3
		Dmax(i) = coord(i)
100	continue
	Rmax = radius(1)
c
	do 300 i = 2,npoint
		do 200 j = 1,3
			if(Dmax(j).lt.coord(3*(i-1)+j)) then
				Dmax(j) = coord(3*(i-1)+j)
			endif
200		continue
		if(Rmax.lt.radius(i)) Rmax = radius(i)
300	continue
c
c=====================================================================================
c	Now add point(s) with center at Dmax + 3*Rmax, and radius Rmax/20
c=====================================================================================
c
	do 500 i = npoint+1,4
		j = i - npoint
		do 400 k = 1,3
			coord(3*(i-1)+k)=sign(j,k)*(Dmax(k)+2*Rmax)
400		continue
		radius(i) = Rmax/20
500	continue
c
	npoint = 4
c
	return
	end
c=====================================================================================
c=====================================================================================
c	readjust_nsphere.f		Version 1 6/17/2005	Patrice Koehl
c
c	This subroutine removes the artificial spheres, if needed
c=====================================================================================
c=====================================================================================
c
	subroutine readjust_nsphere(nsphere,nred,listred)
c
c=====================================================================================
c	Declaration
c=====================================================================================
c
	integer i,j
	integer	npoints,nvertex,nsphere,nred
c
	integer	listred(nred)
c
c=====================================================================================
c	Information on the vertices
c=====================================================================================
c
	integer*1 redinfo(*)
c
	pointer (ptr_redinfo,redinfo)
c
	common  /vertex_zone/   ptr_redinfo,npoints,nvertex
c
	save
c
c=====================================================================================
c	Do nothing if we already have at least 4 balls
c=====================================================================================
c
	if(nsphere.ge.4) return
c
c=====================================================================================
c	Set artificial spheres as "redundant"
c=====================================================================================
c
	do 100 i = nsphere+5,8
		redinfo(i) = 1
100	continue
c
	npoints = nsphere
	nvertex = npoints + 4
c
	j = 0
	do 200 i = 1,nred
		if(listred(i).le.nsphere) then
			j = j + 1
			listred(j) = listred(i)
		endif
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	Setgmp.f
c===================================================================================
c===================================================================================
c
c	Copyright (C) 2007 Patrice Koehl
c 
c	This program sets the floating point filter for GMP calculation
c
c===================================================================================
	subroutine setgmp
c
c===================================================================================
c	Declaration
c===================================================================================
c
	real*8		SCALE_APA, EPS_APA
	real*8		scale,eps2,eps3,eps4,eps5
c
	common /gmp_info/	scale,eps2,eps3,eps4,eps5
c
	data SCALE_APA /100000000.d0/
	data EPS_APA   /1.0d0/
c
c===================================================================================
c	Set scale and precision for switching to GMP
c===================================================================================
c
	scale  = SCALE_APA
	eps2    = EPS_APA
	eps3    = EPS_APA
	eps4    = EPS_APA
	eps5    = EPS_APA
c
	return
	end
c
c===================================================================================
c===================================================================================
c	Setup.f
c===================================================================================
c===================================================================================
c
c	This subroutine gets the coordinates of the N points considered, and stores
c	these into the structures (common blocks) used in all suite of programs
c	Regular3D.f
c
c===================================================================================
c===================================================================================
c
	subroutine setup(nsphere,coord_sphere,radius_sphere)
c
c===================================================================================
c===================================================================================
c	Input:
c	*******
c
c	- nspheres:	number of points to be triangulated
c	- coord_sph:	coordinates of all points (in real*8)
c	- rad:   	weight of each point; this is the radius
c			of the sphere, while the "weight" usually
c			considered in regular triangulations
c			is the square of the radius
c===================================================================================
c
c===================================================================================
c	Declaration
c===================================================================================
c
	integer		ndigit,NDIGIT_APA
	integer		nsphere, npoint,nvertex,ntetra
	integer		nsize,nfudge,nsize_tet
	integer		i,j,k,ip,jp
c
        integer*8       nsize_tetr,ibyte8,ibyte16,ibyte24
        integer*8       nsize_r,nsize_r8,nsize_r16,nsize_r24
c
        integer*8       ival1,ival2,ival3
c
	integer*1	vertex_info(*)
	integer		ranlist(*)
c
	integer*1	tetra_info(*)
	integer*1	tetra_nindex(*)
c
	integer		tetra(4,*)
	integer		tetra_neighbour(4,*)
c
	real*8		scale,eps2,epsd,c_max,s,eps3,eps4,eps5
	real*8		x,xval,radius2
	real*8		y,z,w,xi,yi,zi,wi
	real*8		ranval(*)
	real*8		coord_sphere(*),coord(*)
	real*8		coord4(*)
	real*8		radius_sphere(*),radius(*)
c
	pointer		(ptr_vertex_info,vertex_info)
	pointer		(ptr_ranlist,ranlist)
	pointer		(ptr_ranval,ranval)
	pointer		(ptr_tetra_info,tetra_info)
	pointer		(ptr_tetra_nindex,tetra_nindex)
	pointer		(ptr_tetra,tetra)
	pointer		(ptr_tetra_neighbour,tetra_neighbour)
	pointer		(ptr_coord,coord)
	pointer		(ptr_radius,radius)
	pointer		(ptr_coord4,coord4)
c
	common	/xyz_vertex/	ptr_coord,ptr_radius,ptr_coord4
	common  /vertex_zone/	ptr_vertex_info,npoint,nvertex
	common /gmp_info/	scale,eps2,eps3,eps4,eps5
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
     	common /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common /tetra_size/	nsize_tet
c
	data NDIGIT_APA /8/
c
	save
c
	ndigit = NDIGIT_APA
c
	call setgmp
c
c===================================================================================
c	Define array sizes
c===================================================================================
c
	nfudge = 10
	nsize     = nsphere + nfudge
	nsize_tet = 10*nsize
        ibyte8    = 8
        ibyte16   = 16
        ibyte24   = 24
        nsize_r    = nsize
        nsize_tetr = nsize_tet
        nsize_tetr = nsize_tetr*ibyte16
        nsize_r24  = nsize_r*ibyte24
        nsize_r8   = nsize_r*ibyte8
c
	ptr_vertex_info     = malloc(nsize)
	ptr_ranlist         = malloc(4*nsize)
	ptr_ranval          = malloc(nsize_r24)
	ptr_tetra_info      = malloc(nsize_tet)
	ptr_tetra_nindex    = malloc(nsize_tet)
	ptr_tetra           = malloc(nsize_tetr)
	ptr_tetra_neighbour = malloc(nsize_tetr)
	ptr_coord	    = malloc(nsize_r24)
	ptr_radius	    = malloc(nsize_r8)
	ptr_coord4	    = malloc(nsize_r8)
c
c===================================================================================
c	Pre-processing:
c===================================================================================
c
c	1. Truncate input to NDIGIT precision
c
	npoint = nsphere
c
	c_max = 0
	do 100 i = 1,npoint
		vertex_info(i) = 0
		vertex_info(i)=ibset(vertex_info(i),0)
		x = radius_sphere(i)
		call truncate_real(x,xval,ndigit)
		radius(i) = xval
		do 50 j = 1,3
			k = 3*(i-1)+j
			x = coord_sphere(k)
			call truncate_real(x,xval,ndigit)
			coord(k) = xval
			if(abs(coord(k)).gt.c_max) c_max = abs(coord(k))
50		continue
100	continue
c
	c_max = max(100.d0,c_max)
c
c	2. Set eps
c
	epsd = 1.0d0
125	epsd = epsd / 2.0d0
	s = 1.0d0 + epsd
	if(s.ne.1.0d0) goto 125
	epsd = 2.0d0 * epsd
c
        eps2 = c_max**2 *epsd
        eps3 = eps2*c_max
        eps4 = eps3*c_max
        eps5 = eps4*c_max
	eps2 = 0.01
	eps3 = 0.01
	eps4 = 0.01
	eps5 = 0.01
c
c	eps5 = max(eps5,0.00000001d0)
c
c	3. Pre-compute all weights (stored in coord4):
c
	do 200 i = 1,npoint
                ival1 = nint(radius(i)*10000.d0)
                ival2 = -ival1*ival1
		do 150 j = 1,3
			k = 3*(i-1)+j
                        ival1 = nint(coord(k)*10000.d0)
                        ival2 = ival2 + ival1*ival1
150 		continue
                coord4(i) = dfloat(ival2)/100000000.d0
200	continue
c
c===================================================================================
c	Check for trivial redundant points: twice the same point
c===================================================================================
c
	do 300 i = 1,3*npoint
		ranval(i) = coord(i)
300	continue
	call hpsort_three(ranval,ranlist,npoint)
c
	jp = ranlist(1)
	x = coord(3*jp-2)
	y = coord(3*jp-1)
	z = coord(3*jp)
	w = radius(jp)
	do 400 i = 2,npoint
		ip = ranlist(i)
		xi = coord(3*ip-2)
		yi = coord(3*ip-1)
		zi = coord(3*ip)
		wi = radius(ip)
		if((xi-x)**2+(yi-y)**2+(zi-z)**2.le.100.d0*epsd) then
			if(wi.le.w) then
				vertex_info(ip) = ibclr(vertex_info(ip),0)
			else
				vertex_info(jp) = ibclr(vertex_info(jp),0)
				jp = ip
				w = wi
			endif
		else
			x = xi
			y = yi
			z = zi
			w = wi
			jp = ip
		endif
400	continue
c
c
c===================================================================================
c	Initialisation:
c===================================================================================
c
c	1. Add four infinite points
c
	do 500 i = npoint,1,-1
		coord4(i+4) = coord4(i)
		vertex_info(i+4) = vertex_info(i)
		radius(i+4) = radius(i)
500     continue
c
	do 600 i = 3*npoint,1,-1
		coord(i+12) = coord(i)
600	continue
c
	nvertex = npoint + 4
c
	do 700 i = 1,12
		coord(i) = 0
700	continue
c
	do 800 i = 1,4
		vertex_info(i) = 0
		vertex_info(i) = ibset(vertex_info(i),0)
		radius(i) = 0.d0
		coord4(i) = 0.d0
800	continue
c
c	2. Initialize arrays for GMP calculation
c
	call init_sos_gmp
c
c	3. Initialize tetra for Delaunay calculation
c
	ntetra = 1
	tetra(1,ntetra) = 1
	tetra(2,ntetra) = 2
	tetra(3,ntetra) = 3
	tetra(4,ntetra) = 4
c
	tetra_neighbour(1,ntetra) = 0
	tetra_neighbour(2,ntetra) = 0
	tetra_neighbour(3,ntetra) = 0
	tetra_neighbour(4,ntetra) = 0
c
	tetra_info(ntetra) = 0
c
	tetra_info(ntetra) = ibset(tetra_info(ntetra),1)
c
c	4. orientation is right most bit. bit = 0 means -1, bit = 1 means 1
c	The orientation of the first tetrahedron is -1:
c
	tetra_info(ntetra) = ibclr(tetra_info(ntetra),0)
c
c===================================================================================
c	Free space from temporary arrays
c===================================================================================
c
	call free(ptr_ranval)
	call free(ptr_ranlist)
c
	return
	end
c=====================================================================================
c=====================================================================================
c	resize_tet.f		Version 1 6/17/2005	Patrice Koehl
c
c	This subroutine resizes all arrays related to tetrahedron, if the initial
c	estimate of the number of tetrahedron was wrong
c=====================================================================================
c=====================================================================================
c
	subroutine resize_tet
c
c===================================================================================
c	All declarations
c===================================================================================
c
	integer		i,j
	integer		ntetra
	integer		nsize_tet
c
        integer*8       nsize_tetr,ibyte16
c
	integer*1	tetra_info(*)
	integer*1	tetra_nindex(*)
	integer*1	tetra_info2(*)
	integer*1	tetra_nindex2(*)
c
	integer		tetra(4,*)
	integer		tetra_neighbour(4,*)
	integer		tetra2(4,*)
	integer		tetra_neighbour2(4,*)
c
	pointer		(ptr_tetra_info,tetra_info)
	pointer		(ptr_tetra_nindex,tetra_nindex)
	pointer		(ptr_tetra,tetra)
	pointer		(ptr_tetra_neighbour,tetra_neighbour)

	pointer		(ptr_tetra_info2,tetra_info2)
	pointer		(ptr_tetra_nindex2,tetra_nindex2)
	pointer		(ptr_tetra2,tetra2)
	pointer		(ptr_tetra_neighbour2,tetra_neighbour2)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
     	common /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common /tetra_size/	nsize_tet
c
	save
c
	nsize_tet = 1.5*ntetra
	nsize_tet = max(nsize_tet,ntetra+1000)
c
        ibyte16   = 16
        nsize_tetr = nsize_tet
        nsize_tetr = nsize_tetr*ibyte16
c
	ptr_tetra_info2      = malloc(nsize_tet)
	ptr_tetra_nindex2    = malloc(nsize_tet)
	ptr_tetra2           = malloc(nsize_tetr)
	ptr_tetra_neighbour2 = malloc(nsize_tetr)
c
	do 200 i = 1,ntetra
		tetra_info2(i) = tetra_info(i)
		tetra_nindex2(i) = tetra_nindex(i)
		do 100 j = 1,4
			tetra2(j,i) = tetra(j,i)
			tetra_neighbour2(j,i) = tetra_neighbour(j,i)
100		continue
200	continue
c
	call free(ptr_tetra)
	call free(ptr_tetra_neighbour)
	call free(ptr_tetra_info)
	call free(ptr_tetra_nindex)
c
	ptr_tetra           = ptr_tetra2
	ptr_tetra_neighbour = ptr_tetra_neighbour2
	ptr_tetra_info      = ptr_tetra_info2
	ptr_tetra_nindex    = ptr_tetra_nindex2
c
	return
	end
c=================================================================================
c=================================================================================
c	Delaunay3D.f
c
c	Copyright (C) 2002 Patrice Koehl
c 
c	This program computes the regular triangulation of a set
c	of N weighted points in 3D, using the incremental flipping
c	algorithm of Edelsbrunner.
c
c	This implementation is based on the algorithm published in
c	H. Edelsbrunner and N.R. Shah, Algorithmica (1996) 15: 223-241
c
c	1) Algorithm:
c	*************
c
c	Briefly, the algorithm works as follows:
c
c	- first, a large tetrahedron initialises the program. All
c	four vertices of this tetrahedron are set at "infinite"
c
c	- All N points are added one by one.
c
c	- For each point:
c
c		- localize the tetrahedron in the current regular
c		triangulation that contains this point
c
c		- test if the point is redundant; if yes, remove
c
c		- If the point is not redundant, insert in the
c		tetrahedron : this is a "1-4" flip
c
c		- collect all "link facets" (i.e. all triangles
c		in tetrahedron containing the new point, that face
c		this new point) that are not regular.
c
c		- for each non-regular link facet, check if it
c		is "flippable". If yes, perform a "2-3", "3-2"
c		or "1-4" flip. Add new link facets in the list,
c		if needed.
c
c		- when link facet list if empty, move to next 
c		point
c
c	- Remove "infinite" tetrahedron, i.e. tetrahedron with
c	one vertice at "infinite"
c
c	- collect all remaining tetrahedron, define convex hull,
c	and exit.
c
c	2) Data structure:
c	******************
c
c	I maintain a minimal data structure that includes only 
c	the tetrahedrons of the triangulation (triangles,
c	edges and vertices are implicit).
c
c	For each tetrahedron, I store:
c
c	- the index of its four vertices
c	- pointers to its neighbours (4 maximum).
c		neighbor(i) is the tetrahedron that shares
c		all vertices of the tetrahedron considered, except i
c		(0 if the corresponding face is on the convex hull)
c	- its status: 1 "active" (i.e. part of the triangulation), 0 inactive
c	- its orientation
c
c
c       3/18/07: Modification to the program:
c               To save space, use a bit representation for:
c
c       Vertices:       remove infpoint array, as infinite points
c                       are between 1 and 4, always
c                       make vertex_info array integer*1, to save space
c
c       Tetrahedron:    define two arrays of integer*1:
c                       tetra_info stores:
c                               bit 0: orientation
c                               bit 1+2: status
c                               bit 3-6: surface info 
c					(one for each face of the tetrahedron)
c					surface info is a tag on the
c					face of the tetrahedron considered.
c					This tag can be set to represent
c					the convex hull of the molecule
c					for example, or to indicate if the
c					face belongs to a restricted Delaunay
c					(such as the one used in a skin
c					surface)
c                       tetra_nindex:
c				if tetra is (a,b,c,d), each of its
c				face is shared with another tetrahedron.
c				For example, face (b,c,d) is also a
c				face of tetrahedron (b,c,d,e). 
c				tetra_nindex gives the index of point e
c				in the lsit (b,c,d,e)
c				This index can take value 1, 2, 3, 4.
c				We store (index-1), which has value
c				00, 01, 10, 11 in bit.
c
c                               bit 0+1: face bcd
c                               bit 2+3: face acd
c                               bit 4+5: face abd
c                               bit 6+7: face abc
c
c
c	3) number representation:
c	**************************
c
c	I use double precision floating points. However, if one of
c	the geometricaly test becomes "imprecise", I switch to
c	arbitrary precision arithmetics (using the gmp package).
c
c	All primitives have therefore a floating point filter
c
c=================================================================================
c=================================================================================
c
	subroutine delaunay3D
c
c=================================================================================
c=================================================================================
c	Input:
c	*******
c
c	- npoint  :	number of points to be triangulated
c	- coord   :	coordinates of all points (in real*8)
c	- coord4  :	weight of each point; this is the radius
c			of the sphere, while the "weight" usually
c			considered in regular triangulations
c			is the square of the radius
c
c	Output:
c	********
c
c	- ntetra:	number of tetrahedron in the final
c			regular triangulation
c	- tetra:	for each tetrahedron, gives the position of its
c			four vertices in ascending order
c	- hull:		for each tetrahedron on the convex hull,
c			gives the index (local index, i.e. in 1-4)
c			of the vertex NOT on the convex hull
c
c	Both input and output are exchanged via common blocks
c
c=====================================================================================
c	Now declare all variables
c=====================================================================================
c
	integer		i,j,k
	integer		ival,iredundant,tetra_loc,tetra_last
	integer		npoint,nvertex
	integer		ntetra,iflag
	integer		iseed
	integer		nfree,nkill,n_new
	integer		iweight
	integer		nlink_facet
	integer		nfree_max,nkill_max,nlink_max,new_max
	integer		nsize_tet
c
	integer*1	ival1
c
	integer*1	vertex_info(*)
c
	integer*1	tetra_info(*)
	integer*1	tetra_nindex(*)
c
	integer		list_new(*)
	integer		tetra(4,*)
	integer		tetra_neighbour(4,*)
	integer		free_space(*),kill_space(*)
	integer		link_facet(*)
	integer		link_index(*)
c
	real*4		ran2,r
c
	real*8		scale,eps2,eps3,eps4,eps5
c
	real*8		coord(*)
	real*8		radius(*)
	real*8		coord4(*)
c
	pointer		(ptr_free,free_space)
	pointer		(ptr_kill,kill_space)
	pointer		(ptr_list,list_new)
	pointer		(ptr_link_facet, link_facet)
	pointer		(ptr_link_index,link_index)
c
        pointer         (ptr_vertex_info,vertex_info)
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
        pointer         (ptr_coord,coord)
        pointer         (ptr_radius,radius)
        pointer         (ptr_coord4,coord4)
c
c	for simplicity, information are stored in common blocks
c
        integer ntestred,nfailred,ntestreg,nfailreg
c
        common  /testmp/        ntestred,nfailred,ntestreg,nfailreg
c
	common	/xyz_vertex/	ptr_coord,ptr_radius,ptr_coord4
	common  /vertex_zone/	ptr_vertex_info,npoint,nvertex
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
     	common /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common /freespace/	ptr_free,ptr_kill,nfree,nkill
	common /gmp_info/	scale,eps2,eps3,eps4,eps5
	common /update/		ptr_list,n_new
	common /flags/		iweight
	common /link_zone/	ptr_link_facet,
     1				ptr_link_index,nlink_facet
	common /tetra_size/	nsize_tet
c
	save
c
        ntestred = 0
        nfailred = 0
        ntestreg = 0
        nfailreg = 0
c
c=====================================================================================
c	Set space for local arrays
c=====================================================================================
c
	nfree_max = 20000
	nkill_max = 20000
	new_max   = 20000
	nlink_max = 20000
c
	ptr_free       = malloc(4*nfree_max)
	ptr_kill       = malloc(4*nkill_max)
	ptr_list       = malloc(4*new_max)
	ptr_link_facet = malloc(4*2*nlink_max)
	ptr_link_index = malloc(4*2*nlink_max)
c
c	Initialise "free" space to 0
c
	nfree = 0
	n_new = 0
c
c=====================================================================================
c	Build regular triangulation
c=====================================================================================
c
c	Now loop over all points:
c
	tetra_last = -1
	iweight = 1
	iseed = -1
c
	do 400 i = 1,npoint
c
                ival = i+4
c
		n_new = 0
c
		if(.not.btest(vertex_info(ival),0)) goto  400
c
c		first locate the point in the list of known tetrahedra
c
		tetra_loc = tetra_last
		call locate_jw(iseed,ival,tetra_loc,iredundant)
c
c		If the point is redundant, move to next point
c
		if(iredundant.eq.1) then
			vertex_info(ival) = ibclr(vertex_info(ival),0)
			goto 400
		endif
c
c		Otherwise, add point to tetrahedron : 1-4 flip
c
		call flipjw_1_4(ival,tetra_loc,tetra_last)
c
c		Now scan link_facet list, and flip till list is empty
c
		call flipjw(tetra_last)
c
c		At this stage, I should have a regular triangulation
c		of the i+4 points (i-th real points+4 "infinite" points)
c		Now add another point
c
		if(ntetra.gt.0.9*nsize_tet) then
			call resize_tet
		endif
c
400	continue
c
c=====================================================================================
c	Reorder the tetrahedra, such that vertices are in increasing order
c=====================================================================================
c
	iflag = 1
	call reorder_tetra(iflag,n_new,list_new)
c
c=====================================================================================
c	I have the regular triangulation: I need to remove the
c	simplices including infinite points, and define the
c	convex hull
c=====================================================================================
c
	call remove_inf
c
c=====================================================================================
c	Now I peel off flat tetrahedra at the boundary of the DT
c=====================================================================================
c
c	call peel
c
c=====================================================================================
c	Define list of redundant points
c=====================================================================================
c
c	nredundant = 0
c	do 500 i = 1,npoint
c		if(.not.btest(vertex_info(i+4),0)) then
c			nredundant = nredundant + 1
c			list_redundant(nredundant) = i
c		endif
c500	continue
c
c=====================================================================================
c	Free space from local arrays
c=====================================================================================
c
	call free(ptr_list)
	call free(ptr_free)
	call free(ptr_kill)
	call free(ptr_link_facet)
	call free(ptr_link_index)
c
	return
	end
c
c=====================================================================================
c=====================================================================================
c	Locate_jw.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine locates the tetrahedron containing a new
c	point to be added in the triangulation
c
c	This implementation of the point location scheme
c	uses a "jump-and-walk" technique: first, N active
c	tetrahedra are chosen at random. The "distances" between
c	these tetrahedra and the point to be added are computed,
c	and the tetrahedron closest to the point is chosen as
c	a starting point. The program then "walks" from that tetrahedron
c	to the point, till we find a tetrahedron that contains
c	the point.
c	It also checks if the point is redundant in the current
c	tetrahedron. It it is, the search terminates.
c
c=====================================================================================
c
	subroutine locate_jw(iseed,ival,tetra_loc,iredundant)
c
c=====================================================================================
c	Input:
c
c	- ival:	index of the points to be located
c
c	Output:
c
c	- tetra_loc:	tetrahedron containing the point
c	- iredundant:	flag for redundancy: 0 is not redundant,
c			1 otherwise
c
c=====================================================================================
c	Define array size
c=====================================================================================
c
c	Declare variables
c
	integer	ival,itetra,iorient,idx,iseed,i
	integer	ntetra,tetra_in
	integer	tetra_loc,iredundant
	integer	a,b,c,d
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
c
	logical	test_in,test_red
c
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
c
	common /tetra_zone/ ptr_tetra,ptr_tetra_neighbour,ntetra
	common /tetra_stat/ ptr_tetra_info,ptr_tetra_nindex
c
	save
c
c=====================================================================================
c	Start at the root of the history dag: tetra(1)
c=====================================================================================
c
	iredundant = 0
c
	if(ntetra.eq.1) then
		tetra_loc = 1
		return
	endif
c
        if(tetra_loc.le.0) then
                do 50 i=ntetra,1,-1
                        if(btest(tetra_info(i),1)) then
                                itetra = i
                                goto 75
                        endif
50              continue
75              continue
        else
                itetra = tetra_loc
        endif
c
100	continue
c
	a = tetra(1,itetra)
	b = tetra(2,itetra)
	c = tetra(3,itetra)
	d = tetra(4,itetra)
c
	iorient = -1
	if(btest(tetra_info(itetra),0)) iorient = 1
c
	call inside_tetra_jw(ival,a,b,c,d,iorient,test_in,
     1		test_red,idx)
c
	if(test_in) goto 200
c
	itetra = tetra_neighbour(idx,itetra)
	goto 100
c
200	continue
c
	tetra_loc = itetra
c
c	Now that we have the tetrahedron (at a given layer of the
c	history dag), check if point is redundant
c
	if(test_red) iredundant = 1
c
	return
	end
c
c=====================================================================================
c=====================================================================================
c	Jump.f		Version 1 3/6/2002	Patrice Koehl
c
c	This subroutine picks N active tetrahedra at random,
c	computes the distance of the points to be inserted
c	to each of this tetrahedron, and selects the closest
c	one
c=====================================================================================
c=====================================================================================
c
	subroutine jump(iseed,ival,itetra)
c
c=====================================================================================
c	Declare variables
c=====================================================================================
c
	integer	i,j,k,ival
	integer	Ntry,N,Nkeep
	integer	iseed
	integer	ntetra,itetra
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
	integer	list(100)
c
	real*4	ran2,x
c
	real*8	dist,distmin
	real*8	coord(*),coord4(*)
	real*8	radius(*)
	real*8	xval(3),xa(3)
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
        pointer         (ptr_coord,coord)
        pointer         (ptr_radius,radius)
        pointer         (ptr_coord4,coord4)
c
	common  /xyz_vertex/ptr_coord,ptr_radius,ptr_coord4
	common /tetra_zone/ ptr_tetra,ptr_tetra_neighbour,ntetra
	common /tetra_stat/ ptr_tetra_info,ptr_tetra_nindex
c
	save
c
	Ntry = 20
c
	N = 200
c
50	continue
c
	Nkeep = 0
	do 100 i = 1,N
c
		x = ran2(iseed)
		j = int(ntetra*x) + 1
		j = min(j,ntetra)
		if(.not.btest(tetra_info(j),1)) goto 100
		do 75 k = 1,Nkeep
			if(j.eq.list(k)) goto 100
75		continue
		Nkeep = Nkeep + 1
		list(Nkeep) = j
		if(Nkeep.eq.Ntry) goto 125
c
100	continue
c
125	continue
c
	if(Nkeep.eq.0) goto 50
c
	do 200 i = 1,3
		xval(i) = coord(3*(ival-1)+i)
200	continue
c
	do 275 j = 1,4
		if(tetra(j,list(1)).le.4) goto 275
		do 250 i = 1,3
			xa(i) = coord(3*(tetra(j,list(1))-1)+i)
250		continue
		goto 285
275	continue
285	continue
c
	distmin = 0
	do 300 i = 1,3
		distmin = distmin + (xval(i)-xa(i))*
     1		(xval(i)-xa(i))
300	continue
c
	itetra = list(1)
c
	do 600 i = 2,Nkeep
c
		do 425 k = 1,4
			if(tetra(k,list(i)).le.4) goto 425
			do 400 j = 1,3
				xa(j) = coord(3*(tetra(k,list(i))-1)+j)
400			continue
			goto 450
425		continue
450		continue
c
		dist = 0
		do 500 j = 1,3
			dist = dist + (xval(j)-xa(j))*
     1			(xval(j)-xa(j))
500		continue
c
		if(dist.lt.distmin) then
			distmin = dist
			itetra = list(i)
		endif
600	continue
c
	return
	end
c=====================================================================================
c=====================================================================================
c	Inside_tetra_jw.f	Version 1 2/11/2002	Patrice Koehl
c
c	This subroutine tests if a point p is inside a tetrahedron
c	defined by four points (a,b,c,d) (with orientation "iorient")
c	If p is found inside the tetrahedron, it also checks if it
c	is redundant
c
c	Computation is done in floating point, but it is switched to
c	multiple precision if the result is imprecise
c
c=====================================================================================
c=====================================================================================
c
	subroutine inside_tetra_jw(p,a,b,c,d,iorient,is_in,redundant,
     1				ifail)
c
c=====================================================================================
c	Input:
c		- p	:	index of point to be checked
c		- a,b,c,d:	the four vertices of the tetrahedron
c		- iorient:	orientation of the tetrahedron
c
c	Output:
c		- is_in:	"true" if p in the tetrahedron, "false"
c				otherwise
c		- redundant	"true" is p is redundant
c		- ifail:	In case p is not inside the tetrahedron,
c				ifail gives the index of the face that
c				fails the orientation test
c
c=====================================================================================
c
c=====================================================================================
c	Declare all variables
c=====================================================================================
c
	integer	i,j,k,l,m
	integer	p,a,b,c,d
	integer	ia,ib,ic,id,ie,idx
	integer	ic1,ic5,ic1_k,ic1_l,sign,sign5,sign_k,sign_l
	integer	nswap,iswap,ninf
	integer	iorient,ifail
	integer	val
	integer iweight
c
	integer	list(4)
	integer inf4_1(4),sign4_1(4)
	integer inf4_2(4,4),sign4_2(4,4)
	integer	sign4_3(4)
	integer	inf5_2(4,4),sign5_2(4,4)
	integer inf5_3(4),sign5_3(4)
	integer	order1(3,4),order2(2,6),order3(2,6)
c
	integer*1	infpoint(4)
c
	real*8	Sij_1,Sij_2,Sij_3,Skl_1,Skl_2,Skl_3
	real*8	det_pijk,det_pjil,det_pkjl,det_pikl,det_pijkl
	real*8	eps2,eps3,eps4,eps5,scale
	real*8	detij(3)
	real*8	coordp(3),i_p(4),j_p(4),k_p(4),l_p(4)
c
	real*8	coord(*)
        real*8	coord4(*)
	real*8	radius(*)
c
        pointer         (ptr_coord,coord)
        pointer         (ptr_radius,radius)
        pointer         (ptr_coord4,coord4)
c
	logical	test_pijk,test_pjil,test_pkjl,test_pikl
	logical is_in,redundant
c
        integer ntestred,nfailred,ntestreg,nfailreg
c
	common  /xyz_vertex/	ptr_coord,ptr_radius,ptr_coord4
	common  /gmp_info/	scale,eps2,eps3,eps4,eps5
	common  /flags/		iweight
        common  /testmp/        ntestred,nfailred,ntestreg,nfailreg
c
	save
c
	data	inf4_1 /2,2,1,1/
	data	sign4_1 /-1,1,1,-1/
c
	data	inf4_2 /0,2,3,3,
     1                  2,0,3,3,
     2                  3,3,0,1,
     3                  3,3,1,0/
c
        data	sign4_2 / 0,1,-1,1,
     1                   -1,0,1,-1,
     2                    1,-1,0,1,
     3                   -1,1,-1,0/
c
	data	sign4_3 /-1,1,-1,1/
c
	data inf5_2	/0,2,1,1,
     1 			2,0,1,1,
     2			1,1,0,1,
     3			1,1,1,0/
c
	data sign5_2    /0,-1,-1,1,
     1			1,0,-1,1,
     2			1,1,0,1,
     3			-1,-1,-1,0/
c
	data inf5_3	/1,1,3,3/
	data sign5_3	/1,1,-1,1/
c
	data order1 /3,2,4,1,3,4,2,1,4,1,2,3/
	data order2/3,4,4,2,2,3,1,4,3,1,1,2/
	data order3/1,2,1,3,1,4,2,3,2,4,3,4/
c
c	If (i,j,k,l) is the tetrahedron in positive orientation, we need
c	to test:
c		(p,i,j,k)
c		(p,j,i,l)
c		(p,k,j,l)
c		(p,i,k,l)
c	If all four are positive, than p is inside the tetrahedron.
c	All four tests relies on the sign of the corresponding 4x4
c	determinant. Interestingly, these four determinants share
c	some common lines, which can be used to speed up the computation.
c
c	Let us consider or example:
c
c	det(p,i,j,k) = | p(1) p(2) p(3) 1|
c		       | i(1) i(2) i(3) 1|
c		       | j(1) j(2) j(3) 1|
c		       | k(1) k(2) k(3) 1|
c
c	p appears in each determinant. The corresponding line can therefore
c	be substraced from all 3 other lines . Using the example above,
c	we find:
c
c	det(i,j,k,l) = - |ip(1) ip(2) ip(3)|
c		         |jp(1) jp(2) jp(3)|
c			 |kp(1) kp(2) kp(3)|
c
c	where :	xp(m) = x(m) - p(m) for x = i,j,k and m = 1,2,3
c
c	Now we notice that the first two lines of det(p,i,j,k) and
c	det(p,i,j,l) are the same.
c
c	Let us define: Sij_3 = |ip(1) ip(2)| Sij_2 = |ip(1) ip(3)| 
c			       |jp(1) jp(2)|         |jp(1) jp(3)|
c	and Sij_1 = |ip(2) ip(3)|
c		    |jp(2) jp(3)|
c
c	We find:
c		det(p,i,j,k) = - kp(1)*Sij_1 + kp(2)*Sij_2 - kp(3)*Sij_3
c	and:
c		det(p,j,i,l) =   lp(1)*Sij_1 - lp(2)*Sij_2 + lp(3)*Sij_3
c
c	Similarly, if we define: 
c
c	Skl_3 = |kp(1) kp(2)|	Skl_2 = |kp(1) kp(3)|	Skl_1 = |kp(2) kp(3)|
c		|lp(1) lp(2)|		|lp(1) lp(3)|		|lp(2) lp(3)|
c
c	We find:
c		det(p,k,j,l) = jp(1)*Skl_1 - jp(2)*Skl_2 + jp(3)*Skl_3
c	and:
c		det(p,i,k,l) = -ip(1)*Skl_1 + ip(2)*Skl_2 - ip(3)*Skl_3
c
c	Furthermore:
c
c	det(p,i,j,k,l) = -ip(4)*det(p,k,j,l)-jp(4)*det(p,i,k,l)
c			 -kp(4)*det(p,j,i,l)-lp(4)*det(p,i,j,k)
c
c	The equations above hold for the general case; special care is
c	required to take in account infinite points (see below)
c
	is_in = .false.
	redundant = .false.
	
	list(1) = a
	list(2) = b
	list(3) = c
	list(4) = d
c
	infpoint(1) = 0
	infpoint(2) = 0
	infpoint(3) = 0
	infpoint(4) = 0
c
	if(a.le.4) infpoint(1) = 1
	if(b.le.4) infpoint(2) = 1
	if(c.le.4) infpoint(3) = 1
	if(d.le.4) infpoint(4) = 1
c
	ninf = infpoint(1) + infpoint(2) + infpoint(3) + infpoint(4)
c
c	"General case" : no infinite point
c
	do 100 m = 1,3
		coordp(m) = coord(3*p-3+m)
100	continue
c
	if(ninf.eq.0) then
c
c		Define coordinates (with i=a, j=b, k=c and l=d)
c		(no need to change notation, just bad habit to use i,j,k,l
c		instead of a,b,c,d !)
c
		do 200 m = 1,3
			i_p(m) = coord(3*a-3+m) - coordp(m)
			j_p(m) = coord(3*b-3+m) - coordp(m)
			k_p(m) = coord(3*c-3+m) - coordp(m)
			l_p(m) = coord(3*d-3+m) - coordp(m)
200		continue
c
c		Now compute 2x2 determinants Sij and Skl
c
		Sij_1 = i_p(2)*j_p(3) - i_p(3)*j_p(2)
		Sij_2 = i_p(1)*j_p(3) - i_p(3)*j_p(1)
		Sij_3 = i_p(1)*j_p(2) - i_p(2)*j_p(1)
c
		Skl_1 = k_p(2)*l_p(3) - k_p(3)*l_p(2)
		Skl_2 = k_p(1)*l_p(3) - k_p(3)*l_p(1)
		Skl_3 = k_p(1)*l_p(2) - k_p(2)*l_p(1)
c
c	Now perform tests
c
c	Start with is_in = .false. :
c
		is_in = .false.
c
		det_pijk = -k_p(1)*Sij_1 + k_p(2)*Sij_2 - k_p(3)*Sij_3
		det_pijk = det_pijk*iorient
		test_pijk = abs(det_pijk).gt.eps4
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = 4
			return
		endif
c
c	We check all other four determinants
c
		det_pjil = l_p(1)*Sij_1 - l_p(2)*Sij_2 + l_p(3)*Sij_3
		det_pjil = det_pjil*iorient
		test_pjil = abs(det_pjil).gt.eps4
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = 3
			return
		endif
c
		det_pkjl = j_p(1)*Skl_1 - j_p(2)*Skl_2 + j_p(3)*Skl_3
		det_pkjl = det_pkjl*iorient
		test_pkjl = abs(det_pkjl).gt.eps4
		if(test_pkjl.and.det_pkjl.gt.0) then
			ifail = 1
			return
		endif
c
		det_pikl = -i_p(1)*Skl_1 + i_p(2)*Skl_2 - i_p(3)*Skl_3
		det_pikl = det_pikl*iorient
		test_pikl = abs(det_pikl).gt.eps4
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = 2
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c		In this case, since we may need SoS, we have to rank
c		the indices
c
		if(.not.test_pijk) then
			call valsort4(p,a,b,c,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 4
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort4(p,b,a,d,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 3
				return
			endif
		endif
c
		if(.not.test_pkjl) then
			call valsort4(p,c,b,d,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 1
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort4(p,a,c,d,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 2
				return
			endif
		endif
c
c		If we have gone that far, p is inside the tetrahedron
c
		is_in = .true.
		if(iweight.eq.0) return
c
c		Now we check if p is redundant
c
		i_p(4) = coord4(a) - coord4(p)
		j_p(4) = coord4(b) - coord4(p)
		k_p(4) = coord4(c) - coord4(p)
		l_p(4) = coord4(d) - coord4(p)
c
		det_pijkl = -i_p(4)*det_pkjl - j_p(4)*det_pikl
     1			   -k_p(4)*det_pjil - l_p(4)*det_pijk
c
c	No need to multiply by iorient, since all minors contains iorient...
c
		if(abs(det_pijkl).lt.eps5) then
                        ntestred = ntestred + 1
			call valsort5(p,a,b,c,d,ia,ib,ic,id,
     1			ie,nswap)
			call sos_minor5_gmp(coord,radius,scale,ia,ib,ic,
     1			id,ie,val)
			det_pijkl = val*nswap*iorient
		endif
		redundant = det_pijkl.lt.0
c
	elseif(ninf.eq.1) then
c
c		We know that one of the 4 vertices a,b,c or d is 
c		infinite
c		To find which one it is, we use a map between
c		(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
c		is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4
c		if a,b,c or d are infinite, respectively.
c		A good mapping function is:
c		X = 3 - inf(a) - inf(a) -inf(b) + inf(d)
c
		idx=3-infpoint(1)-infpoint(1)-infpoint(2)+infpoint(4)
		l = list(idx)
c
c		The three finite points:
c
		i = list(order1(1,idx))
		j = list(order1(2,idx))
		k = list(order1(3,idx))
c
		ic1 = inf4_1(l)
		sign = sign4_1(l)
c
c	let us look at the four determinant we need to compute:
c
c	det_pijk	: unchanged
c	det_pjil	: 1 infinite point (l), becomes det3_pji
c			  where det3_pij = |p(ic1) p(ic2) 1|
c					   |i(ic1) i(ic2) 1|
c					   |j(ic1) j(ic2) 1|
c			  and ic1 and ic2 depends on which infinite
c			  (ic2 is always 3)
c			  point is considered
c	det_pkjl	: 1 infinite point (l), becomes det3_pkj
c	det_pikl	: 1 infinite point (l), becomes det3_pik
c
c	Get coordinates
c
		do 300 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			j_p(m) = coord(3*j-3+m) - coordp(m)
			k_p(m) = coord(3*k-3+m) - coordp(m)
300		continue
c
		detij(1) = i_p(1)*j_p(3) - i_p(3)*j_p(1)
		detij(2) = i_p(2)*j_p(3) - i_p(3)*j_p(2)
		detij(3) = i_p(1)*j_p(2) - i_p(2)*j_p(1)
c
c	Now perform tests
c
c	Start with is_in = .false. :
c
		is_in = .false.
c
		det_pijk = -k_p(1)*detij(2)+k_p(2)*detij(1)
     1				- k_p(3)*detij(3)
		det_pijk = det_pijk*iorient
		test_pijk = abs(det_pijk).gt.eps4
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = idx
			return
		endif
c
		det_pjil = -detij(ic1)*sign*iorient
		test_pjil = abs(det_pjil).gt.eps3
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = order1(3,idx)
			return
		endif
c
		det_pkjl = k_p(ic1)*j_p(3) - k_p(3)*j_p(ic1)
		det_pkjl = sign*det_pkjl*iorient
		test_pkjl = abs(det_pkjl).gt.eps3
		if(test_pkjl.and.det_pkjl.gt.0) then
			ifail = order1(1,idx)
			return
		endif
c
		det_pikl = i_p(ic1)*k_p(3) - i_p(3)*k_p(ic1)
		det_pikl = sign*det_pikl*iorient
		test_pikl = abs(det_pikl).gt.eps3
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = order1(2,idx)
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c
c
		if(.not.test_pijk) then
			call valsort4(p,i,j,k,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = idx
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort3(p,j,i,ia,ib,ic,nswap)
			call sos_minor3_gmp(coord,scale,ia,ib,ic,
     1				ic1,3,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order1(3,idx)
				return
			endif
		endif
c
		if(.not.test_pkjl) then
			call valsort3(p,k,j,ia,ib,ic,nswap)
			call sos_minor3_gmp(coord,scale,ia,ib,ic,
     1				ic1,3,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order1(1,idx)
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort3(p,i,k,ia,ib,ic,nswap)
			call sos_minor3_gmp(coord,scale,ia,ib,ic,
     1				ic1,3,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order1(2,idx)
				return
			endif
		endif
c
c	If we have gone so far, p is inside the tetrahedron
c
		is_in = .true.
c
c		Now we check if p is redundant
c
c		since det_pijkl = det_pijk >1
c		p cannot be redundant !
c
		redundant = .false.
c
	elseif(ninf.eq.2) then
c
c		We know that two of the 4 vertices a,b,c or d are
c		infinite
c		To find which one it is, we use a map between
c		(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
c		is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4,5,6
c		if (a,b), (a,c), (a,d), (b,c), (b,d), or (c,d) are
c		infinite, respectively
c		A good mapping function is:
c		X = 3 - inf(a) - inf(a) +inf(c) + inf(d) + inf(d)
c
		idx = 3 -infpoint(1) -infpoint(1) +infpoint(3)
     1			+ infpoint(4) + infpoint(4)
c
c		The two infinite points :
c
		k = list(order3(1,idx))
		l = list(order3(2,idx))
c
c		The two finite points
c
		i = list(order2(1,idx))
		j = list(order2(2,idx))
c
		ic1_k = inf4_1(k)
		ic1_l = inf4_1(l)
		sign_k = sign4_1(k)
		sign_l = sign4_1(l)
		ic1 = inf4_2(k,l)
		sign = sign4_2(k,l)
c
c	Get coordinates
c
		do 400 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			j_p(m) = coord(3*j-3+m) - coordp(m)
400		continue
c
c	Perform test; first set is_in .false.
c
		is_in = .false.
c
c	det_pijk is now det3_pij with k as infinite point
c
		det_pijk = i_p(ic1_k)*j_p(3)-i_p(3)*j_p(ic1_k)
		det_pijk = det_pijk*sign_k*iorient
		test_pijk = abs(det_pijk).gt.eps3
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = order3(2,idx)
			return
		endif
c
c	det_pjil is now det3_pji with l as infinite point
c
		det_pjil = i_p(3)*j_p(ic1_l)-i_p(ic1_l)*j_p(3)
		det_pjil = det_pjil*sign_l*iorient
		test_pjil = abs(det_pjil).gt.eps3
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = order3(1,idx)
			return
		endif
c
c	det_pkjl is now -det2_pj (k,l infinite)
c
		det_pkjl = j_p(ic1)*sign*iorient
		test_pkjl = abs(det_pkjl).gt.eps2
		if(test_pkjl.and.det_pkjl.gt.0) then
			ifail = order2(1,idx)
			return
		endif
c
c	det_pikl is now det2_pi (k,l infinite)
c
		det_pikl = -i_p(ic1)*sign*iorient
		test_pikl = abs(det_pikl).gt.eps2
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = order2(2,idx)
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c
		if(.not.test_pijk) then
			call valsort3(p,i,j,ia,ib,ic,nswap)
			call sos_minor3_gmp(coord,scale,ia,ib,ic,
     1				ic1_k,3,val)
			val = val*sign_k*nswap*iorient
			if(val.eq.1) then
				ifail = order3(2,idx)
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort3(p,j,i,ia,ib,ic,nswap)
			call sos_minor3_gmp(coord,scale,ia,ib,ic,
     1				ic1_l,3,val)
			val = val*sign_l*nswap*iorient
			if(val.eq.1) then
				ifail = order3(1,idx)
				return
			endif
		endif
c
		if(.not.test_pkjl) then
			call valsort2(p,j,ia,ib,nswap)
			call sos_minor2_gmp(coord,scale,ia,ib,ic1,val)
			val = -val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order2(1,idx)
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(coord,scale,ia,ib,ic1,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order2(2,idx)
				return
			endif
		endif
c
c	Again, if we have gone so far, p is inside the tetrahedron
c
		is_in = .true.
		redundant = .false.
		if(iweight.eq.0) return
c
c		Now we check if p is redundant
c
c		det_pijkl becomes det3_pij
c
		ic5 = inf5_2(k,l)
		sign5 = sign5_2(k,l)
		det_pijkl = i_p(ic5)*j_p(3)-i_p(3)*j_p(ic5)
		if(abs(det_pijkl).lt.eps3) then
			call valsort3(p,i,j,ia,ib,ic,nswap)
			call sos_minor3_gmp(coord,scale,ia,ib,ic,
     1			ic5,3,val)
			det_pijkl = val*nswap
		endif
		det_pijkl = det_pijkl*sign5*iorient
c
		redundant = det_pijkl.lt.0
c
	elseif(ninf.eq.3) then
c
c		We know that three of the 4 vertices a,b,c or d are
c		infinite
c		To find which one is finite, we use a map between
c		(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
c		is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4
c		if a,b,c or d are finite, respectively.
c		A good mapping function is:
c		X = 1 + inf(a) + inf(a) +inf(b) - inf(d)
c
		idx=1+infpoint(1)+infpoint(1)+infpoint(2)-infpoint(4)
		i = list(idx) 
		j = list(order1(1,idx))
		k = list(order1(2,idx))
		l = list(order1(3,idx))
c
c	Index of the "missing" infinite point (i.e. the fourth infinite
c	point)
c
		call missinf_sign(j,k,l,ie,iswap)
c
c	Get coordinates
c
		do 500 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
500		continue
c
c	Perform test; first set is_in to .false.
c
		is_in = .false.
c
c	det_pijk is now - det2_pi (missing j,k)
c
		det_pijk = i_p(inf4_2(j,k))*iorient*sign4_2(j,k)
		test_pijk = abs(det_pijk).gt.eps2
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = order1(3,idx)
			return
		endif
c
c	det_pjil is now det2_pi (missing j,l)
c
		det_pjil = -i_p(inf4_2(j,l))*iorient*sign4_2(j,l)
		test_pjil = abs(det_pjil).gt.eps2
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = order1(2,idx)
			return
		endif
c
c	det_pkjl is now det1_p
c
		det_pkjl = iorient*iswap*sign4_3(ie)
		if(det_pkjl.gt.0) then
			ifail = idx
			return
		endif
c
c	det_ikl is now - det2_pi (missing k,l)
c
		det_pikl = i_p(inf4_2(k,l))*iorient*sign4_2(k,l)
		test_pikl = abs(det_pikl).gt.eps2
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = order1(1,idx)
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c
c
		if(.not.test_pijk) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(coord,scale,ia,ib,
     1				inf4_2(j,k),val)
			val = -val*sign4_2(j,k)*iorient*nswap
			if(val.eq.1) then
				ifail = order1(3,idx)
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(coord,scale,ia,ib,
     1				inf4_2(j,l),val)
			val = val*sign4_2(j,l)*iorient*nswap
			if(val.eq.1) then
				ifail = order1(2,idx)
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(coord,scale,ia,ib,
     1				inf4_2(k,l),val)
			val = -val*sign4_2(k,l)*iorient*nswap
			if(val.eq.1) then
				ifail = order1(1,idx)
				return
			endif
		endif
c
		is_in = .true.
		redundant = .false.
		if(iweight.eq.0) return
c
c	Now check for redundancy
c
c		det_pijkl becomes -det2_pi
c
		ic1 = inf5_3(ie)
		sign5 = sign5_3(ie)
		det_pijkl = -i_p(ic1)
		if(abs(det_pijkl).lt.eps2) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(coord,scale,ia,ib,ic1,val)
			det_pijkl = val*nswap
		endif
		det_pijkl = - iorient*det_pijkl*sign5*iswap
c
		redundant = det_pijkl.lt.0
c
	else
c
c	In the case all four points ia,ib,ic, and id are infinite,
c	then is_in = .true. and redundant = .false.
c
		is_in = .true.
		redundant = .false.
c
	endif
c
	return
	end
c===========================================================================================
c===========================================================================================
c	Regular_convex.f	Version 1 2/25/2002	Patrice Koehl
c
c	This subroutine checks if a link facet (a,b,c) is locally
c	regular, as well as if the union of two tetrahedra
c	(a,b,c,p) and (a,b,c,o) that connects to the facet is convex. 
c
c	Computation is done in floating point, but it is switched to
c	multiple precision if the result is imprecise
c
c	In floating point, there is no need to order points in lexicographic
c	order prior to computing a determinant. This simplification is no
c	more true if we need to apply SoS: consequently, as soon as
c	the program swtiches to GMP (i.e. multi-precision), I first
c	order the points, using a series of routines valsort*, where
c	* can be 2,3,4 or 5, all contained in the file valsort
c
c===========================================================================================
c===========================================================================================
c
	subroutine regular_convex(a,b,c,p,o,itest_abcp,
     1			regular,convex,test_abpo,test_bcpo,test_capo)
c
c===========================================================================================
c
c	Input:
c		- a,b,c:	the three points defining the
c				link facet
c		- p:		the current point inserted in the
c				DT
c		- o:		the fourth point of the tetrahedron
c				that attaches to (a,b,c), opposite
c				to the tetrahedron (a,b,c,p,o)
c		- itest_abcp	orientation of the tetrahedron
c				(abcp)
c
c	Output:
c		- convex	"true" of (abcp)U(abco) is convex, "false"
c				otherwise
c		- regular	"true" if (abc) is locally regular, in which
c				case it does not matter if convex!
c
c===========================================================================================
c
c===========================================================================================
c	Declare al variables
c===========================================================================================
c
	integer	p,a,b,c,o,i,j,k,l,m
	integer	ia,ib,ic,id,ie
	integer	ninf,infp,info,iswap,iswap2,idx,val
	integer	icol1,sign1,icol2,sign2,icol4,sign4,icol5,sign5
	integer	itest_abcp
c
	integer	list(3)
	integer	order(2,3)
	integer inf4_1(4),sign4_1(4)
	integer inf4_2(4,4),sign4_2(4,4)
	integer sign4_3(4)
	integer	inf5_2(4,4),sign5_2(4,4)
	integer inf5_3(4),sign5_3(4)
	integer	order1(3,3)
c
	integer*1	infpoint(4)
c
	real*8	eps2,eps3,eps4,eps5,scale
	real*8	det_abpo,det_bcpo,det_capo,det_abcpo,det_abpc
	real*8	a_p(4),b_p(4),c_p(4),o_p(4)
	real*8	i_p(3),j_p(3)
	real*8	Mbo(3), Mca(3),Mjo(3),Mio(3)
	real*8	coordp(3)
c
	real*8	coord(*),radius(*)
        real*8	coord4(*)
c
	logical convex,regular,test_abpo,test_bcpo,test_capo
	logical testc(3)
c
        pointer         (ptr_coord,coord)
        pointer         (ptr_radius,radius)
        pointer         (ptr_coord4,coord4)
c
        integer ntestred,nfailred,ntestreg,nfailreg
c
        common  /testmp/        ntestred,nfailred,ntestreg,nfailreg
c
	common  /xyz_vertex/	ptr_coord,ptr_radius,ptr_coord4
	common  /gmp_info/	scale,eps2,eps3,eps4,eps5
c
	save
c
	data	inf4_1 /2,2,1,1/
	data	sign4_1 /-1,1,1,-1/
c
	data	inf4_2 /0,2,3,3,
     1                  2,0,3,3,
     2                  3,3,0,1,
     3                  3,3,1,0/
c
        data	sign4_2 / 0,1,-1,1,
     1                   -1,0,1,-1,
     2                    1,-1,0,1,
     3                   -1,1,-1,0/
c
	data    sign4_3 /-1,1,-1,1/
c
	data inf5_2	/0,2,1,1,
     1 			2,0,1,1,
     2			1,1,0,1,
     3			1,1,1,0/
c
	data sign5_2    /0,-1,-1,1,
     1			1,0,-1,1,
     2			1,1,0,1,
     3			-1,-1,-1,0/
c
	data inf5_3     /1,1,3,3/
	data sign5_3    /1,1,-1,1/
c
	data order1	/1,2,3,3,1,2,2,3,1/
	data order	/2,3,3,1,1,2/
c
c	To test if the union of the two tetrahedron is convex, we check the
c	position of o with respect to the three faces (a,b,p), (b,c,p)
c	and (c,a,p) of (a,b,c,p). 
c	To do that, we evaluate the three determinants:
c		det(a,b,p,o)
c		det(b,c,p,o)
c		det(c,a,p,o)
c	If the three determinants are positive, and det(a,b,c,p) is negative,
c	then the union is convex
c	Also, if the three determinants are negative, and det(a,b,c,p) is 
c	positive, then the union is convex
c	In all other cases, the union is non convex
c
c	The regularity is tested by computing det(a,b,c,p,o) 
c
c	The computations required are very similar to those used in
c	inside_tetra.f . Look at comments there to decipher this subroutine
c
c	Let us first count how many infinite points we have:
c	(except o)
c	only a and/or b and/or c can be infinite:
c
	list(1) = a
	list(2) = b
	list(3) = c
c
	infpoint(1) = 0
	infpoint(2) = 0
	infpoint(3) = 0
c
	if(a.le.4) infpoint(1) = 1
	if(b.le.4) infpoint(2) = 1
	if(c.le.4) infpoint(3) = 1
c
	ninf = infpoint(1) + infpoint(2) + infpoint(3)
c
	do 100 m = 1,3
		coordp(m) = coord(3*p-3+m)
100	continue
c
c	"General case" : no infinite point
c
	if(ninf.eq.0) then
c
c		First, a simple case: if o is infinite, then
c		det(a,b,c,p,o) = -det(a,b,c,p) and consequently 
c		(a,b,c,p,o) is regular:nothing to do!
c
		if(o.le.4) then
			regular = .true.
			return
		endif
c
c	The three determinants det(a,b,p,o), det(b,c,p,o), and det(c,a,p,o)
c	are "real" 4x4 determinants. 
c	First, we substract the row corresponding to p from the other row,
c	and develop with respect to p. The determinants become:
c
c	det(a,b,p,o)= -| ap(1) ap(2) ap(3) |
c		       | bp(1) bp(2) bp(3) |
c		       | op(1) op(2) op(3) |
c
c	det(b,c,p,o)= -| bp(1) bp(2) bp(3) |
c		       | cp(1) cp(2) cp(3) |
c		       | op(1) op(2) op(3) |
c
c	det(c,a,p,o)= -| cp(1) cp(2) cp(3) |
c		       | ap(1) ap(2) ap(3) |
c		       | op(1) op(2) op(3) |
c
c	where ip(j) = i(j) - p(j) for all i in {a,b,c,o} and j in {1,2,3}
c
c	We compute two types of minors:
c
c		Mbo_ij = bp(i)op(j) - bp(j)op(i)
c	and
c		Mca_ij = cp(i)ap(j) - cp(j)op(i)
c
c	We store Mbo_12 in Mbo(3), Mbo_13 in Mbo(2),...
c
c	Get coordinates
c
		do 200 m = 1,3
			a_p(m) = coord(3*a-3+m) - coordp(m)
			b_p(m) = coord(3*b-3+m) - coordp(m)
			c_p(m) = coord(3*c-3+m) - coordp(m)
			o_p(m) = coord(3*o-3+m) - coordp(m)
200		continue
c
		a_p(4) = coord4(a) - coord4(p)
		b_p(4) = coord4(b) - coord4(p)
		c_p(4) = coord4(c) - coord4(p)
		o_p(4) = coord4(o) - coord4(p)
c
c	Now compute 2x2 determinants Mbo and Mca
c
		Mbo(1) = b_p(2)*o_p(3) - b_p(3)*o_p(2)
		Mbo(2) = b_p(1)*o_p(3) - b_p(3)*o_p(1)
		Mbo(3) = b_p(1)*o_p(2) - b_p(2)*o_p(1)
c
		Mca(1) = c_p(2)*a_p(3) - c_p(3)*a_p(2)
		Mca(2) = c_p(1)*a_p(3) - c_p(3)*a_p(1)
		Mca(3) = c_p(1)*a_p(2) - c_p(2)*a_p(1)
c
c	Now,
c
		det_abpo = - a_p(1)*Mbo(1)+a_p(2)*Mbo(2)
     1				-a_p(3)*Mbo(3)
		det_bcpo = c_p(1)*Mbo(1) - c_p(2)*Mbo(2) 
     1				+ c_p(3)*Mbo(3)
		det_capo = - o_p(1)*Mca(1) + o_p(2)*Mca(2) 
     1				- o_p(3)*Mca(3)
c
c	We also compute:
c
		det_abpc = - b_p(1)*Mca(1) + b_p(2)*Mca(2) 
     1				- b_p(3)*Mca(3)
c
c	Now we compute:
c		det(a,b,c,p,o) = | a(1) a(2) a(3) a(4) 1 |
c				 | b(1) b(2) b(3) b(4) 1 |
c				 | c(1) c(2) c(3) c(4) 1 |
c				 | p(1) p(2) p(3) p(4) 1 |
c				 | o(1) o(2) o(3) o(4) 1 |
c	We first substract row p :
c
c		det(a,b,c,p,o) = - | ap(1) ap(2) ap(3) ap(4) |
c				   | bp(1) bp(2) bp(3) bp(4) |
c				   | cp(1) cp(2) cp(3) cp(4) |
c			 	   | op(1) op(2) op(3) op(4) |
c
c	By developping with respect to the last column, we get:
c
		det_abcpo=-a_p(4)*det_bcpo-b_p(4)*det_capo 
     1			- c_p(4)*det_abpo + o_p(4)*det_abpc
c
c	Test if (a,b,c,p,o) regular, in which case there is no need
c	to flip
c
		if(abs(det_abcpo).lt.eps5) then
                        ntestreg = ntestreg + 1
			call valsort5(a,b,c,p,o,ia,ib,ic,id,ie,iswap)
			call sos_minor5_gmp(coord,radius,scale,ia,ib,
     1			ic,id,ie,val)
			det_abcpo = val*iswap
		endif
c
		if(det_abcpo*itest_abcp.lt.0) then
			regular = .true.
			return
		endif
		regular = .false.
c
c	If not regular, we test for convexity
c
		if(abs(det_abpo).lt.eps4) then
			call valsort4(a,b,p,o,ia,ib,ic,id,iswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			det_abpo = val*iswap
		endif
		if(abs(det_bcpo).lt.eps4) then
			call valsort4(b,c,p,o,ia,ib,ic,id,iswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			det_bcpo = val*iswap
		endif
		if(abs(det_capo).lt.eps4) then
			call valsort4(c,a,p,o,ia,ib,ic,id,iswap)
			call sos_minor4_gmp(coord,scale,ia,ib,ic,id,val)
			det_capo = val*iswap
		endif
c
		test_abpo = det_abpo.gt.0
		test_bcpo = det_bcpo.gt.0
		test_capo = det_capo.gt.0
c
		convex = .false.
		if((itest_abcp*det_abpo).gt.0) return
		if((itest_abcp*det_bcpo).gt.0) return
		if((itest_abcp*det_capo).gt.0) return
		convex = .true.
c
c	Second case: either a,b or c is infinite:
c
	elseif(ninf.eq.1) then
c
c	Let us define as X the infinite point, and (i,j) the pair of finite
c	points.
c	If X = a, (i,j) = (b,c)
c	If X = b, (i,j) = (c,a)
c	If X = c, (i,j) = (a,b)
c	If we define inf(a) = 1 if a infinite, 0 otherwise,
c	then idx_X  = 2 - inf(a) + inf(c)
c
		idx = 2 -infpoint(1) + infpoint(3)
		infp = list(idx)
		i = list(order(1,idx))
		j = list(order(2,idx))
c
c	Get the coordinates
c
		do 300 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			j_p(m) = coord(3*j-3+m) - coordp(m)
			o_p(m) = coord(3*o-3+m) - coordp(m)
300		continue
c
c	First case:	o is finite
c
		if(o.gt.4) then
c
			icol1 = inf4_1(infp)
			sign1 = sign4_1(infp)
c
c	The three 4x4 determinants become:
c
c		-det(i,p,o) [X missing]
c		det(j,p,o) [X missing]
c		det(i,j,p,o)
c
c	And the 5x5 determinant becomes:
c
c		- det(i,j,p,o)
c
			Mjo(1) = j_p(1)*o_p(3) - j_p(3)*o_p(1)
			Mjo(2) = j_p(2)*o_p(3) - j_p(3)*o_p(2)
			Mjo(3) = j_p(1)*o_p(2) - j_p(2)*o_p(1)
c 
c	The correspondence between a,b,c and i,j is not essential
c	We use here the corresponce for a infinite; in the
c	two other cases (b infinite or c infinite), we would
c	have computed the same determinants, but they would
c	not come in the same order
c
			det_abpo = i_p(icol1)*o_p(3)-i_p(3)*o_p(icol1)
			if(abs(det_abpo).lt.eps3) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(coord,scale,ia,ib,
     1				ic,icol1,3,val)
				det_abpo = -val*iswap
			endif
			det_abpo = det_abpo*sign1
			det_capo = - Mjo(icol1)
			if(abs(det_capo).lt.eps3) then
				call valsort3(j,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(coord,scale,ia,ib,
     1				ic,icol1,3,val)
				det_capo = val*iswap
			endif
			det_capo = det_capo*sign1
			det_bcpo =-i_p(1)*Mjo(2)+i_p(2)*Mjo(1)
     1					-i_p(3)*Mjo(3)
			if(abs(det_bcpo).lt.eps3) then
				call valsort4(i,j,p,o,ia,ib,
     1					ic,id,iswap)
				call sos_minor4_gmp(coord,scale,ia,ib,
     1				ic,id,val)
				det_bcpo = val*iswap
			endif
			det_abcpo = -det_bcpo
c
		else
c
c	Second case: o is infinite
c
			info = o
c
c	The three 4x4 determinants become:
c
c		-det(i,p) [o,X missing]
c		det(j,p) [o,X missing]
c		det(i,j,p) [o missing]
c
c	And the 5x5 determinant becomes:
c
c		det(i,j,p) [o,X missing]
c
			icol1 = inf4_2(info,infp)
			sign1 = sign4_2(info,infp)
c
			icol2 = inf4_1(info)
			sign2 = sign4_1(info)
c
			icol5 = inf5_2(info,infp)
			sign5 = sign5_2(info,infp)
c
			det_abpo =-i_p(icol1)*sign1
			if(abs(det_abpo).lt.eps2) then
				call valsort2(i,p,ia,ib,iswap)
				call sos_minor2_gmp(coord,scale,ia,ib,
     1				icol1,val)
				det_abpo = -val*iswap*sign1
			endif
			det_capo =j_p(icol1)*sign1
			if(abs(det_capo).lt.eps2) then
				call valsort2(j,p,ia,ib,iswap)
				call sos_minor2_gmp(coord,scale,ia,ib,
     1				icol1,val)
				det_capo = val*iswap*sign1
			endif
			det_bcpo =i_p(icol2)*j_p(3)
     1				-i_p(3)*j_p(icol2)
			if(abs(det_bcpo).lt.eps3) then
				call valsort3(i,j,p,ia,ib,ic,iswap)
				call sos_minor3_gmp(coord,scale,ia,ib,
     1				ic,icol2,3,val)
				det_bcpo = val*iswap
			endif
			det_bcpo =det_bcpo*sign2
			det_abcpo=i_p(icol5)*j_p(3)
     1				-i_p(3)*j_p(icol5)
			if(abs(det_abcpo).lt.eps3) then
				call valsort3(i,j,p,ia,ib,ic,iswap)
				call sos_minor3_gmp(coord,scale,ia,ib,
     1				ic,icol5,3,val)
				det_abcpo = val*iswap
			endif
			det_abcpo= det_abcpo*sign5
c
		endif
c
c	Test if (a,b,c,p,o) regular, in which case there is no need
c	to flip
c
		if(det_abcpo*itest_abcp.lt.0) then
			regular = .true.
			return
		endif
		regular = .false.
c
c	If not regular, we test for convexity
c
		testc(1) = det_abpo.gt.0
		testc(2) = det_bcpo.gt.0
		testc(3) = det_capo.gt.0
		test_abpo = testc(order1(1,idx))
		test_bcpo = testc(order1(2,idx))
		test_capo = testc(order1(3,idx))
c
		convex = .false.
		if((itest_abcp*det_abpo).gt.0) return
		if((itest_abcp*det_bcpo).gt.0) return
		if((itest_abcp*det_capo).gt.0) return
		convex = .true.
c
c	Now we consider the case where two points are infinite
c
	elseif(ninf.eq.2) then
c
c	Let us define as (k,l) the two infinite points, and i the
c	point that is finite
c	If i = a, (k,l) = (b,c)
c	If i = b, (k,l) = (c,a)
c	If i = c, (k,l) = (a,b)
c
c	Again: i = 2 + inf(a) - inf(c)
c
		idx = 2 + infpoint(1) - infpoint(3)
		i = list(idx)
		k = list(order(1,idx))
		l = list(order(2,idx))
c
c	Get the coordinates
c
		do 400 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			o_p(m) = coord(3*o-3+m) - coordp(m)
400		continue
c
c	First case: o is finite
c
		if(o.gt.4) then
c
c	The three 4x4 determinants become:
c
c		det(i,p,o) [k missing]
c		-det(i,p,o) [l missing]
c		S*det(p,o) [k,l missing, with S =1 if k<l, -1 otherwise]
c
c	The 5x5 determinants become:
c
c		S*det(i,p,o) [k,l missing, with S=1 if k<l, -1 otherwise]
c
			icol1 = inf4_1(k)
			sign1 = sign4_1(k)
			icol2 = inf4_1(l)
			sign2 = sign4_1(l)
			icol4 = inf4_2(k,l)
			sign4 = sign4_2(k,l)
			icol5 = inf5_2(k,l)
			sign5 = sign5_2(k,l)
c
			Mio(1) = i_p(1)*o_p(3) - i_p(3)*o_p(1)
			Mio(2) = i_p(2)*o_p(3) - i_p(3)*o_p(2)
			Mio(3) = i_p(1)*o_p(2) - i_p(2)*o_p(1)
c 
c	The correspondence between a,b,c and i,j,k is not essential
c	We use here the correspondence for a finite; in the
c	two other cases (b finite or c finite), we would
c	have computed the same determinants, but they would
c	not come in the same order
c
			det_abpo = -Mio(icol1)*sign1
			if(abs(det_abpo).lt.eps3) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(coord,scale,ia,ib,
     1					ic,icol1,3,val)
				det_abpo = val*iswap*sign1
			endif
			det_capo =  Mio(icol2)*sign2
			if(abs(det_capo).lt.eps3) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(coord,scale,ia,ib,
     1					ic,icol2,3,val)
				det_capo = -val*iswap*sign2
			endif
			det_bcpo = -o_p(icol4)*sign4
			if(abs(det_bcpo).lt.eps3) then
				call valsort2(p,o,ia,ib,iswap)
				call sos_minor2_gmp(coord,scale,ia,ib,
     1				icol4,val)
				det_bcpo = val*sign4*iswap
			endif
			det_abcpo = - Mio(icol5)*sign5
			if(abs(det_abcpo).lt.eps3) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(coord,scale,ia,ib,
     1				ic,icol5,3,val)
				det_abcpo = val*iswap*sign5
			endif
c
		else
c
c	Second case: o is infinite
c
			info = o
c
c	The three 4x4 determinants become:
c
c		det(i,p) [o,k missing]
c		-det(i,p) [o,l missing]
c		Const [o,k,l missing]
c
c	The 5x5 determinants become:
c
c		Const*det(i,p) [o,k,l missing]
c	
			icol1 = inf4_2(info,k)
			sign1 = sign4_2(info,k)
			icol2 = inf4_2(info,l)
			sign2 = sign4_2(info,l)
c
			call missinf_sign(info,k,l,icol4,iswap)
c
			det_abpo = i_p(icol1)*sign1
			if(abs(det_abpo).lt.eps2) then
				call valsort2(i,p,ia,ib,iswap2)
				call sos_minor2_gmp(coord,scale,
     1				ia,ib,icol1,val)
				det_abpo = val*iswap2*sign1
			endif
			det_capo = -i_p(icol2)*sign2
			if(abs(det_capo).lt.eps2) then
				call valsort2(i,p,ia,ib,iswap2)
				call sos_minor2_gmp(coord,scale,
     1				ia,ib,icol2,val)
				det_capo = -val*iswap2*sign2
			endif
			det_bcpo = sign4_3(icol4)*iswap
			det_abcpo = sign5_3(icol4)*iswap
     1				*i_p(inf5_3(icol4))
			if(abs(det_abcpo).lt.eps2) then
				call valsort2(i,p,ia,ib,iswap2)
				call sos_minor2_gmp(coord,scale,ia,ib,
     1				inf5_3(icol4),val)
				det_abcpo = val*iswap2*iswap*
     1					sign5_3(icol4)
			endif
c
		endif
c
c	Test if (a,b,c,p,o) regular, in which case there is no need
c	to flip
c
		if(det_abcpo*itest_abcp.lt.0) then
			regular = .true.
			return
		endif
		regular = .false.
c
c	If not regular, we test for convexity
c
		testc(1) = det_abpo.gt.0
		testc(2) = det_bcpo.gt.0
		testc(3) = det_capo.gt.0
		test_abpo = testc(order1(1,idx))
		test_bcpo = testc(order1(2,idx))
		test_capo = testc(order1(3,idx))
c
		convex = .false.
		if((itest_abcp*det_abpo).gt.0) return
		if((itest_abcp*det_bcpo).gt.0) return
		if((itest_abcp*det_capo).gt.0) return
		convex = .true.
c
c	We cannot have all three points a,b,c infinite: in this case,
c	the facet a,b,c would be on the convex hull!
c
	elseif(ninf.eq.3) then
		write(6,*) 'This must be an error...'
	endif
c
	return
	end
c
c==========================================================================================
c==========================================================================================
c	Missinf_sign.f	Version1 3/1/2002	Patrice Koehl
c
c	This subroutine takes as input the indices of three "infinite"
c	points (between 1 and 4), finds the index of the "missing"
c	infinite point, and gives the signature of the permutation
c	required to put the three infinite point in order
c
c==========================================================================================
c==========================================================================================
c
	subroutine missinf_sign(i,j,k,l,sign)
c
c==========================================================================================
c	Input:	i,j,k:		the three known infinite points
c
c	Output:	l		the "missing" infinite point
c		sign		the signature of the permutation
c				that orders i,j,k
c==========================================================================================
c
c==========================================================================================
c	Declare all variables
c==========================================================================================
c
	integer	i,j,k,l,sign
	integer	a,b,c,d
c
	save
c
	l = 10 -i -j -k
c
	a = i
	b = j
	c = k
c
	sign = 1
c
	if(a.gt.b) then
		d = a
		a = b
		b = d
		sign = -sign
	endif
c
	if(a.gt.c) then
		d = a
		a = c
		c = d
		sign = -sign
	endif
c
	if(b.gt.c) then
		sign = -sign
	endif
c
	return
	end
c
c==========================================================================================
c==========================================================================================
C	valsort.f	Version 1 3/2/2002	Patrice Koehl
c
c	This file contains a series of routines that sorts integer
c	values in ascending orders, and keep track of the number of
c	flip required (such as to define the signature of the permutation
c	that transforms the un-ordered data into an ordered array)
c
c	In all these routines, the input values are kept unaffected,
c	and new sorted output values are generated
c==========================================================================================
c==========================================================================================
c
c==========================================================================================
c	1. Sort two numbers a and b
c==========================================================================================
c
	subroutine valsort2(a,b,ia,ib,iswap)
c
	integer	a,b,ia,ib,iswap
c
	save
c
	iswap = 1
	if(a.gt.b) then
		ia = b
		ib = a
		iswap = -iswap
	else
		ia = a
		ib = b
	endif
c
	return
	end
c
c==========================================================================================
c	2. Sort three numbers a, b and c
c==========================================================================================
c
	subroutine valsort3(a,b,c,ia,ib,ic,iswap)
c
	integer	a,b,c,ia,ib,ic,iswap,temp
c
	save
c
	call valsort2(a,b,ia,ib,iswap)
c
	ic = c
c
	if(ib.gt.ic) then
		temp = ib
		ib = ic
		ic = temp
		iswap = -iswap
		if(ia.gt.ib) then
			temp = ia
			ia = ib
			ib = temp
			iswap = -iswap
		endif
	endif
c
	return
	end
c
c==========================================================================================
c	3. Sort four numbers a, b, c and d
c==========================================================================================
c
	subroutine valsort4(a,b,c,d,ia,ib,ic,id,iswap)
c
	integer	a,b,c,d,ia,ib,ic,id,iswap,temp
c
	save
c
	call valsort3(a,b,c,ia,ib,ic,iswap)
c
	id = d
c
	if(ic.gt.id) then
		temp = ic
		ic = id
		id = temp
		iswap = -iswap
		if(ib.gt.ic) then
			temp = ib
			ib = ic
			ic = temp
			iswap = -iswap
			if(ia.gt.ib) then
				temp = ia
				ia = ib
				ib = temp
				iswap = -iswap
			endif
		endif
	endif
c
	return
	end
c
c==========================================================================================
c	4. Sort five numbers a, b, c, d and e
c==========================================================================================
c
	subroutine valsort5(a,b,c,d,e,ia,ib,ic,id,ie,iswap)
c
	integer	a,b,c,d,e,ia,ib,ic,id,ie,iswap,temp
c
	save
c
	call valsort4(a,b,c,d,ia,ib,ic,id,iswap)
c
	ie = e
c
	if(id.gt.ie) then
		temp = id
		id = ie
		ie = temp
		iswap = -iswap
		if(ic.gt.id) then
			temp = ic
			ic = id
			id = temp
			iswap = -iswap
			if(ib.gt.ic) then
				temp = ib
				ib = ic
				ic = temp
				iswap = -iswap
				if(ia.gt.ib) then
					temp = ia
					ia = ib
					ib = temp
					iswap = -iswap
				endif
			endif
		endif
	endif
c
	return
	end
c==========================================================================================
c==========================================================================================
c	Flipjw.f	Version 1 12/17/2001	Patrice Koehl
c
c	After a point has been inserted, this subroutine goes over the 
c	link_facet list to restore regularity. When a link_facet is found
c	non_regular and "flippable" (see below), the program attempts
c	to flip it. If the flip is successful, new link_facets are added
c	on the queue.
c	The subroutine ends when the link facet is empty
c
c	This version of flip calls the "jw" flip subroutines, i.e. flip
c	subroutines that do not store a dag
c==========================================================================================
c==========================================================================================
c
	subroutine flipjw(tetra_last)
c
c==========================================================================================
c	Input:
c		- nlink_facet:  when this program is called, the
c				link_facet queue contains four triangles,
c				derived from the tetrahedron in which the
c				new point is added. Each triangle is defined
c				by two tetrahedra, defined by link_facet
c		- link_facet:	the four link facets
c==========================================================================================
c
c	Include array dimensions
c
c==========================================================================================
c	Define variables
c==========================================================================================
c
	integer	j
	integer	p,o,a,b,c
	integer	ierr,ifind,nlink_facet
	integer	itetra,jtetra
	integer	tetra_ab,tetra_ac,tetra_bc
	integer	iorder,ntetra
	integer	ireflex,iflip
	integer	idx_p,idx_o,itest_abcp
	integer	idx_a,idx_b,idx_c
	integer	nfree,nkill,nkill_top,ns
	integer	nfreemax
        integer tetra_last
c
	integer	idxi,idxj,idxk,idxl
	integer	ia,ib,ic,ii,ij
c
	integer*1 ival
c
	integer	facei(3),facej(3),edgei(2),edgej(2),edgek(2)
	integer	edge_val(2,3)
	integer	tetra_flip(3),list_flip(3)
c
	integer table32(3,3),table32_2(2,3),table41(3,3)
	integer	table41_2(2,3)
	integer	vert_flip(5)
c
	integer free(*),kill(*)
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer	link_facet(2,*)
	integer	link_index(2,*)
c
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
c
	logical test,test_or(2,3),regular,convex
	logical test_abpo,test_abpc,test_capo,test_acpb
	logical test_bcpo,test_bcpa,test_acpo
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
	pointer		(ptr_link_facet,link_facet)
	pointer		(ptr_link_index,link_index)
	pointer		(ptr_free,free)
	pointer		(ptr_kill,kill)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common  /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common  /link_zone/	ptr_link_facet,
     1				ptr_link_index,nlink_facet
	common  /freespace/	ptr_free,ptr_kill,nfree,nkill
c
	save
c
	data table32 /1,2,3,1,3,2,3,1,2/
	data table32_2/1,2,1,3,2,3/
	data table41 /2,1,3,1,2,3,1,3,2/
	data table41_2/1,1,2,1,2,2/
c
	nfreemax  = 10000
	nkill_top = nint(nfreemax*0.9)
c
c	Loop over all link facets
c
	j = 0
c
100	if(j.eq.nlink_facet) goto 200
c
	if(nkill.ge.nkill_top) then
c		
		nkill =  nkill_top
c		ns = nfree
c		nfree = min(nfree+nkill,nkill_top)
c		do 150 j = ns+1,nfree
c			free(j) = kill(j-ns)
c150		continue
c		nkill = 0
c
	endif
c
	j = j + 1
c
c	First defined the two tetrahedra that contains the link facet as
c	itetra and jtetra
c
	itetra = link_facet(1,j)
	jtetra = link_facet(2,j)
	idx_p  = link_index(1,j)
	idx_o  = link_index(2,j)
c
c	If the link facet is on the convex hull, discard
c
	if(itetra.eq.0.or.jtetra.eq.0) goto 100
c
c	If these tetrahedra have already been discarded, discard this
c	link facet
c
	if(.not.btest(tetra_info(itetra),1)) then
		if(.not.btest(tetra_info(jtetra),1)) then
			goto 100
		else
			itetra = tetra_neighbour(idx_o,jtetra)
			ival = ibits(tetra_nindex(itetra),2*(idx_o-1),
     1				2)
			idx_p = ival + 1
		endif
	endif
	if(.not.btest(tetra_info(jtetra),1)) then
		jtetra = tetra_neighbour(idx_p,itetra)
		ival = ibits(tetra_nindex(itetra),2*(idx_o-1),2)
		idx_o = ival + 1
	endif
c
c	Let us define the vertices of the two tetrahedra:
c	itetra:		a,b,c,p
c	jtetra:		a,b,c,o
c
	a = tetra(1,itetra)
	b = tetra(2,itetra)
	c = tetra(3,itetra)
	p = tetra(4,itetra)
c
	o = tetra(idx_o,jtetra)
c
	itest_abcp = -1
	if(btest(tetra_info(itetra),0)) itest_abcp = 1
c
c	Check for local regularity (and convexity, at very little
c	extra cost)
c
	call regular_convex(a,b,c,p,o,itest_abcp,regular,convex,
     1	test_abpo,test_bcpo,test_capo)
c
c	if the link facet is locally regular, discard
c
	if(regular) goto 100
c
c	Define neighbors of the facet on itetra and jtetra
c
	call define_facet(itetra,jtetra,idx_o,facei,facej)
c
	test_abpc = itest_abcp.ne.1
c
c	After discarding the trivial case, we now test if the tetrahedra
c	can be flipped. 
c
c	At this stage, I know that the link facet is not locally
c	regular. I still don t know if it is "flippable"
c
c	I first check if {itetra} U {jtetra} is convex. If it is, I
c	perform a 2-3 flip (this is the convexity test performed
c	at the same time as the regularity test)
c
	if(convex) then
		vert_flip(1) = a
		vert_flip(2) = b
		vert_flip(3) = c
		vert_flip(4) = p
		vert_flip(5) = o
		call flipjw_2_3(itetra,jtetra,vert_flip,facei,facej,
     2		test_abpo,test_bcpo,test_capo,ierr,tetra_last)
		goto 100
	endif
c
c	The union of the two tetrahedra is not convex...
c	I now check the edges of the triangle in the link facet, and
c	check if they are "reflexes" (see definition in Edelsbrunner and
c	Shah, Algorithmica (1996), 15:223-241)
c
	ireflex = 0
	iflip = 0
c
c	First check edge (ab): 
c		- (ab) is reflex iff o and c lies on opposite sides of
c		the hyperplane defined by (abp). We therefore test the
c		orientation of (abpo) and (abpc): if they differ (ab)
c		is reflex
c		- if (ab) is reflex, we test if it is of degree 3.
c		(ab) is of degree 3 if it is shared by 3 tetrahedra,
c		namely (abcp), (abco) and (abpo). The first two are itetra
c		and jtetra, so we only need to check if (abpo) exists.
c		since (abpo) contains p, (abp) should then be a link facet
c		of p, so we test all tetrahedra that define link facets
c
c
	if(test_abpo.neqv.test_abpc) then
c
		ireflex = ireflex + 1
c
		call find_tetra(itetra,3,a,b,o,ifind,
     1		tetra_ab,idx_a,idx_b)
c
		if(ifind.eq.1) then
			iflip = iflip + 1
			tetra_flip(iflip) = tetra_ab
			list_flip(iflip) = 1
			edge_val(1,iflip) = idx_a
			edge_val(2,iflip) = idx_b
			test_or(1,iflip) = test_bcpo
			test_or(2,iflip) = .not.test_capo
		endif
c
	endif
c
c	Now check edge (ac): 
c		- (ac) is reflex iff o and b lies on opposite sides of
c		the hyperplane defined by (acp). We therefore test the
c		orientation of (acpo) and (acpb): if they differ (ac)
c		is reflex
c		- if (ac) is reflex, we test if it is of degree 3.
c		(ac) is of degree 3 if it is shared by 3 tetrahedra,
c		namely (abcp), (abco) and (acpo). The first two are itetra
c		and jtetra, so we only need to check if (acpo) exists.
c		since (acpo) contains p, (acp) should then be a link facet
c		of p, so we test all tetrahedra that define link facets
c
	test_acpo = .not.test_capo
	test_acpb = .not.test_abpc
c
	if(test_acpo.neqv.test_acpb) then
c
		ireflex = ireflex + 1
c
		call find_tetra(itetra,2,a,c,o,ifind,
     1		tetra_ac,idx_a,idx_c)
c
		if(ifind.eq.1) then
			iflip = iflip + 1
			tetra_flip(iflip) = tetra_ac
			list_flip(iflip) = 2
			edge_val(1,iflip) = idx_a
			edge_val(2,iflip) = idx_c
			test_or(1,iflip) = .not.test_bcpo
			test_or(2,iflip) = test_abpo
		endif
c
	endif
c
c	Now check edge (bc): 
c		- (bc) is reflex iff o and a lies on opposite sides of
c		the hyperplane defined by (bcp). We therefore test the
c		orientation of (bcpo) and (bcpa): if they differ (bc)
c		is reflex
c		- if (bc) is reflex, we test if it is of degree 3.
c		(bc) is of degree 3 if it is shared by 3 tetrahedra,
c		namely (abcp), (abco) and (bcpo). The first two are itetra
c		and jtetra, so we only need to check if (bcpo) exists.
c		since (bcpo) contains p, (bcp) should then be a link facet
c		of p, so we test all tetrahedra that define link facets
c
	test_bcpa = test_abpc
c
	if(test_bcpo.neqv.test_bcpa) then
c
		ireflex = ireflex + 1
c
		call find_tetra(itetra,1,b,c,o,ifind,
     1		tetra_bc,idx_b,idx_c)
c
		if(ifind.eq.1) then
			iflip = iflip + 1
			tetra_flip(iflip) = tetra_bc
			list_flip(iflip) = 3
			edge_val(1,iflip) = idx_b
			edge_val(2,iflip) = idx_c
			test_or(1,iflip) = test_capo
			test_or(2,iflip) = .not.test_abpo
		endif
c
	endif
c
	if(ireflex.ne.iflip) goto 100
c
	if(iflip.eq.1) then
c
c		Only one edge is "flippable": we do a 3-2 flip
c
		iorder = list_flip(iflip)
		ia = table32(1,iorder)
		ib = table32(2,iorder)
		ic = table32(3,iorder)
		vert_flip(ia) = a
		vert_flip(ib) = b
		vert_flip(ic) = c
		vert_flip(4) = p
		vert_flip(5) = o
		ia = table32_2(1,iorder)
		ib = table32_2(2,iorder)
		edgei(1) = ia
		edgei(2) = ib
		edgej(1) = facej(ia)
		edgej(2) = facej(ib)
		edgek(1) = edge_val(1,iflip)
		edgek(2) = edge_val(2,iflip)
		call flipjw_3_2(itetra,jtetra,tetra_flip(1),vert_flip,
     1		edgei,edgej,edgek,test_or(1,iflip),test_or(2,iflip),
     2          ierr,tetra_last)
c
	elseif(iflip.eq.2) then
c
c		In this case, one point is redundant: the point common to
c		the two edges that can be flipped. We then perform a 4-1
c		flip
c
		iorder = list_flip(1) + list_flip(2) - 2
		vert_flip(table41(1,iorder)) = a
		vert_flip(table41(2,iorder)) = b
		vert_flip(table41(3,iorder)) = c
		vert_flip(4) = p
		vert_flip(5) = o
		ii = table41_2(1,iorder)
		ij = table41_2(2,iorder)
		idxi = iorder
		idxj = facej(iorder)
		idxk = edge_val(ii,1)
		idxl = edge_val(ij,2)
		if(iorder.eq.1) then
			test = test_bcpo
		elseif(iorder.eq.2) then
			test = .not.test_capo
		else
			test = test_abpo
		endif
		call flipjw_4_1(itetra,jtetra,tetra_flip(1),
     1		tetra_flip(2),vert_flip,idxi,idxj,idxk,idxl,test,ierr,
     2          tetra_last)
c
	else
c
c	This case should not occur...
c
		write(6,*) 'Problem...three edges flippable!!'
c
	endif
c
	goto 100
c
200	continue
c
c	Add all "killed" tetrahedra in the free zone
c
	ns = nfree
	nfree = min(ns+nkill,nfreemax)
	do 300 j = ns+1,nfree
		free(j) = kill(j-ns)
300	continue
	nkill = 0
c
	return
	end
c===========================================================================================
c===========================================================================================
c	Flipjw_1_4.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 4->1 flip in 3D for regular triangulation
c
c	a 4->1 flip is a transformation in which a tetrahedron and a single
c	vertex included in the tetrahedron are transformed to 4 tetrahedra,
c	defined from the 4 four faces of the initial tetrahedron, connected
c	to the new point. Each of the faces are then called "link facet",
c	and stored on a queue
c
c	This version of flip_1_4 does not save the old tetrahedron
c	(i.e. no history dag) in order to save space. As a consequence,
c	it cannot be used with a point location scheme that uses the
c	history dag
c
c===========================================================================================
c===========================================================================================
c
	subroutine flipjw_1_4(ipoint,itetra,tetra_last)
c
c===========================================================================================
c	Input:
c		- ipoint:	index of the point p to be included
c		- itetra:	index of the tetrahedra (a,b,c,d) considered
c
c	Output:
c		- nlink_facet:	4
c		- link_facet:	Add the four faces of the initial tetrahedron
c		- link_index:	A link_facet is a triangle defined from its
c				two neighboring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c===========================================================================================
c
c===========================================================================================
c	Define variables
c===========================================================================================
c
	integer	i,j,k,ntetra,newtetra,n_new
	integer	ipoint,itetra,jtetra,nlink_facet
	integer	fact,idx
	integer	nfree,nkill
        integer tetra_last
c
	integer	idx_list(3,4)
c
	integer*1 ival,ikeep
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer	list_new(*)
	integer	link_facet(2,*)
	integer	link_index(2,*)
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
	integer	free(*),kill(*)
c
	integer	vertex(4),neighbour(4),nindex(4)
	integer	position(4)
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
	pointer		(ptr_link_facet,link_facet)
	pointer		(ptr_link_index,link_index)
	pointer		(ptr_free,free)
	pointer		(ptr_kill,kill)
	pointer		(ptr_list_new,list_new)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common  /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common  /link_zone/	ptr_link_facet,
     1				ptr_link_index,nlink_facet
	common  /freespace/	ptr_free,ptr_kill,nfree,nkill
	common  /update/	ptr_list_new,n_new
c
	save
c
	data idx_list /1,1,1,1,2,2,2,2,3,3,3,3/
c
c	Store information about "old" tetrahedron
c
	ikeep = tetra_info(itetra)
c
	do 50 i = 1,4
		vertex(i)    = tetra(i,itetra)
		neighbour(i) = tetra_neighbour(i,itetra)
		ival         = ibits(tetra_nindex(itetra),2*(i-1),2)
		nindex(i)    = ival+1
50	continue
c
	fact = -1
	if(btest(tetra_info(itetra),0)) fact = 1
c
c	The four new tetrahedra are going to be stored
c	in : any free space in the tetrahedron list,
c	and at the end of the list of known tetrahedra
c
	k = 0
c
	do 100 i = nfree,max(nfree-3,1),-1
		k = k + 1
		position(k) = free(i)
100	continue
	nfree = max(nfree-4,0)
c
	do 150 i = k+1,4
		ntetra = ntetra + 1
		position(i) = ntetra
150	continue
        tetra_last = position(4)
c
c	itetra is set to 0, and added to the "kill" list
c
	tetra_info(itetra)= ibclr(tetra_info(itetra),1)
	nkill = 1
	kill(nkill) = itetra
c
c	The tetrahedron is defined as (ijkl); four new tetrahedra are
c	created:	jklp, iklp, ijlp, and ijkp, where p is the new
c	point to be included
c
c	For each new tetrahedron, define all four neighbours:
c	For each neighbour, I store the index of the vertex opposite to 
c	the common face in array tetra_nindex
c
c	tetrahedron jklp : neighbours are iklp, ijlp, ijkp and neighbour
c			   of (ijkl) on face jkl
c	tetrahedron iklp : neighbours are jklp, ijlp, ijkp and neighbour
c			   of (ijkl) on face ikl
c	tetrahedron ijlp : neighbours are jklp, iklp, ijkp and neighbour
c			   of (ijkl) on face ijl
c	tetrahedron ijkp : neighbours are jklp, iklp, ijlp and neighbour
c			   of (ijkl) on face ijk
c
	do 250 i = 1,4
c
		newtetra = position(i)
		n_new = n_new + 1
		list_new(n_new) = newtetra
c
		tetra_info(newtetra) = 0
		tetra_nindex(newtetra) = 0
c
		k = 0
		do 200 j = 1,4
			if(j.eq.i) goto 200
			k = k+1
			tetra(k,newtetra)=vertex(j)
			tetra_neighbour(k,newtetra) = position(j)
			ival = idx_list(k,i) - 1
			call mvbits(ival,0,2,tetra_nindex(newtetra),
     1			2*(k-1))
200		continue
c
		jtetra = neighbour(i)
		idx = nindex(i)
		tetra(4,newtetra) = ipoint
		tetra_neighbour(4,newtetra) = jtetra
		ival = idx-1
		call mvbits(ival,0,2,tetra_nindex(newtetra),6)
c
c		Update surface information
c
		call mvbits(ikeep,2+i,1,tetra_info(newtetra),2+i)
c
c		I must update the neighbors of the neighbour of itetra!
c
		if(jtetra.ne.0.and.idx.ne.0) then
			tetra_neighbour(idx,jtetra) = newtetra
			ival = 3
			call mvbits(ival,0,2,tetra_nindex(jtetra),
     1			2*(idx-1))
		endif
c
		tetra_info(newtetra)=ibset(tetra_info(newtetra),1)
c
c		I also store the orientation of the tetrahedron
c		First, (jklp) and (ijlp) are cw, while (iklp) and (ijkp)
c		are ccw. 
c
		fact = -fact
		if(fact.eq.1) tetra_info(newtetra)=
     1				ibset(tetra_info(newtetra),0)
c
250	continue
c
c	Now add all fours faces of itetra in the link_facet queue.
c	Each link_facet (a triangle) is implicitly defined as the
c	intersection of two tetrahedra
c
c	link_facet:	jkl	tetrahedra:	jklp and neighbour of (ijkl)
c						on jkl
c	link_facet:	ikl	tetrahedra:	iklp and neighbour of (ijkl)
c						on ikl
c	link_facet:	ijl	tetrahedra:	ijlp and neighbour of (ijkl)
c						on ijl
c	link_facet:	ijk	tetrahedra:	ijkp and neighbour of (ijkl)
c						on ijk
c
	nlink_facet = 0
c
	do 300 i = 1,4
		newtetra = position(i)
		nlink_facet = nlink_facet + 1
		link_facet(1,nlink_facet) = newtetra
		link_facet(2,nlink_facet) = tetra_neighbour(4,newtetra)
		link_index(1,nlink_facet) = 4
		ival = ibits(tetra_nindex(newtetra),6,2)
		link_index(2,nlink_facet) = ival+1
300	continue
c
c	end of flip_1_4
c
	return
	end
c========================================================================================
c========================================================================================
c	Flipjw_2_3.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 2->3 flip in 3D for regular triangulation
c
c	a 2->3 flip is a transformation in which two tetrahedrons are
c	flipped into three tetrahedra. The two tetrahedra (abcp) and
c	(abco) shares a triangle (abc) which is in the link_facet of the
c	current point p added to the triangulation. 
c	This flip is only possible if the union of the two tetrahedron is
c	convex, and if their shared triangle is not locally regular. 
c	We assume here that these tests have been performed and are
c	true. 
c	Once the flip has been performed, three new tetrahedra are added
c	and three new "link facet" are added to the link
c	facet queue
c
c	This version of flip_2_3 does not save the old tetrahedra
c	(i.e. no history dag) in order to save space. As a consequence,
c	it cannot be used with a point location scheme that uses a
c	history dag
c========================================================================================
c========================================================================================
c
	subroutine flipjw_2_3(itetra,jtetra,vertices,facei,facej,
     1	test_abpo,test_bcpo,test_capo,ierr,tetra_last)
c
c========================================================================================
c	Input:
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- vertices:	the five vertices a,b,c,o,p
c		- facei		indices of the vertices a,b,c in (a,b,c,p)
c		- facej		indices of the vertices a,b,c in (a,b,c,o)
c		- test_abpo:	orientation of the four points a,b,p,o
c		- test_bcpo:	orientation of the four points b,c,p,o
c		- test_capo:	orientation of the four points c,a,p,o
c
c	Output:
c		- nlink_facet:	3 new link facets are added
c		- link_facet:	the three faces of the initial tetrahedron
c				(a,b,c,o) containing the vertex o are added
c				as link facets
c		- link_index:   A link_facet is a triangle defined from its
c				two neighbouring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c		- ierr:		1 if flip was not possible
c========================================================================================
c
c========================================================================================
c	Define variables
c========================================================================================
c
	integer	i,j,k,p,o
	integer	ierr,ntetra,nlink_facet
	integer	itetra,jtetra
	integer	it,jt,idx,jdx,pos,opos,o_place
	integer	nfree,nkill,newtetra,n_new
        integer tetra_last
c
	integer*1	ival,ikeep,jkeep
c
	integer	jtetra_touch(3),itetra_touch(3)
	integer	jtetra_idx(3),itetra_idx(3)
	integer	idx_list(2,3)
	integer	vertices(5),face(3)
	integer	facei(3),facej(3)
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer	list_new(*)
	integer	link_facet(2,*)
	integer	link_index(2,*)
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
	integer	free(*),kill(*)
c
	integer tests(3),position(3)
c
	logical test_abpo,test_bcpo,test_capo
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
	pointer		(ptr_link_facet,link_facet)
	pointer		(ptr_link_index,link_index)
	pointer		(ptr_free,free)
	pointer		(ptr_kill,kill)
	pointer		(ptr_list_new,list_new)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common  /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common  /link_zone/	ptr_link_facet,
     1				ptr_link_index,nlink_facet
	common  /freespace/	ptr_free,ptr_kill,nfree,nkill
	common  /update/	ptr_list_new,n_new
c
	save
c
	data idx_list/1,1,1,2,2,2/
c
	ierr = 0
c
c	If itetra or jtetra are inactive, cannot flip
c
	if(.not.btest(tetra_info(itetra),1).or..not.
     1		btest(tetra_info(jtetra),1)) then
		ierr = 1
		return
	endif
c
c	Define
c		- itetra_touch: the three tetrahedra that touches itetra on the
c				faces opposite to the 3 vertices a,b,c
c		- itetra_idx:	for the three tetrahedra defined by itetra_touch,
c				index of the vertex opposite to the face
c				common with itetra
c		- jtetra_touch: the three tetrahedra that touches jtetra on the
c				faces opposite to the 3 vertices a,b,c
c		- jtetra_idx:	for the three tetrahedra defined by jtetra_touch,
c				index of the vertex opposite to the face
c				common with jtetra
c
	do 10 i = 1,3
		itetra_touch(i) = tetra_neighbour(facei(i),itetra)
		ival = ibits(tetra_nindex(itetra),2*(facei(i)-1),2)
		itetra_idx(i) = ival + 1
		jtetra_touch(i) = tetra_neighbour(facej(i),jtetra)
		ival = ibits(tetra_nindex(jtetra),2*(facej(i)-1),2)
		jtetra_idx(i) = ival + 1
10	continue
c
c	First three vertices define triangle that is removed
c
	face(1) = vertices(1)
	face(2)	= vertices(2)
	face(3) = vertices(3)
c
	p = vertices(4)
	o = vertices(5)
c
c	The three new tetrahedra are going to be stored
c       in : any free space in the tetrahedron list,
c       and at the end of the list of known tetrahedra if needed
c
	k = 0
	do 50 i = nfree,max(nfree-2,1),-1
		k = k + 1
		position(k) = free(i)
50	continue
	nfree = max(nfree-3,0)
c
	do 100 i = k+1,3
		ntetra = ntetra + 1
		position(i) = ntetra
100	continue
        tetra_last = position(3)
c
c	Set itetra and jtetra to 0, and add them to kill list
c
	ikeep = tetra_info(itetra)
	jkeep = tetra_info(jtetra)
c
	tetra_info(itetra)= ibclr(tetra_info(itetra),1)
	tetra_info(jtetra)= ibclr(tetra_info(jtetra),1)
	kill(nkill+1) = itetra
	kill(nkill+2) = jtetra
	nkill = nkill + 2
c
c	I need :
c		- the vertices a,b,c are the first vertices of itetra
c		- the vertices p and o
c		- for each vertex in the triangle, define the opposing
c		faces in the two tetrahedra itetra and jtetra, and
c		the tetrahedra that share that faces with itetra and
c		jtetra, respectively. This information is stored
c		in two arrays, itetra_touch and jtetra_touch
c
c	These information are provided by the calling program
c
c	For bookkeeping reasons, I always store p as the last vertex
c
c	Now I define the three new tetrahedra: (bcop), (acop) and (abop)
c	as well as their neighbours
c
c	tetrahedron bcop : neighbours are acop, abop, neighbour of (abcp)
c			   on face bcp, and neighbour of (abco) on face bco
c	tetrahedron acop : neighbours are bcop, abop, neighbour of (abcp)
c			   on face acp, and neighbour of (abco) on face aco
c	tetrahedron abop : neighbours are bcop, acop, neighbour of (abcp)
c			   on face abp, and neighbour of (abco) on face abo
c
c
	tests(1) = 1
	if(test_bcpo) tests(1) = -1
	tests(2) = -1
	if(test_capo) tests(2) = 1
	tests(3) = 1
	if(test_abpo) tests(3) = -1
c
	do 200 i = 1,3
c
		newtetra = position(i)
		n_new = n_new + 1
		list_new(n_new) = newtetra
c
		tetra_info(newtetra) = 0
		tetra_nindex(newtetra)  = 0
c
		k = 0
		do 150 j = 1,3
			if(j.eq.i) goto 150
			k = k+1
			tetra(k,newtetra)=face(j)
			tetra_neighbour(k,newtetra) = position(j)
			ival = idx_list(k,i)-1
			call mvbits(ival,0,2,
     1			tetra_nindex(newtetra),2*(k-1))
150		continue
c
		tetra(3,newtetra) = o
		it = itetra_touch(i)
		idx = itetra_idx(i)
		tetra_neighbour(3,newtetra) = it
		ival = idx-1
		call mvbits(ival,0,2,
     1			tetra_nindex(newtetra),4)
		call mvbits(ikeep,2+facei(i),1,
     1			tetra_info(newtetra),5)
		if(idx.ne.0.and.it.ne.0) then
			tetra_neighbour(idx,it) = newtetra
			ival = 2
			call mvbits(ival,0,2,
     1			tetra_nindex(it),2*(idx-1))
		endif
c
		tetra(4,newtetra) = p
		jt = jtetra_touch(i)
		jdx = jtetra_idx(i)
		tetra_neighbour(4,newtetra) = jt
		ival = jdx-1
		call mvbits(ival,0,2,
     1			tetra_nindex(newtetra),6)
		call mvbits(jkeep,2+facej(i),1,
     1			tetra_info(newtetra),6)
		if(jdx.ne.0.and.jt.ne.0) then
			tetra_neighbour(jdx,jt) = newtetra
			ival = 3
			call mvbits(ival,0,2,
     1			tetra_nindex(jt),2*(jdx-1))
		endif
c
		tetra_info(newtetra)=ibset(tetra_info(newtetra),1)
c
		if(tests(i).eq.1) then
			tetra_info(newtetra)=ibset(
     1				tetra_info(newtetra),0)
		endif
c
200	continue
c
c
c	Now add all three faces of jtetra containing o in the link_facet queue.
c	Each link_facet (a triangle) is implicitly defined as the
c	intersection of two tetrahedra
c
c	link_facet:	bco	tetrahedra:	bcop and neighbour of (abco)
c						on bco
c	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
c						on aco
c	link_facet:	abo	tetrahedra:	abop and neighbour of (abco)
c						on abo
c
	do 250 i = 1,3
		newtetra = position(i)
		nlink_facet = nlink_facet + 1
		link_facet(1,nlink_facet) = newtetra
		link_facet(2,nlink_facet) = tetra_neighbour(4,newtetra)
		link_index(1,nlink_facet) = 4
		ival = ibits(tetra_nindex(newtetra),6,2)
		link_index(2,nlink_facet) = ival + 1
250	continue
c
	return
	end
c
c=========================================================================================
c=========================================================================================
c	Flipjw_3_2.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 3->2 flip in 3D for regular triangulation
c
c	a 3->2 flip is a transformation in which three tetrahedrons are
c	flipped into two tetrahedra. The two tetrahedra (abpo), (abcp) and
c	(abco) shares an edge (ab) which is in the link_facet of the
c	current point p added to the triangulation. 
c	This flip is only possible if the edge ab is reflex, with degree 3
c	We assume here that these tests have been performed and are
c	true. 
c	Once the flip has been performed, two new tetrahedra are added
c	and two new "link facet" are added to the link facet queue
c
c	This version of flip_3_2 does not save the old tetrahedron
c	(i.e. no history dag) in order to save space. As a consequence,
c	it cannot be used with a point location scheme that uses the
c	history dag
c=========================================================================================
c=========================================================================================
c
	subroutine flipjw_3_2(itetra,jtetra,ktetra,vertices,
     1			edgei,edgej,edgek,test_bcpo,test_acpo,ierr,
     2                  tetra_last)
c
c=========================================================================================
c	Input:
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- ktetra:	index of the tetrahedra (a,b,o,p) considered
c		- vertices:	the five vertices a,b,c,p,o
c		- edgei		indices of a,b in (a,b,c,p)
c		- edgej		indices of a,b in (a,b,c,o)
c		- edgek		indices of a,b in (a,b,o,p)
c		- test_bcpo   : orientation of the four points b,c,p,o
c		- test_acpo   : orientation of the four points a,c,p,o
c
c	Output:
c		- nlink_facet:	2 new link facets are added
c		- link_facet:	the two faces of the initial tetrahedron
c				(a,b,o,p) containing the edge op are added
c				as link facets
c		- link_index:   A link_facet is a triangle defined from its
c				two neighbouring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c		- ierr:		1 if flip was not possible
c=========================================================================================
c
c=========================================================================================
c	Define variables
c=========================================================================================
c
	integer	i,j,k,p,o,c
	integer	ierr,ntetra,nlink_facet,n_new
	integer	itetra,jtetra,ktetra
	integer	it,jt,kt,idx,jdx,kdx
	integer	nfree,nkill,newtetra,nswap
        integer tetra_last
c
	integer*1 ival
	integer*1 ikeep,jkeep,kkeep
c
	integer	edge(2),tests(2),vertices(5)
	integer	itetra_touch(2),jtetra_touch(2),ktetra_touch(2)
	integer	itetra_idx(2),jtetra_idx(2),ktetra_idx(2)
	integer	position(2)
	integer edgei(2),edgej(2),edgek(2)
c
	integer	free(*),kill(*)
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer	list_new(*)
	integer	link_facet(2,*)
	integer	link_index(2,*)
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
c
	logical test_bcpo,test_acpo
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
	pointer		(ptr_link_facet,link_facet)
	pointer		(ptr_link_index,link_index)
	pointer		(ptr_free,free)
	pointer		(ptr_kill,kill)
	pointer		(ptr_list_new,list_new)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common /link_zone/	ptr_link_facet,
     1				ptr_link_index,nlink_facet
	common /freespace/	ptr_free,ptr_kill,nfree,nkill
	common /update/		ptr_list_new,n_new
c
	save
c
	tests(1) = 1
	if(test_bcpo) tests(1) = -1
	tests(2) = 1
	if(test_acpo) tests(2) = -1
c
	ierr = 0
c
c	If itetra, jtetra or ktetra are inactive, cannot flip
c
	if(.not.btest(tetra_info(itetra),1).or..not.
     1	btest(tetra_info(jtetra),1).or..not.btest(tetra_info(ktetra),1))
     2	then
		ierr = 1
		return
	endif
c
c	Store old info
c
	ikeep = tetra_info(itetra)
	jkeep = tetra_info(jtetra)
	kkeep = tetra_info(ktetra)
c
c	Define
c		- itetra_touch:	indices of the two tetrahedra that share the
c				faces opposite to a and b in itetra,
c				respectively
c		- itetra_idx:   for the two tetrahedra defined by itetra_touch,
c				index position of the vertex opposite to the face
c				common with itetra
c		- jtetra_touch:	indices of the two tetrahedra that share the
c				faces opposite to a and b in jtetra,
c				respectively
c		- jtetra_idx:   for the two tetrahedra defined by jtetra_touch,
c				index position of the vertex opposite to the face
c				common with jtetra
c		- ktetra_touch:	indices of the two tetrahedra that share the
c				faces opposite to a and b in ktetra,
c				respectively
c		- ktetra_idx:   for the two tetrahedra defined by ktetra_touch,
c				index position of the vertex opposite to the face
c				common with ktetra
c
	do 10 i = 1,2
		itetra_touch(i) = tetra_neighbour(edgei(i),itetra)
		jtetra_touch(i) = tetra_neighbour(edgej(i),jtetra)
		ktetra_touch(i) = tetra_neighbour(edgek(i),ktetra)
		ival = ibits(tetra_nindex(itetra),2*(edgei(i)-1),2)
		itetra_idx(i) = ival + 1
		ival = ibits(tetra_nindex(jtetra),2*(edgej(i)-1),2)
		jtetra_idx(i) = ival + 1
		ival = ibits(tetra_nindex(ktetra),2*(edgek(i)-1),2)
		ktetra_idx(i) = ival + 1
10	continue
c
c
	edge(1) = vertices(1)
	edge(2) = vertices(2)
	c       = vertices(3)
	p       = vertices(4)
	o       = vertices(5)
c
c	The two new tetrahedra are going to be stored "free" space, or 
c	at the end of the list
c
	k = 0
	do 50 i = nfree,max(nfree-1,1),-1
		k = k + 1
		position(k) = free(i)
50	continue
	nfree = max(nfree-2,0)
c
	do 100 i = k+1,2
		ntetra = ntetra + 1
		position(i) = ntetra
100	continue
        tetra_last = position(2)
c
c	itetra, jtetra and ktetra becomes "available"; they are added to the
c	"kill" list
c
	tetra_info(itetra)=ibclr(tetra_info(itetra),1)
	tetra_info(jtetra)=ibclr(tetra_info(jtetra),1)
	tetra_info(ktetra)=ibclr(tetra_info(ktetra),1)
c
	kill(nkill+1) = itetra
	kill(nkill+2) = jtetra
	kill(nkill+3) = ktetra
	nkill = nkill + 3
c
c	I need :
c		- the two vertices that define their common edge (ab)
c		 these vertices are stored in the array edge
c		- the vertices c, p and o that form the new triangle
c		- for each vertex in the edge (ab), define the opposing
c		faces in the three tetrahedra itetra, jtetra and ktetra, and
c		the tetrahedron that share these faces with itetra, jtetra and
c		ktetra, respectively. This information is stored
c		in three arrays, itetra_touch, jtetra_touch and ktetra_touch
c
c	These information are given by the calling program
c
c	For bookkeeping reasons, I always set p to be the last vertex
c	of the new tetrahedra
c
c	Now I define the two new tetrahedra: (bcop) and (acop)
c	as well as their neighbours
c
c	tetrahedron bcop : neighbours are acop, neighbour of (abop)
c			   on face bpo, neighbour of (abcp) on face bcp
c			   and neighbour of (abco) on face (bco)
c	tetrahedron acop : neighbours are bcop, neighbour of (abop)
c			   on face apo, neighbour of (abcp) on face acp
c			   and neighbour of (abco) on face (aco)
c
	do 200 i = 1,2
c
		newtetra = position(i)
		n_new = n_new + 1
		list_new(n_new) = newtetra
c
		tetra_info(newtetra)   = 0
		tetra_nindex(newtetra) = 0
c
		k = 0
		do 150 j = 1,2
			if(j.eq.i) goto 150
			k = k+1
			tetra(k,newtetra)=edge(j)
			tetra_neighbour(k,newtetra) = position(j)
150		continue
c
		tetra(2,newtetra) = c
		kt = ktetra_touch(i)
		kdx = ktetra_idx(i)
		tetra_neighbour(2,newtetra) = kt
		ival = kdx - 1
		call mvbits(ival,0,2,tetra_nindex(newtetra),2)
		call mvbits(kkeep,2+edgek(i),1,
     1		tetra_info(newtetra),4)
		if(kdx.ne.0.and.kt.ne.0) then
			tetra_neighbour(kdx,kt) = newtetra
			ival = 1
			call mvbits(ival,0,2,tetra_nindex(kt),2*(kdx-1))
		endif
c
		tetra(3,newtetra) = o
		it = itetra_touch(i)
		idx = itetra_idx(i)
		tetra_neighbour(3,newtetra) = it
		ival = idx - 1
		call mvbits(ival,0,2,tetra_nindex(newtetra),4)
		call mvbits(ikeep,2+edgei(i),1,
     1		tetra_info(newtetra),5)
		if(idx.ne.0.and.it.ne.0) then
			tetra_neighbour(idx,it) = newtetra
			ival = 2
			call mvbits(ival,0,2,tetra_nindex(it),2*(idx-1))
		endif
c
		tetra(4,newtetra) = p
		jt = jtetra_touch(i)
		jdx = jtetra_idx(i)
		tetra_neighbour(4,newtetra) = jt
		ival = jdx - 1
		call mvbits(ival,0,2,tetra_nindex(newtetra),6)
		call mvbits(jkeep,2+edgej(i),1,
     1		tetra_info(newtetra),6)
		if(jdx.ne.0.and.jt.ne.0) then
			tetra_neighbour(jdx,jt) = newtetra
			ival = 3
			call mvbits(ival,0,2,tetra_nindex(jt),2*(jdx-1))
		endif
c
		tetra_info(newtetra) = ibset(tetra_info(newtetra),1)
c
		if(tests(i).eq.1) then
			tetra_info(newtetra) = 
     1				ibset(tetra_info(newtetra),0)
		endif
c
200	continue
c
c	Now add the two faces of ktetra containing (co) in the link_facet 
c	queue.
c	Each link_facet (a triangle) is implicitly defined as the
c	intersection of two tetrahedra
c
c	link_facet:	bco	tetrahedra:	bcop and neighbour of (abco)
c						on bco
c	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
c						on aco
c
	do 700 i = 1,2
		newtetra = position(i)
		nlink_facet = nlink_facet + 1
		link_facet(1,nlink_facet) = newtetra
		link_facet(2,nlink_facet) = tetra_neighbour(4,newtetra)
		link_index(1,nlink_facet) = 4
		ival = ibits(tetra_nindex(newtetra),6,2)+1
		link_index(2,nlink_facet) = ival
700	continue
c
	return
	end
c
c==================================================================================================
c==================================================================================================
c	Flipjw_4_1.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 4->1 flip in 3D for regular triangulation
c
c	a 4->1 flip is a transformation in which four tetrahedra are
c	flipped into one tetrahedron. The two tetrahedra (abop), (bcop),
c	(abcp) and (abco) shares a vertex (b) which is in the link_facet of the
c	current point p added to the triangulation. After the flip, b
c	is set to redundant. 
c	This flip is only possible if the two edges (ab) and (bc)
c	are reflex of order 3.
c	We assume here that these tests have been performed and are
c	true. 
c	Once the flip has been performed, one tetrahedron is added
c	and one new "link facet" is added to the link facet queue
c
c       This version of flip_4_1 does not save the old tetrahedron
c       (i.e. no history dag) in order to save space. As a consequence,
c       it cannot be used with a point location scheme that uses the
c       history dag
c==================================================================================================
c==================================================================================================
c
	subroutine flipjw_4_1(itetra,jtetra,ktetra,ltetra,vertices,
     1			idp,jdp,kdp,ldp,test_acpo,ierr,tetra_last)
c
c==================================================================================================
c	Input:
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- ktetra:	index of the tetrahedra (a,b,o,p) considered
c		- ltetra:	index of the tetrahedra (b,c,o,p) considered
c		- vertices:	index of a,b,c,p,o
c		- idp		index of b in (a,b,c,p)
c		- jdp		index of b in (a,b,c,o)
c		- kdp		index of b in (a,b,o,p)
c		- ldp		index of b in (b,c,o,p)
c		- test_acpo:	orientation of the 4 points (a,c,p,o)
c
c	Output:
c		- nlink_facet:	1 new link facet is added
c		- link_facet:	the face of the initial tetrahedron
c				(a,b,c,o) opposite to the vertex b is added
c				as link facet
c		- link_index:   A link_facet is a triangle defined from its
c				two neighbouring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c		- ierr:		1 if flip was not possible
c
c==================================================================================================
c
c==================================================================================================
c	Define variables
c==================================================================================================
c
	integer	i,j
	integer	p,o,a,b,c
	integer	ierr,npoints,nvertex,nlink_facet,n_new
	integer	ntetra,itetra,jtetra,ktetra,ltetra
	integer	ishare,jshare,kshare,lshare
	integer	idx,jdx,kdx,ldx
	integer	idp,jdp,kdp,ldp
	integer	test1,newtetra
	integer	nfree,nkill,nswap
        integer tetra_last
c
	integer*1 ival
	integer*1 ikeep,jkeep,kkeep,lkeep
	integer	vertices(5)
c
	integer*1 vertex_info(*)
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer	list_new(*)
	integer	link_facet(2,*)
	integer	link_index(2,*)
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
	integer	free(*),kill(*)
c
	logical test_acpo
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
	pointer		(ptr_link_facet,link_facet)
	pointer		(ptr_link_index,link_index)
	pointer		(ptr_free,free)
	pointer		(ptr_kill,kill)
	pointer		(ptr_list_new,list_new)
	pointer		(ptr_vertex_info,vertex_info)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common	/tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common /vertex_zone/	ptr_vertex_info,npoints,nvertex
	common /link_zone/	ptr_link_facet,
     1				ptr_link_index,nlink_facet
	common /freespace/	ptr_free,ptr_kill,nfree,nkill
	common /update/		ptr_list_new,n_new
c
	save
c
	ierr = 0
c
	test1 = 1
	if(test_acpo) test1 = -1
c
c	If itetra, jtetra, ktetra, ltetra are inactive, cannot flip
c
	if(.not.btest(tetra_info(itetra),1).or..not.
     1	btest(tetra_info(jtetra),1).or..not.btest(tetra_info(ktetra),1)
     2  .or..not.btest(tetra_info(ltetra),1)) then
		ierr = 1
		return
	endif
c
c	Store "old" info
c
	ikeep = tetra_info(itetra)
	jkeep = tetra_info(jtetra)
	kkeep = tetra_info(ktetra)
	lkeep = tetra_info(ltetra)
c	
c	Define
c		- ishare:	index of tetrahedron sharing the face 
c				opposite to b in itetra
c		- idx		index of the vertex of ishare opposite to the
c				face of ishare shared with itetra
c		- jshare:	index of tetrahedron sharing the face 
c				opposite to b in jtetra
c		- jdx		index of the vertex of jshare opposite to the
c				face of jshare shared with jtetra
c		- kshare:	index of tetrahedron sharing the face 
c				opposite to b in ktetra
c		- kdx		index of the vertex of kshare opposite to the
c				face of kshare shared with ktetra
c		- lshare:	index of tetrahedron sharing the face 
c				opposite to b in ltetra
c		- ldx		index of the vertex of lshare opposite to the
c				face of lshare shared with ltetra
c
	ishare = tetra_neighbour(idp,itetra)
	jshare = tetra_neighbour(jdp,jtetra)
	kshare = tetra_neighbour(kdp,ktetra)
	lshare = tetra_neighbour(ldp,ltetra)
c
	ival = ibits(tetra_nindex(itetra),2*(idp-1),2)
	idx  = ival + 1
	ival = ibits(tetra_nindex(jtetra),2*(jdp-1),2)
	jdx  = ival + 1
	ival = ibits(tetra_nindex(ktetra),2*(kdp-1),2)
	kdx  = ival + 1
	ival = ibits(tetra_nindex(ltetra),2*(ldp-1),2)
	ldx  = ival + 1
c
c	The new tetrahedron is going to be store in place of itetra
c
	if(nfree.ne.0) then
		newtetra = free(nfree)
		nfree = nfree -1
	else
		ntetra = ntetra + 1
		newtetra = ntetra
	endif
        tetra_last = newtetra
c
	n_new = n_new + 1
	list_new(n_new) = newtetra
c
	tetra_info(newtetra) = 0
	tetra_nindex(newtetra) = 0
c
c	jtetra, ktetra and ltetra become "available"; they
c	are added to the "kill" zone
c
	kill(nkill+1) = itetra
	kill(nkill+2) = jtetra
	kill(nkill+3) = ktetra
	kill(nkill+4) = ltetra
c
	nkill = nkill + 4
c
	tetra_info(itetra) = ibclr(tetra_info(itetra),1)
	tetra_info(jtetra) = ibclr(tetra_info(jtetra),1)
	tetra_info(ktetra) = ibclr(tetra_info(ktetra),1)
	tetra_info(ltetra) = ibclr(tetra_info(ltetra),1)
c
c	I need :
c		- the vertex b that is shared by all 4 tetrahedra
c		- the vertices a, c, p and o
c		- for each tetrahedron, find neighbour attached to the face
c		oposite to b; this information is stored if *share,
c		where * can be i, j, k or l
c
c	These information are provided by the calling program
c
	a = vertices(1)
	b = vertices(2)
	c = vertices(3)
	p = vertices(4)
	o = vertices(5)
c
c	For bookkeeping reason, p is set to be the last vertex of the
c	new tetrahedron
c
c	Now I define the new tetrahedron: (acop)
c
c	tetrahedron acop : neighbor of (bcop) on face cpo, neighbor of (abop)
c			   on face apo, neighbor of (abcp) on face acp
c			   and neighbor of (abco) on face aco
c
	vertex_info(b) = ibclr(vertex_info(b),0)
c
	tetra(1,newtetra) = a
	tetra_neighbour(1,newtetra) = lshare
	ival = ldx - 1
	call mvbits(ival,0,2,tetra_nindex(newtetra),0)
	call mvbits(lkeep,2+ldp,1,tetra_info(newtetra),3)
	if(lshare.ne.0.and.ldx.ne.0) then
		tetra_neighbour(ldx,lshare) = newtetra
		ival = 0
		call mvbits(ival,0,2,tetra_nindex(lshare),2*(ldx-1))
	endif
c
	tetra(2,newtetra) = c
	tetra_neighbour(2,newtetra) = kshare
	ival = kdx - 1
	call mvbits(ival,0,2,tetra_nindex(newtetra),2)
	call mvbits(kkeep,2+kdp,1,tetra_info(newtetra),4)
	if(kshare.ne.0.and.kdx.ne.0) then
		tetra_neighbour(kdx,kshare) = newtetra
		ival = 1
		call mvbits(ival,0,2,tetra_nindex(kshare),2*(kdx-1))
	endif
c
	tetra(3,newtetra) = o
	tetra_neighbour(3,newtetra) = ishare
	ival = idx - 1
	call mvbits(ival,0,2,tetra_nindex(newtetra),4)
	call mvbits(ikeep,2+idp,1,tetra_info(newtetra),5)
	if(ishare.ne.0.and.idx.ne.0) then
		tetra_neighbour(idx,ishare) = newtetra
		ival = 2
		call mvbits(ival,0,2,tetra_nindex(ishare),2*(idx-1))
	endif
c
	tetra(4,newtetra) = p
	tetra_neighbour(4,newtetra) = jshare
	ival = jdx - 1
	call mvbits(ival,0,2,tetra_nindex(newtetra),6)
	call mvbits(jkeep,2+jdp,1,tetra_info(newtetra),6)
	if(jshare.ne.0.and.jdx.ne.0) then
		tetra_neighbour(jdx,jshare) = newtetra
		ival = 3
		call mvbits(ival,0,2,tetra_nindex(jshare),2*(jdx-1))
	endif
c
	tetra_info(newtetra)=ibset(tetra_info(newtetra),1)
c
	if(test1.eq.1) then
		tetra_info(newtetra)=ibset(tetra_info(newtetra),0)
	endif
c
c	Now add one link facet : 
c
c	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
c						on aco
c
	nlink_facet = nlink_facet + 1
	link_facet(1,nlink_facet) = newtetra
	link_facet(2,nlink_facet) = jshare
	link_index(1,nlink_facet) = 4
	link_index(2,nlink_facet) = jdx
c
c	end of flip
c
	return
	end
c
c=================================================================================================
c=================================================================================================
c	Define_facet.f		Version 1 12/21/2001	Patrice Koehl
c
c	A triangle (or facet) is defined by the intersection of two 
c	tetrahedra itetra and jtetra
c	If we know the position of its three vertices in the first
c	tetrahedron (in fact the first three vertices a,b and c
c	of itetra), we need to find the indices of these vertices
c	in the second tetrahedron.
c	This routine also stores information about the neighbours
c	of the two tetrahedra considered
c
c	The vertices are called a,b,c,p, and o, where (abc) is the
c	common facet
c=================================================================================================
c=================================================================================================
c
	subroutine define_facet(itetra,jtetra,idx_o,facei,facej)
c
c=================================================================================================
c	Input:
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- idx_o:	position of o in the vertices of jtetra
c
c	Output:
c		- itouch	itouch(i) is the tetrahedron sharing
c				the face opposite to i in tetrahedron itetra
c		- idx		idx(i) is the vertex of itouch(i) opposite
c				to the face shared with itetra
c		- jtouch	jtouch(i) is the tetrahedron sharing
c				the face opposite to i in tetrahedron jtetra
c		- jdx		jdx(i) is the vertex of jtouch(i) opposite
c				to the face shared with jtetra
c=================================================================================================
c
c=================================================================================================
c	Define variables
c=================================================================================================
c
	integer	i,k
	integer	ia,ib,ie,if
	integer	ntetra,itetra,jtetra
	integer	idx_o
c
	integer*1	ival
c
	integer	other(3,4),other2(2,4,4)
	integer facei(3),facej(3)
c
	integer*1 tetra_info(*),tetra_nindex(*)
c
	integer	tetra(4,*)
	integer	tetra_neighbour(4,*)
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common  /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
c
	save
c
	data other/2,3,4,
     1		   1,3,4,
     2		   1,2,4,
     3		   1,2,3/
c
	data other2/0,0,3,4,2,4,2,3,3,4,0,0,1,4,1,3,2,4,1,4,0,0,1,2,
     1			2,3,1,3,1,2,0,0/
c
c	I need to :
c		- find the three vertices that define their common face
c		 these vertices are stored in the array triangle
c		- find the vertices p and o
c
c	To define the common face of the two tetrahedra itetra and jtetra,
c	I look at the neighbours of itetra : one of them is jtetra!
c	This also provides p. The same procedure is repeated for jtetra,
c	to get o
c
	do 100 i = 1,3
		facei(i) = i
100	continue
c
	ia = tetra(1,itetra)
	do 200 i = 1,3
		k = other(i,idx_o)
		ie = tetra(k,jtetra)
		if(ia.eq.ie) then
			facej(1) = k
			goto 300
		endif
200	continue
c
300	continue
c
	ib = tetra(2,itetra)
	ie = other2(1,facej(1),idx_o)
	if = other2(2,facej(1),idx_o)
	if(ib.eq.tetra(ie,jtetra)) then
		facej(2) = ie
		facej(3) = if
	else
		facej(2) = if
		facej(3) = ie
	endif
c
	return
	end
c
c============================================================================================
c============================================================================================
c	Find_tetra.f		Version 1 1/8/2002	Patrice Koehl
c
c	This subroutine tests if four given points form an existing
c	tetrahedron in the current Delaunay
c============================================================================================
c============================================================================================
c
	subroutine find_tetra(itetra,idx_c,a,b,o,ifind,tetra_loc,
     1		idx_a,idx_b)
c
c============================================================================================
c	Input:
c		itetra:		index of tetrahedra (abcp)
c		idx_c:		index of c in (abcp)
c		o:		index of vertex o
c
c	Output:
c		ifind:		1 if tetrahedron exists, 0 otherwise
c		tetra_loc:	index of existing tetrahedron, if it exists
c============================================================================================
c
c	We are testing if tetrahedron (abpo) exists. If it exists, it is
c	a neighbour of abcp, on the face opposite to vertex c.
c	We test that tetrahedron and see if it contains o
c
	integer	i,ifind,itetra,tetra_loc
	integer	ot,otx,otest
	integer	idx_c,idx_a,idx_b,o,a,b
	integer	ntetra
c
	integer*1 ival
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer tetra(4,*),tetra_neighbour(4,*)
c
	pointer	(ptr_tetra,tetra)
	pointer (ptr_tetra_neighbour,tetra_neighbour)
	pointer (ptr_tetra_info,tetra_info)
	pointer (ptr_tetra_nindex,tetra_nindex)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common  /tetra_stat/    ptr_tetra_info,ptr_tetra_nindex
c
	save
c
	ot = tetra_neighbour(idx_c,itetra)
	ival = ibits(tetra_nindex(itetra),2*(idx_c-1),2)
	otx = ival + 1
c
	otest = tetra(otx,ot)
c
	if(otest.eq.o) then
		ifind = 1
		tetra_loc = ot
c
c		We found the tetrahedron, let us define the position
c		of a and b in this tetrahedron
c
		do 100 i = 1,4
			if(tetra(i,tetra_loc).eq.a) then
				idx_a = i
			elseif(tetra(i,tetra_loc).eq.b) then
				idx_b = i
			endif
100		continue
c		
	else
		ifind = 0
	endif
c
	return
	end
c
c============================================================================================
c============================================================================================
c	Remove_inf.f		Version 1 1/16/2002	Patrice Koehl
c
c	This subroutine sets to 0 the status of tetrahedron that
c	contains infinite points
c============================================================================================
c============================================================================================
c
	subroutine remove_inf
c
	integer	i,a,b,c,d
	integer	ntetra,ninf
	integer	npoint,nvertex
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer*1 vertex_info(*)
	integer tetra(4,*)
	integer tetra_neighbour(4,*)
c
	pointer	(ptr_tetra,tetra)
	pointer (ptr_tetra_neighbour,tetra_neighbour)
	pointer (ptr_tetra_info,tetra_info)
	pointer (ptr_tetra_nindex,tetra_nindex)
	pointer (ptr_vertex_info,vertex_info)
c
	common	/tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common	/tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common  /vertex_zone/	ptr_vertex_info,npoint,nvertex
c
	save
c
	do 100 i = 1,ntetra
c
		if(.not.btest(tetra_info(i),1)) goto 100
c
		a = tetra(1,i)
		b = tetra(2,i)
		c = tetra(3,i)
		d = tetra(4,i)
c
		if(a.le.4.or.b.le.4.or.c.le.4.or.d.le.4) then
			tetra_info(i)=ibset(tetra_info(i),2)
			tetra_info(i)=ibclr(tetra_info(i),1)
			if(a.le.4) call mark_zero(i,1)
			if(b.le.4) call mark_zero(i,2)
			if(c.le.4) call mark_zero(i,3)
			if(d.le.4) call mark_zero(i,4)
		endif
c
100	continue
c
	do 200 i = 1,4
		vertex_info(i)=ibclr(vertex_info(i),0)
200	continue
c
	return
	end
c
c============================================================================================
c============================================================================================
c	Mark_zero.f	Version 1 3/11/2002	Patrice Koehl
c
c	This subroutine marks the tetrahedron that touches
c	a tetrahedron with infinite point as part of the
c	convex hull (i.e. one of its neighbor is 0)
c============================================================================================
c============================================================================================
c
	subroutine mark_zero(itetra,ivertex)
c
	integer	ntetra
	integer	itetra,ivertex
	integer	jtetra,jvertex
c
	integer*1 ival
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer tetra(4,*)
	integer tetra_neighbour(4,*)
c
	pointer	(ptr_tetra,tetra)
	pointer (ptr_tetra_neighbour,tetra_neighbour)
	pointer (ptr_tetra_info,tetra_info)
	pointer (ptr_tetra_nindex,tetra_nindex)
c
	common	/tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common	/tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
c
	save
c
	jtetra = tetra_neighbour(ivertex,itetra)
c
	if(jtetra.ne.0) then
		ival = ibits(tetra_nindex(itetra),2*(ivertex-1),2)
		jvertex = ival + 1
		tetra_neighbour(jvertex,jtetra) = 0
	endif
c
	return
	end
c============================================================================================
c============================================================================================
c	Isort_indx.f		Version 1 1/10/2002	Patrice Koehl
c
c	This subroutine sorts 4 integers, gives the number of swaps
c	required for sorting, and provides the rank of the initial
c	numbers in the sorted array
c============================================================================================
c============================================================================================
c
	subroutine  isort_indx(list,idx,nswap,n)
c
	integer	i,j,a,n,nswap
	integer	list(n),idx(n)
c
	save
c
	do 100 i = 1,n
		idx(i) = i
100	continue
c
	nswap = 0
c
	do 300 i = 1,n-1
		do 200 j = i+1,n
			if(list(i).gt.list(j)) then
				a = list(i)
				list(i) = list(j)
				list(j) = a
				a = idx(i)
				idx(i) = idx(j)
				idx(j) = a
				nswap = nswap + 1
			endif
200		continue
300	continue
c
	return
	end
c
c============================================================================================
c============================================================================================
c	Isort_swap.f		Version 1 1/10/2002	Patrice Koehl
c
c	This subroutine sorts n integers, and gives the number of swaps
c	required for sorting
c============================================================================================
c============================================================================================
c
	subroutine  isort_swap(list,nswap,n)
c
	integer	i,j,a,n,nswap
c
	integer	list(n)
c
	save
c
	nswap = 0
c
	do 300 i = 1,n-1
		do 200 j = i+1,n
			if(list(i).gt.list(j)) then
				a = list(i)
				list(i) = list(j)
				list(j) = a
				nswap = nswap + 1
			endif
200		continue
300	continue
c
	return
	end
c
c============================================================================================
c============================================================================================
c	isort4_swap.f	Version 1 2/7/2002	Patrice Koehl
c
c	This subroutine sorts 4 numbers, and gives the number of swaps
c	required for sorting
c============================================================================================
c============================================================================================
c
	subroutine isort4_swap(a,b,c,d,nswap)
c
	integer	i,a,b,c,d,nswap
c
	save
c
	nswap = 0
c
	if(a.gt.b) then
		i = a
		a = b
		b = i
		nswap = 1
	endif
c
	if(a.gt.c) then
		i = a
		a = c
		c = i
		nswap = nswap + 1
	endif
c
	if(a.gt.d) then
		i = a
		a = d
		d = i
		nswap = nswap + 1
	endif
c
	if(b.gt.c) then
		i = b
		b = c
		c = i
		nswap = nswap + 1
	endif
c
	if(b.gt.d) then
		i = b
		b = d
		d = i
		nswap = nswap + 1
	endif
c
	if(c.gt.d) then
		i = c
		c = d
		d = i
		nswap = nswap + 1
	endif
c
	return
	end
c============================================================================================
c============================================================================================
c	Peel.f	Version 1 4/1/2002	Patrice Koehl
c
c	This subroutine removes the flat tetrahedra at the boundary 
c	of the DT
c============================================================================================
c============================================================================================
c
	subroutine peel
c
	integer	i,j,k,l
	integer	ia,ib,ic,id,val
	integer	ntetra
c
	integer*1 ival
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer tetra(4,*)
	integer tetra_neighbour(4,*)
c
	real*8  eps2,eps3,eps4,eps5,scale
	real*8	vol
c
	real*8	coord(*),radius(*)
	real*8	coord4(*)
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
        pointer         (ptr_coord,coord)
        pointer         (ptr_radius,radius)
        pointer         (ptr_coord4,coord4)
c
	common  /xyz_vertex/	ptr_coord,ptr_radius,ptr_coord4
	common	/tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
	common	/tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
	common  /gmp_info/	scale,eps2,eps3,eps4,eps5
c
	save
c
c	Loop over all tetrahedra, and test the tetrahedra at
c	the boundary
c
	do 400 i = 1,ntetra
c
		if(.not.btest(tetra_info(i),1)) goto 400
c
		do 200 j = 1,4
			if(tetra_neighbour(j,i).eq.0) goto 300
200		continue
c
c		If we get here, the tetrahedron idx is interior, and
c		cannot be flat (see Edelsbrunner,...)
c
		goto 400
c
300		continue
c
c		This is a tetrahedron at the boundary: we test
c		if it is flat, i.e. if its volume is 0
c
		ia = tetra(1,i)
		ib = tetra(2,i)
		ic = tetra(3,i)
		id = tetra(4,i)
c
		call tetra_vol(coord,ia,ib,ic,id,vol)
c
		if(abs(vol).lt.eps4) then
			call minor4_gmp(coord,scale,ia,ib,ic,id,val)
			if(ival.eq.0) then
				tetra_info(i)=ibset(tetra_info(i),2)
			endif
		endif
c
400	continue
c
c	Now we remove that flat tetrahedra, and update the links
c	to their neighbours
c
	do 600 i = 1,ntetra
c
		if(btest(tetra_info(i),2)) then
		    if(btest(tetra_info(i),1)) then
			tetra_info(i)=ibclr(tetra_info(i),1)
			do 500 j = 1,4
				k = tetra_neighbour(j,i)
				if(k.ne.0) then
					ival=ibits(tetra_nindex(j),
     1					2*(i-1),2)
					l = ival + 1
					tetra_neighbour(l,k)=0
				endif
500			continue
c
		   endif
		endif
c
600	continue
c
	return
	end
c
c============================================================================================
c============================================================================================
c	tetra_vol.f	Version 1 4/1/2002	Patrice Koehl
c
c	This subroutine computes the volume of a tetrahedron
c============================================================================================
c============================================================================================
c
	subroutine tetra_vol(coord,ia,ib,ic,id,vol)
c
c	Input:
c		- coord:	array containing coordinates of all vertices
c		- ia,ib,ic,id	four vertices defining the tetrahedron
c
c	Output:
c		- vol		volume of the tetrahedron (floating
c				point calculation)
c
	integer	i
	integer	ia,ib,ic,id
c
	save
c
	real*8	vol
	real*8	ad(3),bd(3),cd(3)
	real*8	Sbcd(3)
	real*8	coord(*)
c
c	The volume of the tetrahedron is proportional to:
c
c	vol = det | a(1)  a(2)  a(3)  1|
c		  | b(1)  b(2)  b(3)  1|
c		  | c(1)  c(2)  c(3)  1|
c		  | d(1)  d(2)  d(3)  1|
c
c	After substracting the last row from the first 3 rows, and
c	developping with respect to the last column, we obtain:
c
c	vol = det | ad(1)  ad(2)  ad(3) |
c		  | bd(1)  bd(2)  bd(3) |
c		  | cd(1)  cd(2)  cd(3) |
c
c	where ad(i) = a(i) - d(i), ...
c
	do 100 i = 1,3
		ad(i) = coord(3*(ia-1)+i) - coord(3*(id-1)+i)
		bd(i) = coord(3*(ib-1)+i) - coord(3*(id-1)+i)
		cd(i) = coord(3*(ic-1)+i) - coord(3*(id-1)+i)
100	continue
c
	Sbcd(3) = bd(1)*cd(2) - cd(1)*bd(2)
	Sbcd(2) = bd(1)*cd(3) - cd(1)*bd(3)
	Sbcd(1) = bd(2)*cd(3) - cd(2)*bd(3)
c
	vol = ad(1)*Sbcd(1) - ad(2)*Sbcd(2) + ad(3)*Sbcd(3)
        vol = dabs(vol)/6.d0
c
	return
	end
c
c============================================================================================
c============================================================================================
c	sort4_sign.f	Version 1 9/30/2006	Patrice Koehl
c
c	This subroutine sorts the list of 4 numbers, and computes the
c	signature of the permutation
c============================================================================================
c============================================================================================
c
	subroutine sort4_sign(list,index,nswap,n)
c
	integer	i,j,k,nswap,n
	integer	list(n),index(n)
c
	save
c
	do 100 i = 1,n
		index(i) = i
100	continue
c
	nswap = 1
	do 300 i = 1,n-1
		do 200 j = i+1,n
			if(list(i).gt.list(j)) then
				k = list(i)
				list(i) = list(j)
				list(j) = k
				k = index(i)
				index(i) = index(j)
				index(j) = k
				nswap = -nswap
			endif
200		continue
300	continue
c
	return
	end
c
c============================================================================================
c============================================================================================
c	reorder_tetra.f	Version 1 9/30/2006	Patrice Koehl
c
c	This subroutine reorders the vertices of a list of tetrahedron,
c	such that now the indices are in increasing order
c
c	if iflag is set to 1, all tetrahedrons are re-ordered
c	if iflag is set to 0, only "new" tetrahedra are re-ordered: stored in list_tetra
c============================================================================================
c============================================================================================
c
	subroutine reorder_tetra(iflag,new,list_tetra)
c
c
	integer	i,j,idx,nswap,new,iflag
	integer	ntetra,ntot
	integer	index(4)
	integer	vertex(4),neighbour(4)
	integer*1 ival,nsurf(4),nindex(4)
c
	integer		tetra(4,*)
	integer		tetra_neighbour(4,*)
c
	integer*1	tetra_info(*)
	integer*1	tetra_nindex(*)
	integer		list_tetra(*)
c
        pointer         (ptr_tetra_info,tetra_info)
        pointer         (ptr_tetra_nindex,tetra_nindex)
        pointer         (ptr_tetra,tetra)
        pointer         (ptr_tetra_neighbour,tetra_neighbour)
c
	common  /tetra_zone/	ptr_tetra,ptr_tetra_neighbour,ntetra
     	common /tetra_stat/	ptr_tetra_info,ptr_tetra_nindex
c
	save
c
	if(iflag.eq.1) then
		ntot = ntetra
	else
		ntot = new
	endif
	do 400 idx = 1,ntot
c
		if(iflag.eq.1) then
			i = idx
		else
			i = list_tetra(idx)
		endif
		if(btest(tetra_info(i),1)) then
c
			do 100 j = 1,4
				vertex(j) = tetra(j,i)
100			continue
c
			call sort4_sign(vertex,index,nswap,4)
c
			do 200 j = 1,4
				neighbour(j)=tetra_neighbour(index(j),i)
				nindex(j)=ibits(tetra_nindex(i),
     1				2*(index(j)-1),2)
				nsurf(j)=ibits(tetra_info(i),
     1				2+index(j),1)
				if(neighbour(j).ne.0) then
					ival = j-1
					call mvbits(ival,0,2,
     1					tetra_nindex(neighbour(j)),
     2					2*nindex(j))
				endif
200			continue
c
			do 300 j = 1,4
				tetra(j,i)=vertex(j)
				tetra_neighbour(j,i)=neighbour(j)
				call mvbits(nindex(j),0,2,
     1				tetra_nindex(i),2*(j-1))
				call mvbits(nsurf(j),0,1,
     1				tetra_info(i),2+j)
300			continue
c
			if(nswap.eq.-1) then
				if(btest(tetra_info(i),0)) then
				   tetra_info(i)=ibclr(tetra_info(i),0)
				else
				   tetra_info(i)=ibset(tetra_info(i),0)
				endif
			endif

		endif
c
400	continue
c
	return
	end
c============================================================================================
c============================================================================================
c	Delaunay_edges.f		Version 1 6/6/2014	Patrice Koehl
c============================================================================================
c============================================================================================
c
c	This subroutine finds the list of edges in the Delaunay triangulation
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
c============================================================================================
c
	subroutine Delaunay_edges
c
c============================================================================================
c	Output:
c		- nedge:	number of edges
c		- edges:	definition of the edges
c============================================================================================
c
c============================================================================================
c	Declaration of all variables
c============================================================================================
c
	integer i,j
	integer	ntetra,nedge,nsize
	integer	npoints,nvertex
	integer idx
	integer	iedge
	integer	trig1,trig2,trig_in,trig_out,triga,trigb
	integer jtetra,ktetra,npass,ipair
c
	integer*1 ival
c
	integer ncontact(*)
	integer face_info(2,6)
	integer pair(2,6)
c
c============================================================================================
c	Information on the tetrahedra of the regular triangulation
c============================================================================================
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
	integer*1 tetra_mask(*)
c
	integer tetra(4,*)
	integer tetra_neighbour(4,*)
c
c============================================================================================
c	Information on the tetrahedra of the regular triangulation
c============================================================================================
c
	integer	edges(2,*)
c
c============================================================================================
c	Information on the vertices
c============================================================================================
c
	integer*1 vertex_info(*)
c
c============================================================================================
c	Pointers
c============================================================================================
c
	pointer	(ptr_mask,tetra_mask)
	pointer	(ptr_tetra,tetra)
	pointer (ptr_tetra_info,tetra_info)
	pointer (ptr_tetra_nindex,tetra_nindex)
	pointer (ptr_tetra_neighbour,tetra_neighbour)
	pointer (ptr_vertex_info,vertex_info)
	pointer	(ptr_edges,edges)
	pointer	(ptr_ncontact,ncontact)
c
c============================================================================================
c	Common blocks (structures)
c============================================================================================
c
	common  /tetra_zone/    ptr_tetra,ptr_tetra_neighbour,ntetra
	common  /tetra_stat/    ptr_tetra_info,ptr_tetra_nindex
	common  /vertex_zone/   ptr_vertex_info,npoints,nvertex
	common /contacts/       ptr_ncontact,ptr_edges,nsize,nedge
c
	data face_info/1,2,1,3,1,4,2,3,2,4,3,4/
	data pair/3,4,2,4,2,3,1,4,1,3,1,2/
c
	save
c
c============================================================================================
c	Find list of all edges in the alpha-complex
c============================================================================================
c
	ptr_mask = malloc(ntetra)
	do 200 i = 1,ntetra
		tetra_mask(i) = 0
200	continue
c
c	loop over all tetrahedron: if it belongs to the Delaunay triangulation,
c	check its edges; include in edge list if not seen before
c
	do 700 idx = 1,ntetra
c
	    if(.not.btest(tetra_info(idx),1)) goto 700
c
c	    Check all six edges
c					
	    do 600 iedge = 1,6
c
c		If this edge has already been considered (from another tetrahedron),
c		discard
c
		if(btest(tetra_mask(idx),iedge-1)) goto 600
c
c		iedge is the edge number in the tetrahedron idx, with:
c		iedge = 1		(c,d)
c		iedge = 2		(b,d)
c		iedge = 3		(b,c)
c		iedge = 4		(a,d)
c		iedge = 5		(a,c)
c		iedge = 6		(a,b)
c		
c		Define indices of the edge
c
		i = tetra(pair(1,iedge),idx)
		j = tetra(pair(2,iedge),idx)
c
		nedge = nedge + 1
		edges(1,nedge) = i
		edges(2,nedge) = j
c
		if(nedge.gt.0.9*nsize) then
			nsize = max(nint(1.5d0*nedge),nint(1.5d0*nsize))
			call resize_contact
		endif
c
c		trig1 and trig2 are the two faces of idx that share
c		iedge
c
		trig1 = face_info(1,iedge)
		trig2 = face_info(2,iedge)
c
c		Now we look at the star of the edge:
c
		ktetra = idx
		npass = 1
		trig_out = trig1
		jtetra = tetra_neighbour(trig_out,ktetra)
c
300		continue
c
c		Leave this side of the star if we hit the convex hull
c		in this case, the edge is not buried
c
		if(jtetra.eq.0) then
			goto 400
		endif
c
c		Leave the loop completely if we have described the
c		full cycle
c
		if(jtetra.eq.idx) goto 500
c
c		Identify the position of iedge in tetrahedron jtetra
c
		if(i.eq.tetra(1,jtetra)) then
			if(j.eq.tetra(2,jtetra)) then
				ipair = 6
			elseif(j.eq.tetra(3,jtetra)) then
				ipair = 5
			else
				ipair = 4
			endif
		elseif(i.eq.tetra(2,jtetra)) then
			if(j.eq.tetra(3,jtetra)) then
				ipair = 3
			else
				ipair = 2
			endif
		else
			ipair = 1
		endif
c
		tetra_mask(jtetra) = ibset(tetra_mask(jtetra),ipair-1)
c
c		Find out the face we "went in"
c
		ival = ibits(tetra_nindex(ktetra),2*(trig_out-1),2)
		trig_in = ival + 1
c
c		We know the two faces of jtetra that share iedge:
c
		triga = face_info(1,ipair)
		trigb = face_info(2,ipair)
c
		trig_out = triga
		if(trig_in.eq.triga) then
			trig_out = trigb
		endif
c
c
		ktetra = jtetra
		jtetra = tetra_neighbour(trig_out,ktetra)
		if(jtetra.eq.idx) goto 500
c
		goto 300
c
400		continue
c
		if(npass.eq.2) goto 500
		npass = npass + 1
		ktetra = idx
		trig_out = trig2
		jtetra = tetra_neighbour(trig_out,ktetra)
		goto 300
c
500		continue
c
600	    continue
c
700	continue
c
c============================================================================================
c	Sort list of all edges in increasing order
c============================================================================================
c
	call hpsort_two_int(edges,nedge)
c
	call free(ptr_mask)
c
	return
	end
c============================================================================================
c============================================================================================
c	Cleanup.f		Version 1 11/16/2009	Patrice Koehl
c============================================================================================
c============================================================================================
c
c	This subroutine cleans up:
c		- it frees all array space
c		- it frees gmp space
c
c============================================================================================
c============================================================================================
c
	subroutine cleanup
c
c============================================================================================
c	Declaration of all variables
c============================================================================================
c
	integer	ntetra,npoints,nvertex
c
c============================================================================================
c	Information on the tetrahedra of the regular triangulation
c============================================================================================
c
	integer*1 tetra_info(*)
	integer*1 tetra_nindex(*)
c
	integer tetra(4,*)
	integer tetra_neighbour(4,*)
c
c============================================================================================
c	Information on the vertices
c============================================================================================
c
	integer*1 vertex_info(*)
c
	real*8	radius(*),weight(*)
	real*8	coord(*)
c
c============================================================================================
c	Information on the vertices
c============================================================================================
c
	pointer	(ptr_radius,radius)
	pointer (ptr_weight,weight)
	pointer (ptr_coord,coord)
	pointer (ptr_tetra,tetra)
	pointer (ptr_tetra_nindex,tetra_nindex)
	pointer (ptr_tetra_neighbour,tetra_neighbour)
	pointer (ptr_tetra_info, tetra_info)
	pointer (ptr_vertex_info,vertex_info)
c
c============================================================================================
c	Common blocks (structures)
c============================================================================================
c
	common  /xyz_vertex/	ptr_coord,ptr_radius,ptr_weight
	common  /tetra_zone/    ptr_tetra,ptr_tetra_neighbour,ntetra
	common  /tetra_stat/    ptr_tetra_info,ptr_tetra_nindex
	common  /vertex_zone/   ptr_vertex_info,npoints,nvertex
c
	save
c
c============================================================================================
c	Free gmp space
c============================================================================================
c
	call clear_sos_gmp
c
c============================================================================================
c	Free all vertex arrays
c============================================================================================
c
	call free(ptr_vertex_info)
	call free(ptr_radius)
	call free(ptr_coord)
	call free(ptr_weight)
c
c============================================================================================
c	Free all tetra arrays
c============================================================================================
c
	call free(ptr_tetra_info)
	call free(ptr_tetra)
	call free(ptr_tetra_neighbour)
	call free(ptr_tetra_nindex)
c
c	We are done...
c
	return
	end
c
c===================================================================================
c===================================================================================
c
c	Resize_contact.f	Version 1 30/8/2010		Patrice Koehl
c
c	This subroutine resizes the array Resize_contact
c
c===================================================================================
c===================================================================================
c
	subroutine resize_contact
c
	integer	i,j
	integer	ncont,nsize
c
	integer	listcontact(2,*)
	integer	listcontact2(2,*)
	integer	ncontact(*)
c
	pointer	(ptr_listcontact,listcontact)
	pointer	(ptr_listcontact2,listcontact2)
	pointer	(ptr_ncontact,ncontact)
c
	common /contacts/ ptr_ncontact,ptr_listcontact,nsize,ncont
c
	ptr_listcontact2 = malloc(8*nsize)
c
	do 200 i = 1,ncont
		do 100 j = 1,2
			listcontact2(j,i) = listcontact(j,i)
100		continue
200	continue
c
	call free(ptr_listcontact)
c
	ptr_listcontact = ptr_listcontact2
c
	return
	end
