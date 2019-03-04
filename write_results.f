c===================================================================================
c===================================================================================
c	Write_results.f
c
c	This files contains a set of subroutines for writing:
c	- the elastic network
c	- the computed and experimental b-factors
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
c===================================================================================
c===================================================================================
c
c	Write_pml.f
c
c       This program generate a pymol file to visualize the elastic network
c
c===================================================================================
c===================================================================================
c
	subroutine write_pml(lun,natom,coordmol,corresp)
c
c===================================================================================
c===================================================================================
c
        integer i,j,k
	integer	lun,ncont,size,natom
	integer	i1,j1
	integer	corresp(*)
c
	real*8	kmin, kmax
	real*8	hmin,hmax,smin,smax,vmin,vmax
	real*8	alphah,alphas,alphav
	real*8	r,g,b,h,s,v
        real*8  a(3),c(3)
c
	integer	ncontact(*),listcontact(2,*)
c
	real*8	kval(*),radius
	real*8	coordmol(*)
c
c===================================================================================
c	Pointers and common blocks
c===================================================================================
c
	pointer	(ptr_kval,kval)
	pointer (ptr_ncontact,ncontact)
	pointer (ptr_listcontact,listcontact)
c
	common /contacts/       ptr_ncontact,ptr_listcontact,size,ncont
	common /elastic/   	ptr_kval
c
c===================================================================================
c	Formats
c===================================================================================
c
1	format('python')
2	format('python end')
3	format('from pymol.cgo import *')
4	format('from pymol import cmd')
5	format(' ')
6	format('vertex = [')
7	format('   COLOR ,',f4.1,',',f4.1,',',f4.1,',','\')
8	format('   SPHERE,',f8.3,',',f8.3,',',f8.3,',',f8.3,',')
9	format('   END','   ]')
10	format('network = [','   BEGIN, LINES, \')
11	format('   VERTEX,',f8.3,',',f8.3,',',f8.3,',','\')
12	format('cmd.load_cgo(network,''network'')')
c
	write(lun,1)
	write(lun,3)
	write(lun,4)
c
	write(lun,5) 
c
	goto 150
c
c===================================================================================
c	Write vertices of elastic network
c===================================================================================
c
	write(lun,6)
	write(lun,7) 0.0,0.0,1.0
c
	radius = 0.5d0
	do 100 k = 1,natom
		i = corresp(k)
		write(lun,8) (coordmol(3*(i-1)+j),j=1,3),radius
100	continue
	write(lun,9)
	write(lun,5)
c
150	continue
c
c===================================================================================
c	Write edges of elastic network
c===================================================================================
c
	do 200 i = 1,ncont
		if(i.eq.1) then
			kmin = kval(i)
			kmax = kval(i)
		else
			if(kval(i).lt.kmin) kmin = kval(i)
			if(kval(i).gt.kmax) kmax = kval(i)
		endif
200	continue
c
c	Build a gradient in HSV coordinates between Blue and Red
c
	r = 0.0d0
	g = 0.d0
	b = 1.d0
	call rgb2hsv(r,g,b,hmin,smin,vmin)
	r = 1.0d0
	g = 0.d0
	b = 0.d0
	call rgb2hsv(r,g,b,hmax,smax,vmax)
c
	alphah = (hmax - hmin)/(kmax-kmin)
	alphas = (smax - smin)/(kmax-kmin)
	alphav = (vmax - vmin)/(kmax-kmin)
c
	write(lun,10)
c
	do 400 i = 1,ncont
c
		i1 = corresp(listcontact(1,i))
		j1 = corresp(listcontact(2,i))
c
		do 300 j = 1,3
			a(j) = coordmol(3*(i1-1)+j)
			c(j) = coordmol(3*(j1-1)+j)
300		continue
c               
		if(kmin.eq.kmax) then
			r = 1.d0
			g = 0.d0
			b = 0.d0
		else
			h = alphah*(kval(i)-kmin) + hmin
			s = alphas*(kval(i)-kmin) + smin
			v = alphav*(kval(i)-kmin) + vmin
		endif
		call hsv2rgb(h,s,v,r,g,b)
		write(lun,7) r,g,b
		write(lun,11) (a(j),j=1,3)
		write(lun,11) (c(j),j=1,3)
		write(lun,5)
c
400     continue
c
	write(lun,9)
	write(lun,5)
c
c	write(lun,13)
	write(lun,12)
	write(lun,2)
c
	return
	end
c===================================================================================
c===================================================================================
c
c	Write_pml10.f
c
c       This program generate a pymol file to visualize the elastic network
c
c===================================================================================
c===================================================================================
c
	subroutine write_pml10(lun,natom,coordmol,corresp)
c
c===================================================================================
c===================================================================================
c
        integer i,j,npos
	integer	lun,ncont,size,natom
	integer	i1,j1
	integer	corresp(*)
c
	real*8	kmin, kmax
	real*8	hmin,hmax,smin,smax,vmin,vmax
	real*8	alphah,alphas,alphav
	real*8	r,g,b,h,s,v
        real*8  a(3),c(3)
c
	integer	idx(*)
	integer	ncontact(*),listcontact(2,*)
c
	real*8	Pcut,Kcut
	real*8	kval(*),kval2(*)
	real*8	coordmol(*)
c
c===================================================================================
c	Pointers and common blocks
c===================================================================================
c
	pointer	(ptr_kval,kval)
	pointer	(ptr_kval2,kval2)
	pointer	(ptr_idx,idx)
	pointer (ptr_ncontact,ncontact)
	pointer (ptr_listcontact,listcontact)
c
	common /contacts/       ptr_ncontact,ptr_listcontact,size,ncont
	common /elastic/   	ptr_kval
c
c===================================================================================
c	Formats
c===================================================================================
c
1	format('python')
2	format('python end')
3	format('from pymol.cgo import *')
4	format('from pymol import cmd')
5	format(' ')
7	format('   COLOR ,',f4.1,',',f4.1,',',f4.1,',','\')
9	format('   END','   ]')
10	format('network10 = [','   BEGIN, LINES, \')
11	format('   VERTEX,',f8.3,',',f8.3,',',f8.3,',','\')
12	format('cmd.load_cgo(network10,''network10'')')
14	format('cmd.set("cgo_line_width",2.0,''network10'')')
c
	write(lun,1)
	write(lun,3)
	write(lun,4)
c
	write(lun,5) 
c
c===================================================================================
c	Write edges of elastic network
c===================================================================================
c
	do 200 i = 1,ncont
		if(i.eq.1) then
			kmin = kval(i)
			kmax = kval(i)
		else
			if(kval(i).lt.kmin) kmin = kval(i)
			if(kval(i).gt.kmax) kmax = kval(i)
		endif
200	continue
c
	ptr_kval2 = malloc(8*ncont)
	ptr_idx = malloc(4*ncont)
	call dcopy(ncont,kval,1,kval2,1)
	call hpsort_key(kval2,idx,ncont)
c
	Pcut = 0.05d0
	npos = ncont - int(Pcut * ncont)
	Kcut = kval2(npos)
c
	call free(ptr_kval2)
	call free(ptr_idx)
c
c	Build a gradient in HSV coordinates between Blue and Red
c
	r = 0.0d0
	g = 0.d0
	b = 1.d0
	call rgb2hsv(r,g,b,hmin,smin,vmin)
	r = 1.0d0
	g = 0.d0
	b = 0.d0
	call rgb2hsv(r,g,b,hmax,smax,vmax)
c
	alphah = (hmax - hmin)/(kmax-kmin)
	alphas = (smax - smin)/(kmax-kmin)
	alphav = (vmax - vmin)/(kmax-kmin)
c
	write(lun,10)
c
	do 400 i = 1,ncont
c
		i1 = corresp(listcontact(1,i))
		j1 = corresp(listcontact(2,i))
c
		if(kval(i).lt.Kcut) goto 400
c
		do 300 j = 1,3
			a(j) = coordmol(3*(i1-1)+j)
			c(j) = coordmol(3*(j1-1)+j)
300		continue
c               
		if(kmin.eq.kmax) then
			r = 1.d0
			g = 0.d0
			b = 0.d0
		else
			h = alphah*(kval(i)-kmin) + hmin
			s = alphas*(kval(i)-kmin) + smin
			v = alphav*(kval(i)-kmin) + vmin
		endif
		call hsv2rgb(h,s,v,r,g,b)
		write(lun,7) r,g,b
		write(lun,11) (a(j),j=1,3)
		write(lun,11) (c(j),j=1,3)
		write(lun,5)
c
400     continue
c
	write(lun,9)
	write(lun,5)
c
c	write(lun,13)
	write(lun,12)
	write(lun,14)
	write(lun,2)
c
	return
	end
c===================================================================================
c===================================================================================
c
	subroutine write_hessian(lun, natom, sqrtmass, coordmol, bfact,
     1				fname_pdb, flag_method, Rc)
c
c===================================================================================
c===================================================================================
c
	integer	lun,ncont,size,i,j,natom
	integer	i1,j1
	integer	idx,flag_method
c
	integer	tags(*)
	integer	ncontact(*),listcontact(2,*)
c
	real*8	Rc,dx,dy,dz,distsq
	real*8	kval(*)
	real*8	sqrtmass(*),bfact(*)
	real*8	coordmol(*)
	real*8	vect(9)
c
	character	label(*)*80
	character	fname_pdb*100,name*100
	character	method(2)*10
	character	type*15
c
c===================================================================================
c	Pointers and common blocks
c===================================================================================
c
	pointer	(ptr_kval,kval)
	pointer (ptr_ncontact,ncontact)
	pointer (ptr_listcontact,listcontact)
c
	pointer	(ptr_label,label)
	pointer	(ptr_tags,tags)
c
	common /contacts/       ptr_ncontact,ptr_listcontact,size,ncont
	common /springconsts/   ptr_kval
	common /atom_tags/	ptr_label,ptr_tags
c
c===================================================================================
c	Formats
c===================================================================================
c
1	format(i6,2x,i10)
2	format(i6,5(2x,f8.3),2x,a30)
3	format(i6,2x,i6,10(2x,f14.7))
4	format('# File contains                        : ',a15)
5	format('#',/,
     1  '# Original Elastic Network for PDB file: ',a30,/,
     2  '# Computed using method                : ',a10)
6	format(
     1  '# Cutoff value Rc                      : ',f8.3)
7	format(
     1  '#',/,
     2  '# line 1: number N of (active) atoms, and number E of edges',/,
     3  '# line 2 to N+1: informations about atoms:',/,
     4  '#               - index (from the input PDB)',/,
     5  '#               - mass (should not be 0!)',/,
     6  '#               - X, Y, Z coordinates,'/,
     7  '#               - Bfactor (set to 0 if not available)',/,
     8  '#               - All 30 first characters from line in PDB',/,
     9  '# line N+2 to N+E: information on edges of Elastic Network:',/,
     &  '#               - indices of the two atoms (starts at one!)',/,
     &  '#               - kval: strength (if traditional EN)',/,
     &  '#               - Upper triangular part of cov. matrix',/,
     &  '#                 (for Hessian)',/,
     &  '#')
c
c===================================================================================
c	Write header: information on elastic network
c===================================================================================
c
	idx=index(fname_pdb,'.') -1
	name=fname_pdb(1:idx)//'                              '
	method(1)='Cutoff    '
	method(2)='Delaunay  '
	type = 'Hessian        '
	write(lun,4) type
	write(lun,5) name,method(flag_method)
	if(flag_method.eq.1) then
		write(lun,6) Rc
	endif
	write(lun,7)
	write(lun,1) natom, ncont
c
c===================================================================================
c	Write atoms: current active atoms with all necessary information
c===================================================================================
c
	do 100 i = 1,natom
		write(lun,2) i,sqrtmass(i)*sqrtmass(i),
     1			(coordmol(3*(i-1)+j),j=1,3),bfact(i),
     2			label(i)(1:30)
c     1			(coordmol(3*(i-1)+j),j=1,3)
100	continue
c
c
c===================================================================================
c	Write edges: atom indices, kval (strength), and covar for Hessian
c===================================================================================
c
	do 200 i = 1,ncont
		i1 = listcontact(1,i)
		j1 = listcontact(2,i)
		dx  = coordmol(3*(i1-1)+1) - coordmol(3*(j1-1)+1)
		dy  = coordmol(3*(i1-1)+2) - coordmol(3*(j1-1)+2)
		dz  = coordmol(3*(i1-1)+3) - coordmol(3*(j1-1)+3)
		distsq = dx*dx+dy*dy+dz*dz
		vect(1) = -kval(i)*dx*dx/distsq
		vect(2) = -kval(i)*dx*dy/distsq
		vect(3) = -kval(i)*dx*dz/distsq
		vect(4) = vect(2)
		vect(5) = -kval(i)*dy*dy/distsq
		vect(6) = -kval(i)*dy*dz/distsq
		vect(7) = vect(3)
		vect(8) = vect(6)
		vect(9) = -kval(i)*dz*dz/distsq
		write(lun,3) i1,j1,kval(i),(vect(j),j=1,9)
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c
	subroutine write_eln(lun, natom, corresp, sqrtmass, coordmol, 
     1		bfact,	fname_pdb, flag_method, Rc)
c
c===================================================================================
c===================================================================================
c
	integer	lun,ncont,size,i,j,natom
	integer	i1,j1,iat
	integer	idx,flag_method
c
	integer	corresp(*)
	integer	ncontact(*),listcontact(2,*)
c
	real*8	Rc
	real*8	kval(*)
	real*8	sqrtmass(*),bfact(*)
	real*8	coordmol(*)
c
	character	label(*)*30
	character	fname_pdb*100,name*100
	character	method(2)*10
	character	type*15
c
c===================================================================================
c	Pointers and common blocks
c===================================================================================
c
	pointer	(ptr_kval,kval)
	pointer (ptr_ncontact,ncontact)
	pointer (ptr_listcontact,listcontact)
c
	pointer	(ptr_label,label)
c
	common /contacts/       ptr_ncontact,ptr_listcontact,size,ncont
	common /elastic/        ptr_kval
	common /tags/		ptr_label
c
c===================================================================================
c	Formats
c===================================================================================
c
1	format(i6,2x,i10)
2	format(i6,5(2x,f8.3),2x,a30)
3	format(i6,2x,i6,10(2x,f10.6))
4	format('# File contains                        : ',a15)
5	format('#',/,
     1  '# Original Elastic Network for PDB file: ',a30,/,
     2  '# Computed using method                : ',a10)
6	format(
     1  '# Cutoff value Rc                      : ',f8.3)
7	format(
     1  '#',/,
     2  '# line 1: number N of (active) atoms, and number E of edges',/,
     3  '# line 2 to N+1: informations about atoms:',/,
     4  '#               - index (from the input PDB)',/,
     5  '#               - mass (should not be 0!)',/,
     6  '#               - X, Y, Z coordinates,'/,
     7  '#               - Bfactor (set to 0 if not available)',/,
     8  '#               - All 30 first characters from line in PDB',/,
     9  '# line N+2 to N+E: information on edges of Elastic Network:',/,
     &  '#               - indices of the two atoms (starts at one!)',/,
     &  '#               - kval: strength (if traditional EN)',/,
     &  '#')
c
c===================================================================================
c	Write header: information on elastic network
c===================================================================================
c
	idx=index(fname_pdb,'.') -1
	name=fname_pdb(1:idx)//'                              '
	method(1)='Cutoff    '
	method(2)='Delaunay  '
	type = 'Elastic Network'
	write(lun,4) type
	write(lun,5) name,method(flag_method)
	if(flag_method.eq.1) then
		write(lun,6) Rc
	endif
	write(lun,7)
	write(lun,1) natom, ncont
c
c===================================================================================
c	Write atoms: current active atoms with all necessary information
c===================================================================================
c
	do 100 i = 1,natom
		iat = corresp(i)
		write(lun,2) i,sqrtmass(iat)*sqrtmass(iat),
     1			(coordmol(3*(iat-1)+j),j=1,3),bfact(iat),
     2			label(iat)(1:30)
100	continue
c
c
c===================================================================================
c	Write edges: atom indices, kval (strength), and covar for Hessian
c===================================================================================
c
	do 200 i = 1,ncont
		i1 = listcontact(1,i)
		j1 = listcontact(2,i)
		write(lun,3) i1,j1,kval(i)
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c
c	rgb2hsv.f
c
c	Convert a color in RGB (with each color in [0,1]) to its corresponding HSV
c	code
c
c===================================================================================
c===================================================================================
c
	subroutine rgb2hsv(r,g,b,h,s,v)
c
	real*8	r,g,b,h,s,v
	real*8	a,c
	real*8	minRGB, maxRGB
c
	minRGB = min(r,min(g,b))
	maxRGB = max(r,max(g,b))
c
	if(minRGB.eq.maxRGB) then
		h = 0.d0
		s = 0.d0
		v = minRGB
		return
	endif
c
	if(minRGB.eq.r) then
		a = g - b
		c = 3
	elseif(minRGB.eq.b) then
		a = r-g
		c = 1
	else
		a = b - r
		c = 5
	endif
c
	h = 60.d0*(c - a/(maxRGB-minRGB))
	h = h/360.d0
	s = (maxRGB - minRGB) / maxRGB
	v = maxRGB
c
	return
	end
c===================================================================================
c===================================================================================
c
c	hsv2rgb.f
c
c	Convert a color in HSV to RGB (with each color in [0,1])
c
c===================================================================================
c===================================================================================
c
	subroutine hsv2rgb(h,s,v,r,g,b)
c
	integer	var_i
c
	real*8	r,g,b,h,s,v
	real*8	var_h,var_1,var_2,var_3
c
	if(s.eq.0.d0) then
		r = v
		g = v
		b = v
		return
	endif
c
	var_h = h*6
	if(var_h.eq.6.d0) var_h = 0.d0
c
	var_i = floor(var_h)
	var_1 = v*(1-s)
	var_2 = v*(1-s*(var_h-var_i))
	var_3 = v*(1-s*(1-(var_h-var_i)))
c
	if(var_i.eq.0) then
		r = v
		g = var_3
		b = var_1
	elseif(var_i.eq.1) then
		r = var_2
		g = v
		b = var_1
	elseif(var_i.eq.2) then
		r = var_1
		g = v
		b = var_3
	elseif(var_i.eq.3) then
		r = var_1
		g = var_2
		b = v
	elseif(var_i.eq.4) then
		r = var_3
		g = var_1
		b = v
	else
		r = v
		g = var_1
		b = var_2
	endif
c
	return
	end
