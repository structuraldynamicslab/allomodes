c===================================================================================
c===================================================================================
c	Hessian_Go.f
c
c       This program computes the Hessian of the Go elastic energy for a given 
c	conformation of the protein, and stores it as a sparse matrix, in compact
c       row format
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
	subroutine prepare_Hessian(nat,spr_jdx,spr_row,spr_hessian)
c
c===================================================================================
c===================================================================================
c
c	This subroutine prepares the storage for the Hessian for the Go elastic 
c	model of the molecule,	storing in to sparse row format.
c
c	Input:
c		nat :	Number of atoms considered
c	Input (from common blocks):
c		ncontact: # of contacts for each atom
c		listcontact: list of edges in elastic network
c
c	Output:
c		spr_jdx:   Index array for bookkeeping on sparse array
c		spr_row:   second index array for bookkeeping on sparse array
c
c
c	The Hessian H is the matrix of second order derivatives of E with
c	respect to the 3*Nat coordinates of the Nat atoms.
c	H is theoretically of size (3*Nat) * (3*Nat), but because
c	we use a cutoff, H is sparse, and it is stored as such
c
c	Let N = 3*Nat
c	Storing scheme:
c	We use three arrays:
c		- spr_hessian:	contains the non zero values of the Hessian
c				matrix (only the upper part, as the
c				Hessian is symmetric)
c		- spr_row:	pointer of the position of the beginning of each row
c				in spr_hessian
c		- spr_jdx:	integer array of pointers that allows
c				fast retrieval of the values along each row (compressed)
c
c	The storage rules are:
c
c	- value j of spr_row, with j = 1 to N, is the index of the position
c	in spr_hessian of the first non zero element of row j of H
c	Note that 
c			spr_row(1) = 1
c			spr_row(N+1) gives the size of spr_hessian and
c					spr_jdx
c	- spr_hessian(k) are the non-zero elements of H, where H is scanned row by row
c	- spr_jdx(k) is the column number of the corresponding element in spr_hessian
c
c	Example : if we consider the matrix:
c
c		2.	3.	0.	0.
c		1.	3.	0.	4.
c		0.	2.	0.	2.
c		1.	0.	0.	4.
c
c	Then:
c	spr_row    =[1  3  6  8]
c	spr_hessian=[2  3  1  3  4  2  2  1  4 ]
c	spr_jdx =   [1  2  1  2  3  2  4  1  4 ]
c
c	where x can be set to any value (not used) 
c
c	For our application, since H is symmetric, we store it as 
c	if it was triangular
c
c	Note that an atom i corresponds to 3 rows in the Hessian (for derivatives
c	with respect to its 3 coordinates.
c	Each contact i,j generate three columns in the sparse representation
c	of H.
c
c===================================================================================
c===================================================================================
c
	integer	nat,ipos,nsize,ncont,npos
	integer	iat,jat,ja
	integer spr_jdx(*),spr_row(*)
	integer ncontact(*),listcontact(2,*)	
c
	real*8	spr_hessian(*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
c
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
c
	ipos = 0
	npos = 0
c
c===================================================================================
c	Build adjacency information
c===================================================================================
c
	do 400 iat = 1,nat
c
c		Initialise index for the three rows 3*iat-2,3*iat-1 and
c		3*iat
c
		spr_row(3*iat-2) = ipos + 1
c
		spr_jdx(ipos + 1) = 3*iat-2
		spr_jdx(ipos + 2) = 3*iat-1
		spr_jdx(ipos + 3) = 3*iat
c
		spr_hessian(ipos+1) = 0
		spr_hessian(ipos+2) = 0
		spr_hessian(ipos+3) = 0
c
		ipos = ipos + 3
c
		do 100 ja = 1,ncontact(iat)
c
			jat = listcontact(2,ja+npos)
c
			spr_jdx(ipos + 1) = 3*jat-2
			spr_jdx(ipos + 2) = 3*jat-1
			spr_jdx(ipos + 3) = 3*jat
c
			ipos = ipos + 3
c
100		continue
c
		spr_row(3*iat-1) = ipos + 1
c
		spr_jdx(ipos+1) = 3*iat-1
		spr_jdx(ipos+2) = 3*iat
c
		spr_hessian(ipos+1) = 0
		spr_hessian(ipos+2) = 0
c
		ipos = ipos + 2
c
		do 200 ja = 1,ncontact(iat)
c
			jat = listcontact(2,npos+ja)
c
			spr_jdx(ipos + 1) = 3*jat-2
			spr_jdx(ipos + 2) = 3*jat-1
			spr_jdx(ipos + 3) = 3*jat
c
			ipos = ipos + 3
c
200		continue
c
		spr_row(3*iat) = ipos + 1
c
		spr_jdx(ipos+1) = 3*iat
		spr_hessian(ipos+1) = 0
c
		ipos = ipos + 1
c
		do 300 ja = 1,ncontact(iat)
c
			jat = listcontact(2,npos+ja)
c
			spr_jdx(ipos + 1) = 3*jat-2
			spr_jdx(ipos + 2) = 3*jat-1
			spr_jdx(ipos + 3) = 3*jat
c
			ipos = ipos + 3
c
300		continue
c
		npos = npos + ncontact(iat)
c
400	continue
c
	spr_row(3*nat+1) = ipos + 1
c
	return
	end
c===================================================================================
c===================================================================================
c	fill_hessian_go.f
c
c	This subroutine computes the Hessian for a Go model
c
c	Copyright (C) 2016 Patrice Koehl
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
	subroutine Fill_hessian_GO(natom,corresp,sqrtmass,coord,
     1		spr_jdx,spr_row,spr_hessian)
c
c===================================================================================
c===================================================================================
c
	integer	i, j, k, l
	integer	i_k, j_k, ip1_k
	integer	jat, irow, add1, add2
	integer	natom,i_start,iat,ipos
	integer	nsize,ncont,npos
        integer spr_row(*),spr_jdx(*)
	integer ncontact(*),listcontact(2,*)	
	integer	corresp(*)
c
	real*8	epsilon, K_bond, K_angle, K_dih1, K_dih3, K_dih
	real*8	dist, sum, dist1, dist2, dist1sq, dist2sq
	real*8	num, den, angle, fact1, fact
	real*8	f1, f2, f3
	real*8	facti,facti2,factj
	real*8	diff(3),diff1(3),diff2(3),diff3(3)
	real*8	vect1(3), vect2(3)
	real*8	v1a(3), v1b(3), v1c(3), v2a(3), v2b(3), v2c(3)
	real*8	d_num(12), d_dist1(12), d_dist2(12)
	real*8	d_angle(12)
	real*8	coord(*), spr_hessian(*),sqrtmass(*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
c
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
c
c===================================================================================
c	Set all constants
c===================================================================================
c
	epsilon = 0.36d0
	K_bond = 100.d0*epsilon
	K_angle = 20.d0*epsilon
	K_dih1  = epsilon
	K_dih3  = 0.5d0*epsilon
c
c===================================================================================
c	First term: bonds (between residues i and i + 1)
c===================================================================================
c
	fact1 = -2.d0*K_bond
	do 400 i = 1,natom-1
c
		do 100 j = 1,3
			diff(j) = coord(3*(i-1)+j)-coord(3*i+j)
100		continue
		call dotvect(diff,diff,dist)
		fact = fact1/dist
c
		do 300 k = 1,3
			i_k = spr_row(3*(i-1)+k) + 3-k
			f1 = fact*diff(k)
			do 200 j = 1,3
				spr_hessian(i_k+j) = f1*diff(j)
200			continue
300		continue
c
400	continue
c
c===================================================================================
c	Second term: angles (between residues i, i + 1, and i+2)
c===================================================================================
c
	fact1 = 2.d0*K_angle
c
	do 1300 i = 1, natom-2
c
		do 500 j = 1,3
			diff1(j) = coord(3*i+j)-coord(3*(i-1)+j)
			diff2(j) = coord(3*i+j)-coord(3*(i+1)+j)
500		continue
		call dotvect(diff1,diff2,num)
		call dotvect(diff1,diff1,dist1sq)
		call dotvect(diff2,diff2,dist2sq)
c
		dist1 = sqrt(dist1sq)
		dist2 = sqrt(dist2sq)
		den = dist1*dist2
		angle = num / den
c
		fact = fact1/(1.d0-angle*angle)
c
		do 600 j = 1,3
			d_num(j) = -diff2(j)
			d_num(j+3) = diff2(j) + diff1(j)
			d_num(j+6) = -diff1(j)
			d_dist1(j) = -diff1(j)-diff1(j)
			d_dist1(j+3) = diff1(j)+diff1(j)
			d_dist1(j+6) = 0.d0
			d_dist2(j) = 0.d0
			d_dist2(j+3) = diff2(j)+diff2(j)
			d_dist2(j+6) = -diff2(j)-diff2(j)
600		continue
c
		f1 = 1.d0/den
		f2 = -0.5d0*angle/dist1sq
		f3 = -0.5d0*angle/dist2sq
		do 700 j = 1,9
			d_angle(j) = f1*d_num(j)+f2*d_dist1(j)+
     1					f3*d_dist2(j)
700		continue
c
		do 1000 k = 1,3
			i_k = spr_row(3*(i-1)+k) + 3-k
			f1 = fact*d_angle(k)
			do 800 j = 1,3
				spr_hessian(i_k+j)= spr_hessian(i_k+j)
     1					+ f1*d_angle(j+3)
800			continue
			do 900 j = 4,6
				spr_hessian(i_k+j)= spr_hessian(i_k+j)
     1					+ f1*d_angle(j+3)
900			continue
1000		continue
c
		do 1200 k = 1,3
			ip1_k = spr_row(3*i+k) + 3-k
			f1 = fact*d_angle(k+3)
			do 1100 j = 1,3
				j_k = ip1_k + j
				spr_hessian(j_k)=spr_hessian(j_k)
     1					+ f1*d_angle(j+6)
1100			continue
1200		continue
c
1300	continue
c
c===================================================================================
c	Third term: dihedral (between residues i, i + 1, i+2, and i+4)
c===================================================================================
c
	K_dih = K_dih1 + 9 * K_dih3
c
	do 2600 i = 1, natom-3
c
		do 1400 j = 1,3
			diff1(j) = coord(3*i+j)-coord(3*(i-1)+j)
			diff2(j) = coord(3*(i+1)+j)-coord(3*i+j)
			diff3(j) = coord(3*(i+2)+j)-coord(3*(i+1)+j)
1400		continue
c
		call crossvect(diff1,diff2,vect1)
		call crossvect(diff2,diff3,vect2)
c
		call dotvect(vect1,vect2,num)
		call dotvect(vect1,vect1,dist1sq)
		call dotvect(vect2,vect2,dist2sq)
c
		dist1 = sqrt(dist1sq)
		dist2 = sqrt(dist2sq)
		den = dist1*dist2
		angle = num / den
c
		call crossvect(vect1,diff1,v1a)
		call crossvect(vect1,diff2,v1b)
		call crossvect(vect1,diff3,v1c)
		call crossvect(vect2,diff1,v2a)
		call crossvect(vect2,diff2,v2b)
		call crossvect(vect2,diff3,v2c)
c
		do 1500 j = 1,3
			d_num(j) = v2b(j)
			d_num(j+3) = -v2b(j) + v1c(j) -v2a(j)
			d_num(j+6) = v2a(j) - v1b(j) - v1c(j)
			d_num(j+9) = v1b(j)
			d_dist1(j) = v1b(j) + v1b(j)
			d_dist1(j+3) = -v1b(j) - v1b(j) -v1a(j) -v1a(j)
			d_dist1(j+6) = v1a(j) + v1a(j)
			d_dist1(j+9) = 0.d0
			d_dist2(j) = 0.d0
			d_dist2(j+3) = v2c(j) + v2c(j)
			d_dist2(j+6) = -v2c(j) - v2c(j) -v2b(j) -v2b(j)
			d_dist2(j+9) = v2b(j) + v2b(j)
1500		continue
c
		f1 = 1.d0/den
		f2 = -0.5d0*angle/dist1sq
		f3 = -0.5d0*angle/dist2sq
		do 1600 j = 1,12
			d_angle(j) = f1*d_num(j)+f2*d_dist1(j)+
     1					f3*d_dist2(j)
1600		continue
c
		fact = K_dih /(1.d0-angle*angle)
c
		do 2000 k = 1,3
			i_k = spr_row(3*(i-1)+k) + 3-k
			f1 = fact*d_angle(k)
			do 1700 j = 1,3
				spr_hessian(i_k+j)=spr_hessian(i_k+j)
     1					+f1*d_angle(j+3)
1700			continue
			do 1800 j = 4,6
				spr_hessian(i_k+j)=spr_hessian(i_k+j)
     1					+f1*d_angle(j+3)
1800			continue
			do 1900 j = 7,9
				spr_hessian(i_k+j)=spr_hessian(i_k+j)
     1					+f1*d_angle(j+3)
1900			continue
2000		continue
c
		do 2300 k = 1,3
			i_k = spr_row(3*i+k) + 3-k
			f1 = fact*d_angle(k+3)
			do 2100 j = 1,3
				spr_hessian(i_k+j)=spr_hessian(i_k+j)
     1					+f1*d_angle(j+6)
2100			continue
			do 2200 j = 4,6
				spr_hessian(i_k+j)=spr_hessian(i_k+j)
     1					+f1*d_angle(j+6)
2200			continue
2300		continue
c
		do 2500 k = 1,3
			i_k = spr_row(3*(i+1)+k) + 3-k
			f1 = fact*d_angle(k+6)
			do 2400 j = 1,3
				spr_hessian(i_k+j)=spr_hessian(i_k+j)
     1					+f1*d_angle(j+9)
2400			continue
2500		continue
c
2600	continue
c
c===================================================================================
c	Fourth term: non-local interactions (vdw + electrostatics)
c===================================================================================
c
	fact1 = -120.d0*epsilon
	npos = 0
	do 3100 i = 1, natom-4 
c
	   do 3000 jat = 1,ncontact(i)
c
		j = listcontact(2,npos+jat)
		if(j-i.le.3) goto 3000
c
		do 2700 k = 1,3
			diff(k) = coord(3*(i-1)+k)-coord(3*(j-1)+k)
2700		continue
c
		call dotvect(diff,diff,dist)
		fact = fact1/(dist*dist)
c
		do 2900 k = 1,3
			i_k = spr_row(3*(i-1)+k) + 3-k +3*(jat-1)
			f1 = fact*diff(k)
			do 2800 l = 1,3
				spr_hessian(i_k+l) = 
     1				spr_hessian(i_k+l)+f1*diff(l)
2800			continue
2900		continue
c
3000	    continue
c
	    npos = npos + ncontact(i)
c
3100	continue
c
c===================================================================================
c	Fill in diagonal blocks
c===================================================================================
c
	do 3500 i = 1,natom
	    do 3400 j = 1,3
c
		irow = 3*(i-1)+j
		i_k = spr_row(irow)
		i_start = i_k + 3-j 
c
		do 3300 k = 1, 3
c
			sum = 0.d0
			do 3200 l = i_start+k,spr_row(irow+1)-1,3
				f1 = spr_hessian(l)
				sum = sum + f1
				j_k = spr_jdx(l)
				if(j.ge.k) then
					j_k = spr_row(j_k) + j - k
					spr_hessian(j_k) = 
     1					spr_hessian(j_k)-f1
				endif
3200			continue
			if(k.ge.j) spr_hessian(i_k+k-j) = 
     1				spr_hessian(i_k+k-j)-sum
c
3300		continue
c
3400	    continue
3500	continue
c
c===================================================================================
c	Add mass information
c===================================================================================
c
	do 3800 iat = 1,natom
c
		facti = 1/sqrtmass(corresp(iat))
		facti2 = facti*facti
c
		ipos = spr_row(3*iat-2)
		add1 = spr_row(3*iat-1) - spr_row(3*iat-2) - 1
		add2 = spr_row(3*iat) - spr_row(3*iat-1) - 1
c
		spr_hessian(ipos) = spr_hessian(ipos)*facti2
		spr_hessian(ipos+1) = spr_hessian(ipos+1)*facti2
		spr_hessian(ipos+2) = spr_hessian(ipos+2)*facti2
		spr_hessian(ipos+add1+1) = 
     1			spr_hessian(ipos+add1+1)*facti2
		spr_hessian(ipos+add1+2) = 
     1			spr_hessian(ipos+add1+2)*facti2
		spr_hessian(ipos+add1+add2+2)=
     1			spr_hessian(ipos+add1+add2+2)*facti2
c
		do 3700 j = spr_row(3*iat-2) + 3,spr_row(3*iat-1)-1,3
c
			jat = (spr_jdx(j)+2)/3
	
			factj = facti/sqrtmass(corresp(jat))
c
			do 3600 l = 1,3
				spr_hessian(j+l-1) = spr_hessian(j+l-1)
     1					*factj
				spr_hessian(j+add1+l-1)=
     1				    spr_hessian(j+add1+l-1)*factj
				spr_hessian(j+add1+add2+l-1)=
     1				    spr_hessian(j+add1+add2+l-1)*factj
3600			continue
c
3700		continue
c
3800	continue
c
	return
	end
c
c===================================================================================
c===================================================================================
c	hessian_go_full.f
c
c	This subroutine computes the Hessian for a Go model
c
c	Copyright (C) 2016 Patrice Koehl
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
	subroutine hessian_go_full(natom,corresp,sqrtmass,coord,hessian)
c
c===================================================================================
c===================================================================================
c
	integer	i, j, k, l, jat
	integer	natom
	integer	nsize,ncont,npos
	integer ncontact(*),listcontact(2,*)	
	integer	corresp(*)
c
	real*8	epsilon, K_bond, K_angle, K_dih1, K_dih3, K_dih
	real*8	dist, sum
	real*8	dist1sq, dist2sq, dist1, dist2
	real*8	num, den, angle, fact1, fact
	real*8	f1, f2, f3
	real*8	diff(3), diff1(3), diff2(3), diff3(3)
	real*8	vect1(3), vect2(3)
	real*8	v1a(3), v1b(3), v1c(3), v2a(3), v2b(3), v2c(3)
	real*8	d_num(12), d_dist1(12), d_dist2(12)
	real*8	d_angle(12)
	real*8	coord(*), hessian(3*natom,3*natom), sqrtmass(*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
c
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
c
c===================================================================================
c	Set all constants
c===================================================================================
c
	epsilon = 0.36d0
	K_bond = 100.d0*epsilon
	K_angle = 20.d0*epsilon
	K_dih1  = epsilon
	K_dih3  = 0.5d0*epsilon
c
c===================================================================================
c	First term: bonds (between residues i and i + 1)
c===================================================================================
c
	fact1 = -2.d0*K_bond
	do 400 i = 1,natom-1
c
		do 100 j = 1,3
			diff(j) = coord(3*(i-1)+j)-coord(3*i+j)
100		continue
		call dotvect(diff,diff,dist)
		fact = fact1/dist
c
		do 300 k = 1,3
			f1 = fact*diff(k)
			do 200 j = 1,3
				hessian(3*(i-1)+k,3*i+j) = f1*diff(j)
200			continue
300		continue
c
400	continue
c
c===================================================================================
c	Second term: angles (between residues i, i + 1, and i+2)
c===================================================================================
c
	fact1 = 2.d0*K_angle
c
	do 1300 i = 1, natom-2
c
		do 500 j = 1,3
			diff1(j) = coord(3*i+j)-coord(3*(i-1)+j)
			diff2(j) = coord(3*i+j)-coord(3*(i+1)+j)
500		continue
		call dotvect(diff1,diff2,num)
		call dotvect(diff1,diff1,dist1sq)
		call dotvect(diff2,diff2,dist2sq)
c
		dist1 = sqrt(dist1sq)
		dist2 = sqrt(dist2sq)
		den = dist1*dist2
		angle = num / den
c
		fact = fact1/(1.d0-angle*angle)
c
		do 600 j = 1,3
			d_num(j) = -diff2(j)
			d_num(j+3) = diff2(j) + diff1(j)
			d_num(j+6) = -diff1(j)
			d_dist1(j) = -diff1(j)-diff1(j)
			d_dist1(j+3) = diff1(j)+diff1(j)
			d_dist1(j+6) = 0.d0
			d_dist2(j) = 0.d0
			d_dist2(j+3) = diff2(j)+diff2(j)
			d_dist2(j+6) = -diff2(j)-diff2(j)
600		continue
c
		f1 = 1.d0/den
		f2 = -0.5d0*angle/dist1sq
		f3 = -0.5d0*angle/dist2sq
		do 700 j = 1,9
			d_angle(j) = f1*d_num(j)+f2*d_dist1(j)+
     1					f3*d_dist2(j)
700		continue
c
		do 1000 k = 1,3
			f1 = fact*d_angle(k)
			do 800 j = 1,6
				hessian(3*(i-1)+k,3*i+j)=
     1				hessian(3*(i-1)+k,3*i+j)
     2					+ f1*d_angle(j+3)
800			continue
1000		continue
c
		do 1200 k = 1,3
			f1 = fact*d_angle(k+3)
			do 1100 j = 1,3
				hessian(3*i+k,3*i+3+j)=
     1				hessian(3*i+k,3*i+3+j)
     2				+ f1*d_angle(j+6)
1100			continue
1200		continue
c
1300	continue
c
c===================================================================================
c	Third term: dihedral (between residues i, i + 1, i+2, and i+4)
c===================================================================================
c
	K_dih = K_dih1 + 9 * K_dih3
c
	do 2600 i = 1, natom-3
c
		do 1400 j = 1,3
			diff1(j) = coord(3*i+j)-coord(3*(i-1)+j)
			diff2(j) = coord(3*(i+1)+j)-coord(3*i+j)
			diff3(j) = coord(3*(i+2)+j)-coord(3*(i+1)+j)
1400		continue
c
		call crossvect(diff1,diff2,vect1)
		call crossvect(diff2,diff3,vect2)
c
		call dotvect(vect1,vect2,num)
		call dotvect(vect1,vect1,dist1sq)
		call dotvect(vect2,vect2,dist2sq)
c
		dist1 = sqrt(dist1sq)
		dist2 = sqrt(dist2sq)
		den = dist1*dist2
		angle = num / den
c
		call crossvect(vect1,diff1,v1a)
		call crossvect(vect1,diff2,v1b)
		call crossvect(vect1,diff3,v1c)
		call crossvect(vect2,diff1,v2a)
		call crossvect(vect2,diff2,v2b)
		call crossvect(vect2,diff3,v2c)
c
		do 1500 j = 1,3
			d_num(j) = v2b(j)
			d_num(j+3) = -v2b(j) + v1c(j) -v2a(j)
			d_num(j+6) = v2a(j) - v1b(j) - v1c(j)
			d_num(j+9) = v1b(j)
			d_dist1(j) = v1b(j) + v1b(j)
			d_dist1(j+3) = -v1b(j) - v1b(j) -v1a(j) -v1a(j)
			d_dist1(j+6) = v1a(j) + v1a(j)
			d_dist1(j+9) = 0.d0
			d_dist2(j) = 0.d0
			d_dist2(j+3) = v2c(j) + v2c(j)
			d_dist2(j+6) = -v2c(j) - v2c(j) -v2b(j) -v2b(j)
			d_dist2(j+9) = v2b(j) + v2b(j)
1500		continue
c
		f1 = 1.d0/den
		f2 = -0.5d0*angle/dist1sq
		f3 = -0.5d0*angle/dist2sq
		do 1600 j = 1,12
			d_angle(j) = f1*d_num(j)+f2*d_dist1(j)+
     1					f3*d_dist2(j)
1600		continue
c
		fact = K_dih /(1.d0-angle*angle)
c
		do 2000 k = 1,3
			f1 = fact*d_angle(k)
			do 1700 j = 1,9
				hessian(3*(i-1)+k,3*i+j)=
     1				hessian(3*(i-1)+k,3*i+j)
     2					+ f1*d_angle(j+3)
1700			continue
2000		continue
c
		do 2300 k = 1,3
			f1 = fact*d_angle(k+3)
			do 2100 j = 1,6
				hessian(3*i+k,3*i+3+j)=
     1				hessian(3*i+k,3*i+3+j)
     2					+ f1*d_angle(j+6)
2100			continue
2300		continue
c
		do 2500 k = 1,3
			f1 = fact*d_angle(k+6)
			do 2400 j = 1,3
				hessian(3*i+3+k,3*i+6+j)=
     1				hessian(3*i+3+k,3*i+6+j)
     2					+ f1*d_angle(j+9)
2400			continue
2500		continue
c
2600	continue
c
c===================================================================================
c	Fourth term: non-local interactions (vdw + electrostatics)
c===================================================================================
c
	fact1 = -120.d0*epsilon
	npos = 0
	do 3100 i = 1, natom-4 
c
	   do 3000 jat = 1,ncontact(i)
c
		j = listcontact(2,npos+jat)
		if(j-i.le.3) goto 3000
c
		do 2700 k = 1,3
			diff(k) = coord(3*(i-1)+k)-coord(3*(j-1)+k)
2700		continue
c
		call dotvect(diff,diff,dist)
		fact = fact1/(dist*dist)
c
		do 2900 k = 1,3
			f1 = fact*diff(k)
			do 2800 l = 1,3
				hessian(3*(i-1)+k,3*(j-1)+l)=
     1				hessian(3*(i-1)+k,3*(j-1)+l)
     2					+ f1*diff(l)
2800			continue
2900		continue
c
3000	    continue
c
	    npos = npos + ncontact(i)
c
3100	continue
c
c===================================================================================
c	Make symmetric, and fill in diagonal
c===================================================================================
c
	do 3300 i = 2,3*natom
		do 3200 j = 1,i-1
			hessian(i,j) = hessian(j,i)
3200		continue
3300	continue
c
	do 3600 i = 1,3*natom
c
		l = i/3 + 1
		if(mod(i,3).eq.0) l = l-1
		do 3500 j = 1,3
			sum = 0.d0
			do 3400 k = j, 3*natom, 3
				sum = sum + hessian(i,k)
3400			continue
			hessian(i,3*(l-1)+j) = -sum
3500		continue
3600	continue
c
c===================================================================================
c===================================================================================
c
	do 3800 i = 1, 3*natom
		k = i/3 + 1
		if(mod(i,3).eq.0) k = k-1
		do 3700 j = 1, 3*natom
			l = j/3 + 1
			if(mod(j,3).eq.0) l = l-1
			fact = 1.d0/(sqrtmass(corresp(k))*
     1				sqrtmass(corresp(l)))
			hessian(i,j) = hessian(i,j)*fact
3700		continue
3800	continue
c
	return
	end
c===================================================================================
c===================================================================================
c
	subroutine contact_GO(natom)
c
c===================================================================================
c===================================================================================
c
	integer	i,j,j1
	integer	natom,ncont,nsize,ncont2,nsize2,nc,npos,ncount
	integer ncontact(*),listcontact(2,*)
	integer	listcontact2(2,*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
	pointer	(ptr_listcontact2,listcontact2)
c
c===================================================================================
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
c===================================================================================
c
	nsize2 = nsize + 4*natom
	ptr_listcontact2 = malloc(4*2*nsize2)
	nsize = nsize2
c
	ncont2 = 0
	nc = 0
	do 300 i = 1,natom-1
c
		do 100 j = 1,min(3,natom-i)
			ncont2 = ncont2 + 1
			listcontact2(1,ncont2) = i
			listcontact2(2,ncont2) = i + j
100		continue
c
		do 200 j = 1,ncontact(i)
			j1 = listcontact(2,nc+j)
			if(j1-i.le.3) goto 200
			ncont2 = ncont2 + 1
			listcontact2(1,ncont2) = i
			listcontact2(2,ncont2) = j1
200		continue
c
		nc = nc + ncontact(i)
c
300	continue
c
	ncont = ncont2
	call free(ptr_listcontact)
	ptr_listcontact = ptr_listcontact2
c
	do 400 i = 1,natom
		ncontact(i) = 0
400	continue
c
	npos = 0
	do 500 j = 1,ncont
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
500	continue
	ncontact(npos) = ncount
c
	return
	end
c===================================================================================
c===================================================================================
c	hessian_elastic_full.f
c
c	This subroutine computes the Hessian for an elastic model
c
c	Copyright (C) 2016 Patrice Koehl
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
	subroutine hessian_elastic_full(natom,corresp,sqrtmass,
     1			coord,hessian)
c
c===================================================================================
c===================================================================================
c
	integer	i, j, k, l
	integer ic
	integer	natom
	integer	nsize,ncont
	integer ncontact(*),listcontact(2,*)	
	integer	corresp(*)
c
	real*8	dist, sum, Kconst
	real*8	fact
	real*8	f1
	real*8	diff(3)
	real*8	coord(*), kval(*), hessian(3*natom,3*natom)
	real*8	sqrtmass(*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
	pointer	(ptr_kval,kval)
c
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
	common /elastic/	ptr_kval
c
c===================================================================================
c	non-local interactions
c===================================================================================
c
	call set_zeros(9*natom*natom, hessian)
c
	do 400 ic = 1,ncont
c
		i = listcontact(1,ic)
		j = listcontact(2,ic)
		Kconst = kval(ic)
c
		do 100 k = 1,3
			diff(k) = coord(3*(i-1)+k)-coord(3*(j-1)+k)
100		continue
c
		call dotvect(diff,diff,dist)
		fact = -Kconst/dist
c
		do 300 k = 1,3
			f1 = fact*diff(k)
			do 200 l = 1,3
				hessian(3*(i-1)+k,3*(j-1)+l)=
     1				hessian(3*(i-1)+k,3*(j-1)+l)
     2					+ f1*diff(l)
200			continue
300		continue
c
400	continue

c
c===================================================================================
c	Make symmetric, and fill in diagonal
c===================================================================================
c
	do 600 i = 2,3*natom
		do 500 j = 1,i-1
			hessian(i,j) = hessian(j,i)
500		continue
600	continue
c
	do 900 i = 1,3*natom
c
		l = i/3 + 1
		if(mod(i,3).eq.0) l = l-1
		do 800 j = 1,3
			sum = 0.d0
			do 700 k = j, 3*natom, 3
				sum = sum + hessian(i,k)
700			continue
			hessian(i,3*(l-1)+j) = -sum
800		continue
900	continue
c
c===================================================================================
c===================================================================================
c
	do 1100 i = 1, 3*natom
		k = i/3 + 1
		if(mod(i,3).eq.0) k = k-1
		do 1000 j = 1, 3*natom
			l = j/3 + 1
			if(mod(j,3).eq.0) l = l-1
			fact = 1.d0/(sqrtmass(corresp(k))*
     1				sqrtmass(corresp(l)))
			hessian(i,j) = hessian(i,j)*fact
1000		continue
1100	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	fill_hessian_elastic.f
c
c	This subroutine computes the Hessian for an Elastic model
c
c	Copyright (C) 2016 Patrice Koehl
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
	subroutine Fill_hessian_elastic(natom,corresp,
     1		sqrtmass,coord,spr_jdx,spr_row,spr_hessian)
c
c===================================================================================
c===================================================================================
c
	integer	i, j, k, l, iat, jat, irow
	integer	ipos,add1,add2
	integer	natom,i_start,i_k,j_k
	integer	nsize,ncont,npos
        integer spr_row(*),spr_jdx(*)
	integer ncontact(*),listcontact(2,*)	
	integer	corresp(*)
c
	real*8	dist, sum, Kconst
	real*8	fact, facti, facti2, factj
	real*8	f1
	real*8	diff(3)
	real*8	coord(*), kval(*), spr_hessian(*), sqrtmass(*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
	pointer	(ptr_kval,kval)
c
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
	common /elastic/	ptr_kval
c
c===================================================================================
c	non-local interactions
c===================================================================================
c
	npos = 0
c
	do 500 i = 1, natom
c
	   do 400 jat = 1,ncontact(i)
c
		j = listcontact(2,npos+jat)
		Kconst = kval(npos+jat)
c
		do 100 k = 1,3
			diff(k) = coord(3*(i-1)+k)-coord(3*(j-1)+k)
100		continue
c
		call dotvect(diff,diff,dist)
		fact = -Kconst/dist
c
		do 300 k = 1,3
			i_k = spr_row(3*(i-1)+k) + 3-k + 3*(jat-1)
			f1 = fact*diff(k)
			do 200 l = 1,3
				spr_hessian(i_k+l) = 
     1				spr_hessian(i_k+l)+f1*diff(l)
200			continue
300		continue
c
400	   continue
c
	   npos = npos + ncontact(i)
c
500	continue
c
c===================================================================================
c	Fill in diagonal blocks
c===================================================================================
c
	do 900 i = 1,natom
	    do 800 j = 1,3
c
		irow = 3*(i-1)+j
		i_k = spr_row(irow)
		i_start = i_k + 3-j 
c
		do 700 k = 1, 3
c
			sum = 0.d0
			do 600 l = i_start+k,spr_row(irow+1)-1,3
				f1 = spr_hessian(l)
				sum = sum + f1
				j_k = spr_jdx(l)
				if(j.ge.k) then
					j_k = spr_row(j_k) + j - k
					spr_hessian(j_k) = 
     1					spr_hessian(j_k)-f1
				endif
600			continue
			if(k.ge.j) spr_hessian(i_k+k-j) = 
     1				spr_hessian(i_k+k-j)-sum
c
700		continue
c
800	    continue
900	continue
c
c===================================================================================
c	Add mass information
c===================================================================================
c
	do 1200 iat = 1,natom
c
		facti = 1/sqrtmass(corresp(iat))
		facti2 = facti*facti
c
		ipos = spr_row(3*iat-2)
		add1 = spr_row(3*iat-1) - spr_row(3*iat-2) - 1
		add2 = spr_row(3*iat) - spr_row(3*iat-1) - 1
c
		spr_hessian(ipos) = spr_hessian(ipos)*facti2
		spr_hessian(ipos+1) = spr_hessian(ipos+1)*facti2
		spr_hessian(ipos+2) = spr_hessian(ipos+2)*facti2
		spr_hessian(ipos+add1+1) = 
     1			spr_hessian(ipos+add1+1)*facti2
		spr_hessian(ipos+add1+2) = 
     1			spr_hessian(ipos+add1+2)*facti2
		spr_hessian(ipos+add1+add2+2)=
     1			spr_hessian(ipos+add1+add2+2)*facti2
c
		do 1100 j = spr_row(3*iat-2) + 3,spr_row(3*iat-1)-1,3
c
			jat = (spr_jdx(j)+2)/3
	
			factj = facti/sqrtmass(corresp(jat))
c
			do 1000 l = 1,3
				spr_hessian(j+l-1) = spr_hessian(j+l-1)
     1					*factj
				spr_hessian(j+add1+l-1)=
     1				    spr_hessian(j+add1+l-1)*factj
				spr_hessian(j+add1+add2+l-1)=
     1				    spr_hessian(j+add1+add2+l-1)*factj
1000			continue
c
1100		continue
c
1200	continue
c
c===================================================================================
c===================================================================================
c
	return
	end
c===================================================================================
c===================================================================================
c	update_hessian_elastic.f
c
c	This subroutine updates the Hessian for an elastic model, when
c	only one elastic constant has been modified
c
c	Copyright (C) 2016 Patrice Koehl
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
	subroutine update_hessian_elastic(nchange,k_old,corresp,
     1		sqrtmass,coord,spr_jdx,spr_row,spr_hessian)
c
c===================================================================================
c===================================================================================
c
	integer	i, j, k, l, jat, irow
	integer	i_start,i_k,j_k
	integer	nsize,ncont,npos
	integer	nchange
	integer	iatom,jatom,icontact
        integer spr_row(*),spr_jdx(*)
	integer ncontact(*),listcontact(2,*)	
	integer	corresp(*)
c
	real*8	dist, Kconst, val
	real*8	fact, fact_old, factij, facti, factj
	real*8	f1, f2, k_old
	real*8	diff(3), H_old(3,3), H_new(3,3)
	real*8	coord(*), kval(*), spr_hessian(*), sqrtmass(*)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
	pointer	(ptr_kval,kval)
c
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
	common /elastic/	ptr_kval
c
c===================================================================================
c	Modify non-local interactions corresponding to "nchange"
c===================================================================================
c
	write(6,*) 'nchange =',nchange
	iatom = listcontact(1,nchange)
	jatom = listcontact(2,nchange)
	write(6,*) 'iatom,jatom:',iatom,jatom
c
	npos = 0
	do 100 i = 1,iatom-1
		npos = npos + ncontact(i)
100	continue
c
	do 200 jat = 1,ncontact(iatom)
		if(listcontact(2,npos+jat).eq.jatom) then
			icontact = jat
			goto 300
		endif
200	continue
300	continue
c	
	Kconst = kval(nchange)
c
	do 400 k = 1,3
		diff(k) = coord(3*(iatom-1)+k)-coord(3*(jatom-1)+k)
400	continue
c
	call dotvect(diff,diff,dist)
	fact = -Kconst/dist
	fact_old = -k_old/dist
	factij=1.d0/(sqrtmass(corresp(iatom))*sqrtmass(corresp(jatom)))
c
	do 600 k = 1,3
		i_k = spr_row(3*(iatom-1)+k) + 3-k + 3*(icontact-1)
		f1 = fact*diff(k)
		f2 = fact_old*diff(k)
		do 500 l = 1,3
			H_new(k,l) = f1*diff(l)
			H_old(k,l) = f2*diff(l)
			spr_hessian(i_k+l) = f1*diff(l)*factij
500		continue
600	continue
c
c===================================================================================
c	Modify diagonal blocks for iatom, jatom
c===================================================================================
c
	facti = 1/sqrtmass(corresp(iatom))
	facti = facti*facti
	factj = 1/sqrtmass(corresp(jatom))
	factj = factj*factj
c
	do 800 j = 1,3
c
		irow = 3*(iatom-1)+j
		i_k = spr_row(irow)
		i_start = i_k + 3-j +3*(icontact-1)
c
		do 700 k = 1, 3
c
			val = H_new(j,k)-H_old(j,k)
			l = i_start+k
			j_k = spr_jdx(l)
			if(j.ge.k) then
				j_k = spr_row(j_k) + j - k
				spr_hessian(j_k) = 
     1				spr_hessian(j_k)-val*factj
			endif
			if(k.ge.j) then
				spr_hessian(i_k+k-j) = 
     1				spr_hessian(i_k+k-j)-val*facti
			endif
c
700		continue
c
800	continue
c
c===================================================================================
c===================================================================================
c
	return
	end
