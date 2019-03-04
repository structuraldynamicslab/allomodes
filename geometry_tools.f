c ==========================================================================================
c ==========================================================================================
c	This file contains a library of fortran routines used to
c	perform simple vector and geometric operations
c
c	Copyright (C) 1990 Patrice Koehl
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
c ==========================================================================================
c ==========================================================================================
c	Vector		Version 1 8/2/1990		Patrice Koehl
c
c	This file contains several subroutine that can be used for any
c	vector operations (vector in 3D cartesian space)
c
c	This includes :	crossvect	: cross vector of two vectors
c			dotvect		: dot product of two vectors
c			normvect	: norm of a vector
c			detvect		: determinant of three vectors
c			diffvect	: substract two vectors
c			addvect		: add two vectors
c
c	For each subroutine : u is for vectors (arrays of size 3)
c			      all other value are scalar
c			      calculations are done in double precision
c ==========================================================================================
c ==========================================================================================
c
c	1 . crossvect :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine crossvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	u3(1) = u1(2)*u2(3) - u1(3)*u2(2)
	u3(2) = -u1(1)*u2(3) + u1(3)*u2(1)
	u3(3) = u1(1)*u2(2) - u1(2)*u2(1)
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
c	2. dotvect :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine dotvect(u1,u2,dot)
c
	integer	i
c
	real*8	u1(3),u2(3),dot
c
	dot = 0.d0
	do 100 i = 1,3
		dot = dot + u1(i)*u2(i)
100	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c
c	3. normvect :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine normvect(u1,norm)
c
	real*8	u1(3),norm
c
	call dotvect(u1,u1,norm)
	norm = dsqrt(norm)
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c
c	4. detvect :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine detvect(u1,u2,u3,det)
c
	real*8	u1(3),u2(3),u3(3),det,u4(3)
c
	call crossvect(u2,u3,u4)
	call dotvect(u1,u4,det)
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c
c	5. diffvect :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine diffvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u2(i) - u1(i)
100	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c
c	6. addvect :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine addvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u1(i) + u2(i)
100	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c
c	7. Normalise a vector : given a vector u1, output u1/norm(u1) :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine unitvector(u1,u2)
c
	real*8  u1(3),u2(3),norm
c
	integer i
c
	call normvect(u1,norm)
c
	do 100 i = 1,3
		u2(i) = u1(i)/norm
100	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c
c	8. Normalise a vector : given a vector u1, output norm(u1) and u1/norm(u1) :
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine unitvect(u1,u2,norm)
c
	real*8  u1(3),u2(3),norm
c
	integer i
c
	call normvect(u1,norm)
c
	do 100 i = 1,3
		u2(i) = u1(i)/norm
100	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	nsp3.f		version 1 8/3/1990		Patrice Koehl
c
c	This subroutine adds a proton H4 on a carbon with a sp3
c	conformation
c
c				X
c				!
c			C1 ---- C2 ---- C3
c				!
c				H4
c
c	Input of the program :
c				ang1	: angle C1-C2-H4 
c				ang2	: angle C2-C3-H4 
c				ip	: 1 for H4, 2 for X
c				dist	: length of C2-H4 (or C2-X)
c				p1	: C1
c				p2	: C2
c				p3	: C3
c	Output of the program :
c				p4	: H4 (or X)
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine nsp3(p1,p2,p3,p4,ang1,ang2,ip,dist)
c
	real*8	p1(3),p2(3),p3(3),p4(3),ang1,ang2,dist
	real*8	v1(3),v2(3),v3(3),b(3),u1(3),dd(3),det
	real*8	dist1,dist2,d21(3),d23(3),pn(3),pm(3)
	real*8	alpha,val
c
	integer i,ip
c
	call diffvect(p2,p1,d21)
	call diffvect(p2,p3,d23)
	call normvect(d21,dist1)
	call normvect(d23,dist2)
	call crossvect(d21,d23,u1)
c
	do 10 i = 1,3
		v1(i) = dist*cos(ang1)*d21(i)/dist1
		v2(i) = dist*cos(ang2)*d23(i)/dist2
10	continue
c
	call addvect(v1,p2,pn)
	call addvect(v2,p2,pm)
c
	call dotvect(pn,d21,b(1))
	call dotvect(pm,d23,b(2))
	call dotvect(p2,u1,b(3))
c
	v1(1) = d21(1)
	v2(1) = d21(2)
	v3(1) = d21(3)
c
	v1(2) = d23(1)
	v2(2) = d23(2)
	v3(2) = d23(3)
c
	v1(3) = u1(1)
	v2(3) = u1(2)
	v3(3) = u1(3)
c
	call detvect(v1,v2,v3,det)
c
	call detvect(b,v2,v3,dd(1))
	call detvect(v1,b,v3,dd(2))
	call detvect(v1,v2,b,dd(3))
c
	do 100 i = 1,3
		v2(i) = dd(i)/det
100	continue
c
	call diffvect(p2,v2,v1)
	call normvect(u1,dist1)
	call normvect(v1,dist2)
c
	val = dabs(dist**2-dist2**2)
	if(ip.eq.1) then
		alpha = -dsqrt(val)
	else
		alpha = dsqrt(val)
	endif
c
	do 200 i = 1,3
		u1(i) = (alpha/dist1) * u1(i)
200	continue
c
	call addvect(u1,v2,p4)
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	sp2.f		version 1 8/3/1990		Patrice Koehl
c
c	This subroutine adds a proton H4 on a carbon with a sp2
c	conformation
c
c		           C1		 
c			      \		 
c			        C2 ---- C3
c			      /	
c			   H4
c
c	Input of the program :
c				ang1	: angle H4-c2-c1
c				dist	: length of C2-H4 
c				p1	: C1
c				p2	: C2
c				p3	: C3
c	Output of the program :
c				p4	: H4 (or X)
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine sp2(p1,p2,p3,p4,ang1,dist)
c
	real*8	p1(3),p2(3),p3(3),p4(3),ang1,ang2,dist
	real*8	v1(3),v2(3),v3(3),b(3),u1(3),dd(3),det
	real*8	dist1,dist2,d12(3),d23(3),dot,cos1,ang3
c
	integer i
c
	include 'const_pi.h'
c
	call diffvect(p1,p2,d12)
	call diffvect(p2,p3,d23)
c
	call normvect(d12,dist1)
	call normvect(d23,dist2)
c
	call dotvect(d12,d23,dot)
	cos1 = -dot/(dist1*dist2)
	ang3 = acos(cos1)
	ang2 = 2*pi - ang1 - ang3
c
	b(1) = - dist1*dist*cos(ang1)
	b(2) = dist2*dist*cos(ang2)
	b(3) = 0.d0
c
	v1(1) = d12(1)
	v2(1) = d12(2)
	v3(1) = d12(3)
c
	v1(2) = d23(1)
	v2(2) = d23(2)
	v3(2) = d23(3)
c
	call crossvect(d12,d23,u1)
c
	v1(3) = u1(1)
	v2(3) = u1(2)
	v3(3) = u1(3)
c
	call detvect(v1,v2,v3,det)
c
	call detvect(b,v2,v3,dd(1))
	call detvect(v1,b,v3,dd(2))
	call detvect(v1,v2,b,dd(3))
c
	do 100 i = 1,3 
		u1(i) = dd(i)/det
100	continue
c
	call addvect(u1,p2,p4)
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Buildfour.f		Version 1 8/2/1990	Patrice Koehl
c
c	This subroutine, knowing the position of 3 points and the torsion
c	angle that relates them to a 4th point, get the coordinates of
c	this fourth point
c	For details on how it works, see author
c	input of the program :
c			p1,p2,p3 and p4 are the four points considered;
c			tor is the torsional angle (in radian)
c			dist is the distance between p3 and p4
c			ang is ang(p2p3,p3p4)
c	output :
c			p4, coordinates of the fourth points
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine buildfour(p1,p2,p3,p4,tor,ang,dist)
c
	integer i
	real*8	p1(3),p2(3),p3(3),p4(3),tor,dist,ang
	real*8	d21(3),d23(3),u1(3),u2(3)
	real*8	fact,b(3),v1(3),v2(3),v3(3),fact2
	real*8	norm,norm2,det,detx(3),sin1
	real*8	norm3,dot,cosm1
c
	call diffvect(p1,p2,d21)
	call diffvect(p2,p3,d23)
c
	call crossvect(d21,d23,u1)
	call dotvect(d21,d23,dot)
c
	call normvect(u1,norm)
	call normvect(d21,norm2)
	call normvect(d23,norm3)
c
	cosm1 = acos(dot/(norm2*norm3))
	sin1 = dabs(sin(cosm1))
c
	v1(1) = u1(1)/norm
	v2(1) = u1(2)/norm
	v3(1) = u1(3)/norm
c
	v1(2) = d21(1)/norm2
	v2(2) = d21(2)/norm2
	v3(2) = d21(3)/norm2
c
	v1(3) = d23(1)
	v2(3) = d23(2)
	v3(3) = d23(3)
c
	b(1)  = cos(tor)
	b(2)  = sin(tor)*sin1
	b(3)  = 0.d0
c
	call detvect(v1,v2,v3,det)
c
	call detvect(b,v2,v3,detx(1))
	call detvect(v1,b,v3,detx(2))
	call detvect(v1,v2,b,detx(3))
c
	do 100 i = 1,3
		u2(i) = detx(i)/det
100	continue
c
	fact = dist*dabs(sin(ang))/norm3
	fact2 = dist*dabs(cos(ang))/norm3
c
	call crossvect(u2,d23,v1)
c
	do 200 i = 1,3
		v1(i) = v1(i) * fact
		v2(i) = d23(i)* fact2
200	continue
c
	call addvect(v1,v2,v3)
	call addvect(v3,p3,p4)
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Torsion.f	Version 1 8/02/90	Patrice Koehl
c
c	This program calculates the dihedral (or torsion) angle defined
c	by 4 points
c	Input		- p1,p2,p3,p4	: 4 points defining the angle
c	Output		- val		: angle value (in radian)
c ==========================================================================================
c ==========================================================================================
c
	subroutine torsion(p1,p2,p3,p4,val)
c
	real*8 	val,dot,alpha1,alpha2,det,cosi,eps
	real*8	p1(3),p2(3),p3(3),p4(3)
	real*8	d21(3),d23(3),d34(3),u1(3),u2(3)
c
	eps = 1.d-10
c
c ==========================================================================================
C	Now calculate the dihedral angle cosine using the expression :
C		cosine(p1,p2,p3,p4) = u1*u2/!u1!*!u2!
C	where :
C			u1 = vect(p2,p1)^vect(p2,p3)
C			u2 = vect(p2,p3)^vect(p3,p4)
C	! ! represents the norm of the vector, and * the dot product
C	and 
C		angle = acos(cosine)
C	The sign of the angle is found from the sine expression
C		sine = det(u3,u1,u2)/!u1!*!u2!*!u3!
C	where u3 = vect(p2,p3)
c ==========================================================================================
C
c ==========================================================================================
C	First get u1, and u2 :
c ==========================================================================================
C
	call diffvect(p1,p2,d21)
	call diffvect(p2,p3,d23)
	call diffvect(p3,p4,d34)
c
	call crossvect(d21,d23,u1)
	call crossvect(d23,d34,u2)
c
	call dotvect(u1,u2,dot)
c
	call normvect(u1,alpha1)
	call normvect(u2,alpha2)
c
c ==========================================================================================
c	then calculate the cosine value :
c ==========================================================================================
c
	cosi = dot/(alpha1*alpha2)
c
	if(abs(cosi-1).le.eps) cosi = 1.d0-eps
	if(abs(cosi+1).le.eps) cosi = -1.d0+eps
c
c ==========================================================================================
c	now get the sine to find sign of the angle :
c ==========================================================================================
c
	call detvect(d23,u1,u2,det)
c
	if(det.le.0) then
		val = -dacos(cosi)
	else
		val = dacos(cosi)
	endif
c
	return    
	end
c ==========================================================================================
c ==========================================================================================
c	Complete_peptbond.f	Version 1 2/29/2000	Patrice Koehl
c
c	This subroutine builds one missing atom of the four atoms
c	forming a peptide bond, based on the knowledge of the 3 others
c ==========================================================================================
c ==========================================================================================
c
	subroutine complete_peptbond(a,b,c,d,imiss)
c
c ==========================================================================================
c	Input:
c		a,b,c,d		: the four atoms forming the peptide
c				  bond (one of them has coordinates
c				  0,0,0)
c		imiss		: missing atom number (from 1 to 4)
c	Output:
c		coordinates of the missing atom
c ==========================================================================================
c
	integer	i,ichange
	integer	imiss
c
	real*8	pi,d2r
	real*8	dist1,dist2,dist3
	real*8	tor,ang1,ang2,a_val
	real*8	ab_bd,ab_bc,bc_bd,cb_ac,cb_cd,ac_cd
	real*8	dist_ac2,dist_bd2,dist_ab2,dist_cd2,dist_ad2
	real*8	Acoef,Bcoef,det
	real*8	a(3),b(3),c(3),d(3)
	real*8	vect1(3),vect2(3),vect3(3),vect_cb(3),vect_bc(3)
c
	pi  = dacos(-1.d0)
	d2r = pi/180.d0
c
	ang1 = 114.d0*d2r
	ang2 = 123.d0*d2r
c
	dist1 = 1.53d0
	dist2 = 1.32d0
	dist3 = 1.47d0
c
	if(imiss.eq.1) then
c
c 		============================================================================
c		1. Imiss = 1: atom "a" (i.e. Ca of first residue) 
c		   is missing
c		It is built from the three other atoms, based on the 
c		fact that the torsion angle (a,b,c,d) = 0
c 		============================================================================
c
		tor  = pi
		call buildfour(d,c,b,a,tor,ang1,dist1)
c
	elseif(imiss.eq.2) then
c
c 		============================================================================
c		2. Imiss = 2: atom "b" (i.e. C of first residue) 
c		is missing
c		It is built by saying that vector(b,c) is a linear
c		combination of vector(a,c) and vector(c,d). 
c		The coefficients of this linear combination are 
c		derived by using scalar	product:
c		vector(c,b) = A vector(a,c) + B vector(c,d)
c
c		vector(c,b)*vector(a,c)=A d(a,c)**2 + 
c				B vector(c,d)*vector(a,c)
c		and
c		vector(c,b)*vector(c,d)=A vector(a,c)*vector(c,d) 
c				+ B d(c,d)**2
c
c		from which A and B are computed, and then used to get 
c		coordinates of b
c 		============================================================================
c
		call diffvect(a,c,vect1)
		call dotvect(vect1,vect1,dist_ac2)
c
		call diffvect(c,d,vect2)
		call dotvect(vect1,vect2,ac_cd)
c
		dist_cd2 = dist3*dist3
c
		cb_cd = dist2*dist3*cos(ang2)
c
		cb_ac = - (dist2**2 + dist_ac2 - dist1**2)/2.d0
c
		det = dist_ac2*dist_cd2 - ac_cd*ac_cd
c
		Acoef = (cb_ac*dist_cd2 - cb_cd*ac_cd)/det
		Bcoef = (dist_ac2*cb_cd - ac_cd*cb_ac)/det
c
		do 100 i = 1,3
			vect_cb(i) = Acoef*vect1(i) + Bcoef*vect2(i)
100		continue
c
		call addvect(c,vect_cb,b)
c
		call diffvect(a,d,vect1)
		call dotvect(vect1,vect1,dist_ad2)
c
		call torsion(a,b,c,d,tor)
c
		if(dist_ad2.lt.10) then
			if(abs(tor).gt.170*d2r) then
				ichange = 1
			else
				ichange = 0
			endif
		else
			if(abs(tor).le.10*d2r) then
				ichange = 1
			else
				ichange = 0
			endif
		endif
		if(ichange.eq.1) then
			call unitvector(vect1,vect2)
			call diffvect(a,b,vect1)
			call dotvect(vect1,vect2,a_val)
			a_val = 2*a_val
			do 150 i = 1,3
				vect3(i) = a_val*vect2(i)-vect1(i)
150			continue
			call addvect(a,vect3,b)
		endif
c
	elseif(imiss.eq.3) then
c
c 		============================================================================
c		3. Imiss = 3: atom "c" (i.e. N of last residue) 
c			is missing
c		It is built by saying that vector(b,c) is a linear
c		combination of vector(a,b) and vector(b,d). 
c		The coefficients of this linear combination are 
c		derived by using scalar	product:
c
c		vector(b,c) = A vector(a,b) + B vector(b,d)
c
c		vector(b,c)*vector(a,b)=A d(a,b)**2 
c			+ B vector(b,d)*vector(a,b)
c		and
c		vector(b,c)*vector(b,d)=A vector(a,b)*vector(b,d) 
c				+ B d(b,d)**2
c
c		from which A and B are computed, and then used to get 
c		coordinates of c
c 		============================================================================
c
		call diffvect(a,b,vect1)
		dist_ab2 = dist1*dist1
c
		call diffvect(b,d,vect2)
		call dotvect(vect1,vect2,ab_bd)
c
		call dotvect(vect2,vect2,dist_bd2)
c
		ab_bc = -dist1*dist2*cos(ang1)
c
		bc_bd = (dist2**2 + dist_bd2 - dist3**2)/2.d0
c
		det = dist_ab2*dist_bd2 - ab_bd*ab_bd
c
		Acoef = (ab_bc*dist_bd2 - bc_bd*ab_bd)/det
		Bcoef = (dist_ab2*bc_bd - ab_bc*ab_bd)/det
c
		do 200 i = 1,3
			vect_bc(i) = Acoef*vect1(i) + Bcoef*vect2(i)
200		continue
c
		call addvect(b,vect_bc,c)
		call diffvect(a,d,vect1)
		call dotvect(vect1,vect1,dist_ad2)
c
		call torsion(a,b,c,d,tor)
c
		if(dist_ad2.lt.10) then
			if(abs(tor).gt.170*d2r) then
				ichange = 1
			else
				ichange = 0
			endif
		else
			if(abs(tor).le.10*d2r) then
				ichange = 1
			else
				ichange = 0
			endif
		endif
		if(ichange.eq.1) then
			call unitvector(vect1,vect2)
			call diffvect(d,c,vect1)
			call dotvect(vect1,vect2,a_val)
			a_val = 2*a_val
			do 300 i = 1,3
				vect3(i) = a_val*vect2(i)-vect1(i)
300			continue
			call addvect(d,vect3,c)
		endif
			
c		if(abs(tor).le.10*d2r) then
c			call diffvect(c,d,vect1)
c 			call addvect(b,vect1,c)
c		endif
c
	elseif(imiss.eq.4) then
c
c 		============================================================================
c		4. Imiss = 4: atom "d" (i.e. Ca of last residue) 
c		is missing
c		It is built from the three other atoms, based on the 
c		fact that
c		the torsion angle (a,b,c,d) = pi
c 		============================================================================
c
		tor  = pi
		call buildfour(a,b,c,d,tor,ang2,dist3)
c
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	rotate_torsion.f	Version 1 8/02/90	Patrice Koehl
c
c	This program rotates a point around an axis given by a point and a unit vector
c	Input		- e,r           : unit vector and point representing axis
c			- ang		: angle to be applied (assumed in radian)
c			- p1		: initial position of the point
c	Output		- p2		: rotated point
c ==========================================================================================
c ==========================================================================================
c
	subroutine rotate_torsion(e,r,ang,p1,p2)
c
	integer	i
c
	real*8	ang,cos_ang,sin_ang
	real*8	val1
	real*8	e(3),p1(3),p2(3),r(3)
	real*8	v1(3),v2(3)
c
	cos_ang = cos(ang)
	sin_ang = sin(ang)
	call diffvect(r,p1,v1)
	call dotvect(v1,e,val1)
	val1 = val1*(1.d0-cos_ang)
c
	call crossvect(e,v1,v2)
c
	do 100 i = 1,3
		p2(i) = r(i) + val1*e(i)+cos_ang*v1(i)+sin_ang*v2(i)
100	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	Inverse3x3.f
c
c	This subroutine computes the inverse of a 3x3 matrix
c===================================================================================
c===================================================================================
c
	subroutine inverse3x3(mat,inv)
c
	integer	i,j
c
	real*8	det,detinv,deter
	real*8	mat(3,3)
	real*8	inv(3,3)
c
c===================================================================================
c	Compute determinant
c===================================================================================
c
	det = deter(mat)
c
	if(det.eq.0) then
		write(6,*) 'Problem with one inertia matrix!!'
		stop
	endif
c
	detinv = 1.d0/det
c
	inv(1,1) = mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3)
	inv(1,2) = -(mat(1,2)*mat(3,3)-mat(1,3)*mat(3,2))
	inv(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
	inv(2,1) = -(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))
	inv(2,2) = mat(1,1)*mat(3,3)-mat(3,1)*mat(1,3)
	inv(2,3) = -(mat(1,1)*mat(2,3)-mat(1,3)*mat(2,1))
	inv(3,1) = mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)
	inv(3,2) = -(mat(1,1)*mat(3,2)-mat(1,2)*mat(3,1))
	inv(3,3) = mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2)
c
	do 200 j = 1,3
		do 100 i = 1,3
			inv(i,j)=inv(i,j)*detinv
100		continue
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	MatVect3x3.f
c
c	This subroutine computes the product of a 3x3 matrix by a vector
c===================================================================================
c===================================================================================
c
	subroutine MatVect3x3(mat,v1,v2)
c
	integer	i,j
c
	real*8	mat(3,3)
	real*8	v1(3),v2(3)
c
c===================================================================================
c	Compute determinant
c===================================================================================
c
	do 200 i = 1,3
		v2(i) = 0.d0
		do 100 j = 1,3
			v2(i) = v2(i) + mat(i,j)*v1(j)
100		continue
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	set_zeros.f
c
c	Initializes a vector to 0
c
c===================================================================================
c===================================================================================
c
	subroutine set_zeros(n,vect)
c
	integer	i,n
c
	real*8	vect(*)
c
	do 100 i = 1,n
		vect(i) = 0.d0
100	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	scale.f
c
c	Scale a vector
c
c===================================================================================
c===================================================================================
c
	subroutine scale(n,vect,val)
c
	integer	i,n
c
	real*8	val,vect(*)
c
	do 100 i = 1,n
		vect(i) = val*vect(i)
100	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	copy_vect.f
c
c	Copy a vector into another
c
c===================================================================================
c===================================================================================
c
	subroutine copy_vect(n,vect1,vect2)
c
	integer	i,n
c
	real*8	vect1(*),vect2(*)
c
	do 100 i = 1,n
		vect2(i) = vect1(i)
100	continue
c
	return
	end
c===================================================================================
c===================================================================================
c	scale_vect.f
c
c	scale a vector with a constant
c
c===================================================================================
c===================================================================================
c
	subroutine scale_vect(n,alpha,vect)
c
	integer	i,n
c
	real*8	alpha
	real*8	vect(*)
c
	do 100 i = 1,n
		vect(i) = alpha*vect(i)
100	continue
c
	return
	end
