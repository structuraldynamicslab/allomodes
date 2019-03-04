c=================================================================================
c=================================================================================
c	sorting_tools.f
c
c	This file contains a library of fortran routines used to
c	sort numbers using heapsort
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
c	Hpsort_key.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in ascending order, and
c	provide an index of the ranked element
c=================================================================================
c=================================================================================
c
	subroutine hpsort_key(ra,index,n)
c
	integer	n,i,ir,j,l,idx
	integer	index(n)
c
	real*8	rra
	real*8	ra(n)
c
	save
c
	do 50 i = 1,n
		index(i) = i
50	continue
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		rra = ra(l)
		idx = l
	else
		rra = ra(ir)
		idx = index(ir)
		ra(ir) = ra(1)
		index(ir) = index(1)
		ir = ir -1
		if(ir.eq.1) then
			ra(1) = rra
			index(1) = idx
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(ra(j).lt.ra(j+1)) j = j + 1
		endif
		if(rra.lt.ra(j)) then
			ra(i) = ra(j)
			index(i) = index(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	ra(i) = rra
	index(i) = idx
	goto 100
c
	end
c=================================================================================
c=================================================================================
c	Hpsort_int_key.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in ascending order, and
c	provide an index of the ranked element
c=================================================================================
c=================================================================================
c
	subroutine hpsort_int_key(ra,index,n)
c
	integer	n,i,ir,j,l,idx
	integer	index(n)
c
	integer	rra
	integer ra(n)
c
	save
c
	do 50 i = 1,n
		index(i) = i
50	continue
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		rra = ra(l)
		idx = l
	else
		rra = ra(ir)
		idx = index(ir)
		ra(ir) = ra(1)
		index(ir) = index(1)
		ir = ir -1
		if(ir.eq.1) then
			ra(1) = rra
			index(1) = idx
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(ra(j).lt.ra(j+1)) j = j + 1
		endif
		if(rra.lt.ra(j)) then
			ra(i) = ra(j)
			index(i) = index(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	ra(i) = rra
	index(i) = idx
	goto 100
c
	end
c=================================================================================
c=================================================================================
c	Hpsort_three.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in ascending order, and
c	provide an index of the ranked element
c=================================================================================
c=================================================================================
c
	subroutine hpsort_three(ra,index,n)
c
	integer	n,i,ir,j,l,idx,comp3,k
	integer	index(n)
c
	save
c
	real*8	rra(3)
	real*8	ra(3,n)
c
	do 50 i = 1,n
		index(i) = i
50	continue
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		do 10 k = 1,3
			rra(k) = ra(k,l)
10		continue
		idx = l
	else
		do 20 k = 1,3
			rra(k) = ra(k,ir)
20		continue
		idx = index(ir)
		do 30 k = 1,3
			ra(k,ir) = ra(k,1)
30		continue
		index(ir) = index(1)
		ir = ir -1
		if(ir.eq.1) then
			do 40 k = 1,3
				ra(k,1) = rra(k)
40			continue
			index(1) = idx
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(comp3(ra(1,j),ra(1,j+1)).eq.1) j = j + 1
		endif
		if(comp3(rra,ra(1,j)).eq.1) then
			do 210 k = 1,3
				ra(k,i) = ra(k,j)
210			continue
			index(i) = index(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	do 220 k = 1,3
		ra(k,i) = rra(k)
220	continue
	index(i) = idx
	goto 100
c
	end
c
c=================================================================================
c=================================================================================
c
c	Comp3.f		version 1 10/17/2006	Patrice Koehl
c
c	This function compares 2 arrays of 3 real numbers
c
c=================================================================================
c=================================================================================
c
	function comp3(a,b)
c
	integer	comp3,k
	real*8	a(3),b(3)
c
	save
c
	comp3 = 0
c
	do 100 k = 1,3
		if(a(k).lt.b(k)) then
			comp3 = 1
			return
		elseif(a(k).gt.b(k)) then
			return
		endif
100	continue
c
	return
	end
c=================================================================================
c=================================================================================
c	Hpsort_two_int.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in ascending order, and
c	provide an index of the ranked element
c
c	Elements of the array are integers
c=================================================================================
c=================================================================================
c
	subroutine hpsort_two_int(ra,n)
c
	integer	n,i,ir,j,l,comp2,k
c
	save
c
	integer	rra(2)
	integer	ra(2,n)
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		do 10 k = 1,2
			rra(k) = ra(k,l)
10		continue
	else
		do 20 k = 1,2
			rra(k) = ra(k,ir)
20		continue
		do 30 k = 1,2
			ra(k,ir) = ra(k,1)
30		continue
		ir = ir -1
		if(ir.eq.1) then
			do 40 k = 1,2
				ra(k,1) = rra(k)
40			continue
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(comp2(ra(1,j),ra(1,j+1)).eq.1) j = j + 1
		endif
		if(comp2(rra,ra(1,j)).eq.1) then
			do 210 k = 1,2
				ra(k,i) = ra(k,j)
210			continue
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	do 220 k = 1,2
		ra(k,i) = rra(k)
220	continue
	goto 100
c
	end
c=================================================================================
c=================================================================================
c	Hpsort_two_int_idx.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in ascending order, and
c	provide an index of the ranked element
c
c	Elements of the array are integers
c=================================================================================
c=================================================================================
c
	subroutine hpsort_two_int_idx(ra,index,n)
c
	integer	n,i,ir,j,l,idx,comp2,k
	integer	index(n)
c
	save
c
	integer	rra(2)
	integer	ra(2,n)
c
	do 50 i = 1,n
		index(i) = i
50	continue
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		do 10 k = 1,2
			rra(k) = ra(k,l)
10		continue
		idx = l
	else
		do 20 k = 1,2
			rra(k) = ra(k,ir)
20		continue
		idx = index(ir)
		do 30 k = 1,2
			ra(k,ir) = ra(k,1)
30		continue
		index(ir) = index(1)
		ir = ir -1
		if(ir.eq.1) then
			do 40 k = 1,2
				ra(k,1) = rra(k)
40			continue
			index(1) = idx
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(comp2(ra(1,j),ra(1,j+1)).eq.1) j = j + 1
		endif
		if(comp2(rra,ra(1,j)).eq.1) then
			do 210 k = 1,2
				ra(k,i) = ra(k,j)
210			continue
			index(i) = index(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	do 220 k = 1,2
		ra(k,i) = rra(k)
220	continue
	index(i) = idx
	goto 100
c
	end
c
c=================================================================================
c=================================================================================
c
c	Comp2.f		version 1 10/17/2006	Patrice Koehl
c
c	This function compares 2 arrays of 3 real numbers
c
c=================================================================================
c=================================================================================
c
	function comp2(a,b)
c
	integer	comp2,k
	integer	a(2),b(2)
c
	save
c
	comp2 = 0
c
	do 100 k = 1,2
		if(a(k).lt.b(k)) then
			comp2 = 1
			return
		elseif(a(k).gt.b(k)) then
			return
		endif
100	continue
c
	return
	end
