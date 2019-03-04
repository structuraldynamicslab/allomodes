c===================================================================================
c===================================================================================
c
c	Normal_modes.f		Version 1 7/9/2010		Patrice Koehl
c
c	This code uses ARPACK to find a few eigenvalues 
c	(lambda) and corresponding eigenvectors (x) for the standard 
c	eigenvalue problem:
c          
c                        A*x = lambda*x
c 
c     	where A is an n by n real symmetric matrix.
c
c     The only thing that must be supplied in order to use this routine
c     is a function that performs a matrix-vector product
c
c                         w <-  Av
c
c	This specific implementation is designed to compute a few eigenvalues and the
c	corresponding eigenvectors of the Hessian of the energy ofa biomolecule, 
c	where the Hessian is stored using a sparse format
c
c===================================================================================
c	Subroutine called
c===================================================================================
c
c	From ARPACK:
c
c	dsaupd  ARPACK reverse communication interface routine.
c	dseupd  ARPACK routine that returns Ritz values and (optionally)
c             	Ritz vectors.
c	dnrm2   Level 1 BLAS that computes the norm of a vector.
c	daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c
c	For matrix-vector multiplication:
c
c	sprsym_cr_mult : performs A*x, where A is a sparse, symmetric matrix
c			 stored in Condensed Row format	
c
c===================================================================================
c===================================================================================
c
	subroutine normal_modes(ndim,spr_jdx,spr_row,spr_hessian,
     1		m,nval,d,V,eigenval,eigenvect,work1,select)
c
c===================================================================================
c	Local Arrays
c===================================================================================
c
	integer		ndim,m,nval
	integer		nmvp
c
	real*8		V(ndim,m)
	real*8		workl(*),workd(*),work1(*)
	real*8		d(m,2)
c
	logical		select(*)
c
	integer		iparam(11),ipntr(11)
c
c===================================================================================
c	Local Scalars
c===================================================================================
c
	character	bmat*1,which*2
c
	integer		ido,n,nev,ncv,lworkl,info,ierr
	integer		k,j,nx,ishfts,maxitr,mode1,nconv
c
	logical		rvec
c
	real*8		tol,sigma
c
c===================================================================================
c	Parameters
c===================================================================================
c
	real*8		zero
	parameter	(zero=0.0d0)
c  
c===================================================================================
c	BLAS & LAPACK routines used
c===================================================================================
c
	real*8		dnrm2
	external	dnrm2,daxpy
c
c===================================================================================
c	Intrinsic function
c===================================================================================
c
	intrinsic 	abs
c
c===================================================================================
c	Input for Hessian matrix
c===================================================================================
c
	integer		spr_jdx(*),spr_row(*)
c
	real*8		spr_hessian(*)
c
c===================================================================================
c	Results of the calculation: eigenvalues, eigenvectors
c===================================================================================
c
	real*8		eigenval(*),eigenvect(ndim,*)
c
c===================================================================================
c	Pointers for variable size arrays
c===================================================================================
c
	pointer		(ptr_workl,workl)
	pointer		(ptr_workd,workd)
c
c===================================================================================
c	Common blocks
c===================================================================================
c
	common /space_lanczos/	ptr_workl,ptr_workd
c
c ============================================================================
c 	Set dimensions for this problem
c 	(see header in dsaupd for the meaning of each variable)
c ============================================================================
c
	nx = ndim
	n = nx
c
	nev = m/2
	ncv = m
c
	bmat  = 'I'
	which = 'SM'
c
	lworkl = ncv*(ncv+8)
	tol = zero 
	info = 0
	ido = 0
c
	ishfts = 1
	maxitr = 1000 
	mode1 = 1
c
	iparam(1) = ishfts
	iparam(3) = maxitr
	iparam(7) = mode1
c
c ============================================================================
c 	Main loop to get Lanczos decomposition
c ============================================================================
c
	nmvp = 0
100	continue
c
c               Call dsaupd repeatedly to build Krylov subspace
c               ===================================================
c
         	call dsaupd (ido,bmat,n,which,nev,tol,work1,ncv,V,
     1			ndim,iparam,ipntr,workd,workl,lworkl,info) 
c
c               Take action indicated by parameter IDO
c               ===================================================
c
		if(ido.eq.-1.or.ido.eq.1) then
c
c                       Perform matrix-vector multiplication
c                       ====================================
c
			call sprsym_cr_mult(spr_hessian,spr_row,
     1			spr_jdx,workd(ipntr(1)),workd(ipntr(2)),
     2			ndim)
			nmvp = nmvp + 1
c
			goto 100
c
		endif
c
c       We have either reached convergence, or there is an error
c       =========================================================
c
	if(info.lt.0) then
c
		print *, ' '
		print *, ' Error with dsaupd, info = ', info
		print *, ' Check documentation in dsaupd '
		print *, ' '
c
	else 
c
c       ===============================================================
c	No fatal error: proceed and compute eigenvalues (Ritz values)
c	and eigenvectors (Ritz vectors) of Hessian
c       ===============================================================
c
		rvec = .true.
c
		call dseupd(rvec,'All',select,d,V,ndim,sigma,
     1			bmat,n,which,nev,tol,work1,ncv,V,ndim,
     2			iparam,ipntr,workd,workl,lworkl,ierr)
c
		if(ierr.ne.0) then
c
c		Error condition: check the documentation of dseupd
c       	======================================================
c
			print *, ' '
			print *, ' Error with dseupd, info = ', ierr
             		print *, ' Check the documentation of _deupd. '
             		print *, ' '
c
		else
c
c       	======================================================================
c		Only nconv (=iparam(5)) values have converged; corresponding
c		eigenvalues are returned in the first column of the two dimensional
c		array d and the corresponding eigenvectors are returned in the first
c		nconv columns of the matrix V, if requested
c       	======================================================================
c
			nconv =  iparam(5)
			nval = min(nval,nconv)
c
			do 300 j = 1,nconv
c
				eigenval(j) = d(j,1)
				do 200 k = 1,nx
					eigenvect(k,j) = v(k,j)
200				continue
c
300			continue
c
		endif
c
	endif
c
c ============================================================================
c	Done with program
c ============================================================================
c
	write(6,*) ' '
	write(6,*) 'Number of matrix-vector products: ',nmvp
	write(6,*) ' '
c
	return
	end
c
c===================================================================================
c===================================================================================
c
c	sprsym_cr_mult.f
c
c	This subroutine computes the product of a real, symmetric
c	sparse matrix M, with a vector A.
c		(M A = B)
c
c	For efficiency, the matrix M is stored as a sparse matrix, using
c	the compressed row scheme
c	This storage scheme requires two vectors:
c		- spr_hessian:	contains the non zero values of the Hessian
c				matrix (only the upper part, as the
c				Hessian is symmetric)
c		- spr_row:	pointer of the position of the beginning of each row
c				in spr_hessian
c		- spr_jdx:	integer array of pointers that allows
c				fast retrieval of the values along each row (compressed)
c
c	For info on storage, see toolkit_hessian_cr.f
c	
c===================================================================================
c===================================================================================
c
	subroutine sprsym_cr_mult(spr_M,spr_row,spr_jdx,A,B,N)
c
c	Input:
c		
c		- spr_M	: 	real, non-zero values of matrix M
c		- spr_row:	row pointers for matrix M
c		- spr_jdx:	column pointers for sparse matrix
c		- A:		left vector
c		- N:		M is NxN, and A and B are of size N
c		- nmax:		maximum size of sparse vectors
c
c	Output:
c		- B:		result vector
c
	integer	i,j,k,N
	integer	spr_jdx(*),spr_row(N+1)
c
	real*8	spr_M(*),A(N),B(N)
c
c	Initialize result vector B: 
c
	do 100 i = 1,N
		B(i) = 0
100	continue
c
c	Now loop over all non zero elements of M:
c
	do 300 i = 1,N
c
		do 200 k = spr_row(i),spr_row(i+1)-1
			j = spr_jdx(k)
c
			B(i) = B(i) + spr_M(k)*A(j)
c
c			Recover symetric value
c
			if(j.ne.i) B(j) = B(j) + spr_M(k)*A(i)
c
200		continue
c
300	continue
c
c	Thats it !
c
	return
	end
c
c===================================================================================
c===================================================================================
c
c	Project_nm.f		Version 1 7/9/2010		Patrice Koehl
c
c	This routine projects a difference between two positions of a molecule
c	on the normal mode directions of that molecule
c
c	Input:
c		- ndim		: number of coordinates (3*nat)
c		- m		: twice the number of normal modes
c		- Eigenvect	: matrix of eigenvectors (size ndimxm)
c		- Xpos		: current position
c		- X0		: reference position (on which normal modes
c				  were computed)
c		- Diff		: work array
c
c	Output:
c		- amplitude	: projection of Xpos-X0 on the m/2 eigenvectors
c
c===================================================================================
c===================================================================================
c
	subroutine project_nm(ndim,m,Eigenvect,Xpos,X0,Diff,amplitude)
c
c===================================================================================
c	Local Variables - Arrays
c===================================================================================
c
	integer	m,nval,ndim,i,j
c
	real*8	Xpos(ndim),X0(ndim),Diff(ndim),amplitude(m)
	real*8	Eigenvect(ndim,m)
c
	do 100 i = 1,ndim
		Diff(i) = Xpos(i) - X0(i)
100	continue
c
	nval = m/2
	do 300 i = 1,nval
		amplitude(i) = 0.d0
		do 200 j = 1,ndim
			amplitude(i)=amplitude(i)+Eigenvect(j,i)*Diff(j)
200		continue
300	continue
c
	return
	end
c
