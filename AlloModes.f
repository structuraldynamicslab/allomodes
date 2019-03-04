c===================================================================================
c===================================================================================
c
c	AlloModes.f
c
c	This program computes paths for allostery using the elastic normal modes 
c	model
c
c	For references, see:
c
c	1. Monique Tirion, "Large amplitude elastic motions in proteins
c	from a single parameter, Atomic analysis". Phys. Rev. Letters, 
c	77, 1905-1908 (1996).
c
c	Copyright (C) 2017 Patrice Koehl and Marc Delarue
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
c	Include some constants
c===================================================================================
c
	include 'toolkit.h'
c
c===================================================================================
c	Declarations
c===================================================================================
c
	integer	i,j,idx
	integer	iat,jat
	integer	ierror,lwork,info,iprint
	integer	nchain,natom,nseq,nres,nkeep
	integer	nat,nang,ndihed
	integer	lun_rpdb,lun_rpdb2,lun,lun_bfac,lun_eln,lun_out
	integer	lun_log
	integer	flag_atom,flag_mass
	integer natom_select,noffset,besti
	integer	flag_network
	integer	nsize,ncont,nchange,ncycle,icycle
	integer	nmodes,nmode1,nmode2
	integer	nligand, nactive
	integer	natm
	integer	m
c
	integer	itype(*),listatom(*),listdihed(*),chain(*)
	integer	occup(*),atom_info(*)
	integer	corresp_select(*)
	integer	corresp2(*)
	integer	ncontact(*),listcontact(2,*)
	integer	nbound(*),iwa(*),isave(44)
c
	logical lsave(4)
c
	real*8  fact, sum
	real*8	xmass
	real*8	Kconst, Knew, cutoffE
	real*8	Temp, kb, hbar, t0, fac, fac2
	real*8	pi, kT, facb
	real*8	x, A, B, expA, expmA, ent, entropy
	real*8	rms,correl, besto
	real*8	score, dscore
	real*8	over5, over10, over, val
	real*8	Kmin,Kmax,Kmean
	real*8	factr, pgtol, tol
	real*8	deriv_kval(*)
	real*8	mass_atom(*),mass_sqrt(*)
	real*8	coord(*),coord_target(*),bfact(*),coord_select(*)
	real*8	dihed(*)
	real*8	diff_coord(*)
	real*8	coord_atm(*)
	real*8	hessian(*)
	real*8	overlaps(*),kval(*)
	real*8	work1(*),work2(*)
	real*8	work3(*),work4(*),work5(*)
	real*8	eigenval(*),eigenvect(*)
	real*8	bfact_calc(*)
	real*8	xval(*),lower(*),upper(*),wa(*),dsave(30)
c
	character	firstl*1
	character	fname_rpdb*100,fname_target*100,fname_res*100
	character	fname_ligand*100, fname_active*100
	character	fname*100
	character	seq(*)*4,atom(*)*4
	character	label(*)*30
	character	chname(10)*1
	character	csave(60)*1, task*60
c
c===================================================================================
c	Pointers for dynamically allocated arrays
c===================================================================================
c
	pointer	(ptr_itype,itype)
	pointer	(ptr_chain,chain)
	pointer	(ptr_seq,seq)
	pointer	(ptr_coord,coord)
	pointer	(ptr_coord_atm,coord_atm)
	pointer	(ptr_coord_target,coord_target)
	pointer	(ptr_bfact,bfact)
	pointer	(ptr_overlaps,overlaps)
	pointer	(ptr_bfact_calc,bfact_calc)
	pointer	(ptr_occup,occup)
	pointer	(ptr_listatom,listatom)
	pointer	(ptr_atom,atom)
	pointer	(ptr_label,label)
	pointer	(ptr_atom_info,atom_info)
	pointer	(ptr_dihed,dihed)
	pointer	(ptr_listdihed,listdihed)
	pointer	(ptr_mass_atom,mass_atom)
	pointer	(ptr_mass_sqrt,mass_sqrt)
	pointer	(ptr_corresp_select,corresp_select)
	pointer	(ptr_coord_select,coord_select)
	pointer	(ptr_diff_coord,diff_coord)
	pointer	(ptr_corresp2,corresp2)
c
	pointer	(ptr_hessian,hessian)
	pointer	(ptr_kval,kval)
	pointer	(ptr_deriv_kval,deriv_kval)
c
	pointer	(ptr_eigenval,eigenval)
	pointer	(ptr_eigenvect,eigenvect)
	pointer	(ptr_work1,work1)
	pointer	(ptr_work2,work2)
	pointer	(ptr_work3,work3)
	pointer	(ptr_work4,work4)
	pointer	(ptr_work5,work5)
c
	pointer	(ptr_ncontact,ncontact)
	pointer	(ptr_listcontact,listcontact)
c
	pointer	(ptr_xval,xval)
	pointer	(ptr_lower,lower)
	pointer (ptr_upper,upper)
	pointer	(ptr_nbound,nbound)
	pointer	(ptr_wa,wa)
	pointer	(ptr_iwa,iwa)
c
c===================================================================================
c	Common blocks
c===================================================================================
c
        common /protein/ 	ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /angle/   	ptr_dihed,ptr_listdihed,ndihed
        common /xyz/     	ptr_coord,ptr_occup,ptr_listatom,natom
        common /names/   	ptr_atom
	common /tags/		ptr_label
	common /infos/		ptr_atom_info
	common /contacts/ 	ptr_ncontact,ptr_listcontact,nsize,ncont
	common /elastic/	ptr_kval
c
c ==========================================================================================
c = Input/Output formats       =============================================================
c ==========================================================================================
c
2	format(80('='),/,'  Mode #  ',' Eigenval  ',
     1          '  Frequency        ','  entropy (this mode)'
     2          ' entropy (total)',/,'                     ',
     3          '   (cm^-1)           (kb unit)       ',
     4          '   (kb unit)',/,80('='))
3	format(4x,i4,4x,g12.5,4x,g12.5,4x,g12.5,4x,g12.5)
4	format(80('='))
5	format(60('='),/,'  Cycle #      Score        Kmin        Kmax',
     1  '       Kmean',/,60('='))
6	format(4x,i4,4x,f8.5,4x,f8.3,4x,f8.3,4x,f8.3)
7	format(60('='))
c
c ==================================================================================
c 	Some constants
c ==================================================================================
c
	pi = dacos(-1.d0)
	kT = 0.593d0
	facb = 8.d0*kT*pi*pi/3.d0
c
c===================================================================================
c 	Read all input information 
c===================================================================================
c
        call read_flags(fname_rpdb,Kconst,cutoffE,flag_network,
     1	nkeep,ncycle,flag_atom,flag_mass,fname_res,fname_target,
     2	fname_ligand,fname_active)
c
	idx = index(fname_res,'.')
	if(idx.eq.0) then 
		idx = index(fname_res,' ')
	endif
	idx = idx - 1
c
	nmode1 = 7
	nmode2 = nkeep
c
c===================================================================================
c 	Define log file for the run
c===================================================================================
c
	lun_log = 25
	fname=fname_res(1:idx)//'.log'
	open(unit=lun_log,file=fname,status='unknown')
c
c===================================================================================
c	1.1 Read amino acid topologies and parameters
c===================================================================================
c
	call topology_new
c
c===================================================================================
c      	1.2 count number of residues in protein considered, and 
c	    allocate arrays accordingly
c===================================================================================
c
	lun_rpdb = 11
	open(unit=lun_rpdb,file=fname_rpdb,status='old')
c
	call defchain(lun_rpdb,nchain,chname,nres)
	rewind(unit=lun_rpdb)
c
	nseq = nres + 2*nchain
c
	ptr_itype	= malloc(4*nseq)
	ptr_seq		= malloc(4*nseq)
	ptr_chain	= malloc(4*nseq)
	ptr_listatom	= malloc(4*nseq)
	ptr_listdihed	= malloc(4*nseq)
c
c===================================================================================
c      	1.3 Read in sequence, define numbers of atoms, and allocate arrays accordingly
c===================================================================================
c
	call extractseq(lun_rpdb)
	rewind(unit=lun_rpdb)
c
	call count_atom_dihed(nat,nang)
c
	natom = nat
	ndihed = nang
c
        write(6,*) ' '
        write(6,*) 'nseq      = ',nseq-2*nchain
        write(6,*) 'nat (all) = ',nat
        write(6,*) ' '
c
	ptr_coord	= malloc(8*3*natom)
	ptr_bfact	= malloc(8*natom)
	ptr_occup	= malloc(4*natom)
	ptr_atom	= malloc(4*natom)
	ptr_dihed	= malloc(8*ndihed)
	ptr_mass_atom	= malloc(8*natom)
	ptr_mass_sqrt	= malloc(8*natom)
	ptr_label	= malloc(30*natom)
	ptr_atom_info	= malloc(8*natom)
c
	call prepare_protein
c
c===================================================================================
c       1.4 Read protein (initial and target conformations)
c===================================================================================
c
	call set_zeros(natom,bfact)
c
	lun_rpdb2 = 13
	ptr_coord_target = malloc(8*3*natom)
	open(unit=lun_rpdb2,file=fname_target,status='old')
	call readpdb_new(lun_rpdb2,ierror,bfact)
	call dcopy(3*natom,coord,1,coord_target,1)
	close(unit=lun_rpdb2)
c
	call readpdb_new(lun_rpdb,ierror,bfact)
	close(unit=lun_rpdb)
c
c===================================================================================
c	1.5 Define mass of all atoms
c===================================================================================
c
	do 100 i = 1,natom
		firstl=atom(i)(1:1)
		mass_atom(i) = xmass(firstl)
		if(flag_mass.eq.0) mass_atom(i) = 1.d0
		mass_sqrt(i) = sqrt(mass_atom(i))
100	continue
c
c===================================================================================
c	2. Select atoms for potential
c===================================================================================
c
	ptr_corresp_select = malloc(4*natom)
	if(flag_atom.eq.0) then
		call get_ca(natom_select, corresp_select)
	elseif(flag_atom.eq.1) then
		call get_backbone(natom_select, corresp_select)
	elseif(flag_atom.eq.2) then
		call get_heavy(natom_select, corresp_select)
	endif
c
c===================================================================================
c       2.1 Read ligand pocket and active site (channel)
c	    For efficiency, the array of selected atoms is organized as follows:
c		1-nactive:   active site(s)
c		nactive+1-nactive+nligand:  ligand
c		>nactive+nligand+1:         remaining atoms
c===================================================================================
c
	ptr_coord_atm = malloc(3*8*natom)
	ptr_corresp2  = malloc(4*natom)
c
	lun = 13
	open(unit=lun,file=fname_active,status='unknown')
	call read_atom(lun,natm,coord_atm)
	close(unit=lun)
c
	noffset = 0
	call identify_region(coord,natom_select,corresp_select,
     1		coord_atm,natm,noffset,corresp2,nactive)
c
	open(unit=lun,file=fname_ligand,status='unknown')
	call read_atom(lun,natm,coord_atm)
	close(unit=lun)
c
	noffset = nactive
	call identify_region(coord,natom_select,corresp_select,
     1		coord_atm,natm,noffset,corresp2,nligand)
c
	call free(ptr_coord_atm)
c
	write(6,*) ' '
	write(6,*) 'Number of atoms selected        : ',natom_select
	write(6,*) 'Number of atoms in ligand pocket: ',nligand
	write(6,*) 'Number of atoms in active site  : ',nactive
	write(6,*) ' '
c
        write(lun_log,*) ' '
	write(lun_log,*) 'Input PDB file: start structure  : ',
     1  fname_rpdb(1:index(fname_rpdb,' ')-1)
	write(lun_log,*) 'Input PDB file: target structure : ',
     1  fname_target(1:index(fname_target,' ')-1)
        write(lun_log,*) ' '
        write(lun_log,*) '# of residues in protein         : ',
     1		nseq-2*nchain
        write(lun_log,*) 'Tot. # of atoms (includ. H)      : ',natom
c
	write(lun_log,*) '# of atoms selected              : ',
     1		natom_select
	write(lun_log,*) '# of atoms in ligand pocket      : ',nligand
	write(lun_log,*) '# of atoms in active site        : ',nactive
	if(flag_mass.eq.1) then
		write(lun_log,*) 'Normal mode type                 :',
     1		' mass-weighted'
	else
		write(lun_log,*) 'Normal mode type                 :',
     1		' unweighted'
	endif
	write(lun_log,*) ' '
c
c===================================================================================
c       2.2 corresp_select is from selected -> all; set corresp2 from all -> select2
c===================================================================================
c
	call def_corresp(nat,corresp2,natom_select,corresp_select)
c
c===================================================================================
c       2.3 Get selected atoms, and for those, displacements between start and target
c	    conformations
c	    The displacements for the active sites are normalized
c===================================================================================
c
	ptr_coord_select = malloc(8*3*natom_select)
	ptr_diff_coord = malloc(8*3*natom_select)
	call get_select(natom_select,corresp_select,coord_select)
	call get_diff(natom_select,corresp_select,coord_select,
     1  coord_target,diff_coord)
c
	call normalize(nactive,diff_coord)
c
c===================================================================================
c	3. Define elastic network
c===================================================================================
c
	if(flag_network.eq.1) then
		if(cutoffE.gt.0) then
			call contact(coord_select,natom_select,cutoffE)
		else
			call contact_all(coord_select,natom_select)
		endif
	else
		call listcontact_del(natom_select,coord_select)
	endif
c
	ptr_kval = malloc(8*ncont)
	ptr_deriv_kval = malloc(8*ncont)
c
	do 200 i = 1,ncont
		kval(i) = Kconst
200	continue
c
c===================================================================================
c	5. Compute Hessian of elastic energy
c===================================================================================
c
	ptr_hessian = malloc(8*9*natom_select*natom_select)
c
	call hessian_elastic_full(natom_select,corresp_select,mass_sqrt,
     1		coord_select,hessian)
c
c===================================================================================
c	8. Solve eigenvalue problem:
c		H x = lambda x
c	where H is the Hessian matrix
c
c	compute all modes
c===================================================================================
c
	ptr_eigenval = malloc(8*3*natom_select)
	ptr_eigenvect = malloc(8*3*natom_select*3*natom_select)
c
	call dcopy(9*natom_select*natom_select,hessian,1,eigenvect,1)
c
	lwork = -1
	call dsyev('V','U',3*natom_select,eigenvect,3*natom_select,
     2  eigenval,fact,lwork,info)
	lwork = int(fact)
	ptr_work1 = malloc(8*lwork)
	call dsyev('V','U',3*natom_select,eigenvect,3*natom_select,
     2	eigenval,work1,lwork,info)
	nmodes = 3*natom_select
	if(nkeep.eq.0) nmode2 = nmodes

	call rescale_eigvect(natom_select,corresp_select,mass_sqrt,
     1		eigenvect,nmodes)
c
c ==========================================================================================
c =      9. Compute atomic fluctuations and correlations (for CA only) =====================
c ==========================================================================================
c
	sum = 0.d0
	do 300 i = 1,nat
		sum = sum + bfact(i)
300	continue
c
	if(sum.ne.0.d0) then
c
		ptr_bfact_calc	= malloc(8*natom)
		call compute_bfact(natom_select,nmode1,nmode2,
     1		eigenval,eigenvect,bfact_calc,facb)
		call compare_bfact(natom_select,corresp_select,bfact,
     1		bfact_calc,rms,correl,Knew)
c
		write(6,*) ' '
		write(6,*) 'K (constant) optimized         : ',Knew
		write(6,*) 'RMS between exp and calc Bfact : ',rms
		write(6,*) 'CC between exp and calc Bfact  : ',correl
		write(6,*) ' '
		write(lun_log,*) ' '
		write(lun_log,*) 'K (constant) optimized         : ',
     1			Knew
		write(lun_log,*) 'RMS between exp and calc Bfact : ',
     1			rms
		write(lun_log,*) 'CC between exp and calc Bfact  : ',
     1			correl
		write(lun_log,*) ' '
c
		lun_bfac = 14
		fname=fname_res(1:idx)//'_bfact.dat'
		open(unit=lun_bfac,file=fname,status='unknown')
		call write_bfact(lun_bfac,natom,bfact,bfact_calc,
     1		label,corresp2)
		close(unit=lun_bfac)
c
		call free(ptr_bfact_calc)
c
		do 350 i = 1,ncont
			kval(i) = Knew
350		continue
c
		call hessian_elastic_full(natom_select,corresp_select,
     1		mass_sqrt,coord_select,hessian)
c
		call dcopy(9*natom_select*natom_select,hessian,1,
     1		eigenvect,1)
c
		call dsyev('V','U',3*natom_select,eigenvect,
     1		3*natom_select,eigenval,work1,lwork,info)
		call rescale_eigvect(natom_select,corresp_select,
     1		mass_sqrt,eigenvect,nmodes)
	endif
c
c ==========================================================================================
c	13. Compute overlap between modes and atomic displacements between two structures
c ==========================================================================================
c
	ptr_overlaps = malloc(8*nmodes)
c
	call compute_overlap(nmode1,nmode2,nmodes,natom_select,
     1	corresp_select,eigenvect,coord_target,mass_atom,overlaps,
     2	besto,besti)
c
	over5 = 0.d0
	over10 = 0.d0
	over = 0.d0
	do 375 i = nmode1,nmode2
		val = overlaps(i)*overlaps(i)
		if(i-nmode1.lt.5) then
			over5 = over5 + val
		endif
		if(i-nmode1.lt.10) then
			over10 = over10 + val
		endif
		over = over + val
375	continue
	over5  = sqrt(over5)
	over10 = sqrt(over10)
	over   = sqrt(over)
c
	write(6,*) 'Before refinement:'
	write(6,*) 'Best overlap start/target      : ',besto
	write(6,*) 'Mode # for best overlap        : ',besti
	write(6,*) 'Cumul. overlap 5 modes         : ',over5
	write(6,*) 'Cumul. overlap 10 modes        : ',over10
	write(6,*) 'Cumul. overlap up to nkeep     : ',over
	write(6,*) ' '
	write(lun_log,*) 'Before refinement:'
	write(lun_log,*) 'Best overlap start/target      : ',besto
	write(lun_log,*) 'Mode # for best overlap        : ',besti
	write(lun_log,*) 'Cumul. overlap 5 modes         : ',over5
	write(lun_log,*) 'Cumul. overlap 10 modes        : ',over10
	write(lun_log,*) 'Cumul. overlap up to nkeep     : ',over
	write(lun_log,*) ' '
c
c ==========================================================================================
c =      11. Compute correlation between atomic fluctuations in ligand pocket and in active   
c =         site
c ==========================================================================================
c
	ptr_work2 = malloc(8*3*nactive)
	ptr_work3 = malloc(8*3*nactive)
	ptr_work4 = malloc(8*3*(nactive+nligand))
	ptr_work5 = malloc(8*3*nmodes)
c
	call score_allostery(natom_select,nactive,nligand,
     1	diff_coord,nmode1,nmode2,eigenval,eigenvect,work2,score)
c
	write(6,*) ' '
	write(6,*) 'Initial allosteric score  : ',score
	write(6,*) 'Number of kconst to refine: ',ncont
	write(6,*) 'First/last mode considered : ',7,nkeep
	write(6,*) ' '
	write(lun_log,*) ' '
	write(lun_log,*) 'Initial allosteric score   : ',score
	write(lun_log,*) 'Number of kconst to refine : ',ncont
	write(lun_log,*) 'First/last mode considered : ',7,nkeep
	write(lun_log,*) ' '
c
c ==========================================================================================
c =     10. Compute score and its derivatives for one kval            ======================
c ==========================================================================================
c
	iprint = -1
	factr = 1.0d7
	pgtol = 1.d-5
	tol = 0.01d0
	m = 5
c
	ptr_xval = malloc(8*ncont)
	ptr_lower = malloc(8*ncont)
	ptr_upper = malloc(8*ncont)
	ptr_nbound = malloc(4*ncont)
	ptr_wa = malloc(8*((2*m+5)*ncont + 11*m*m + 8*m))
	ptr_iwa = malloc(4*3*ncont)
c
	do 400 i = 1,ncont
		xval(i) = kval(i)
		lower(i) = 0.d0
		nbound(i) = 1
400	continue
c
	task = 'START'
c
	icycle = 0
c
	write(6,*) ' '
	write(6,5)
	write(lun_log,*) ' '
	write(lun_log,5)
c
	do 500 nchange = 1,ncont
c
		iat = listcontact(1,nchange)
		jat = listcontact(2,nchange)
c
		call dscore_allostery(iat,jat,natom_select,
     1		coord_select,nactive,nligand,diff_coord,
     2		nmodes,nmode1,nmode2,eigenval,eigenvect,
     3		work2,work3,work1,work4,work5,score,dscore)
c
		deriv_kval(nchange) = dscore
c
500	continue
c
600	continue
c
		call setulb(ncont, m, xval, lower, upper, nbound, 
     1		score, deriv_kval, factr, pgtol, wa, iwa, task, iprint, 
     2		csave, lsave, isave,dsave)
c
		if(task(1:2).eq.'FG') then
c
			Kmean = 0.d0
			do 700 nchange = 1,ncont
				kval(nchange) = xval(nchange)
				if(nchange.eq.1) then
					Kmin = kval(nchange)
					Kmax = kval(nchange)
				else
					if(kval(nchange).lt.Kmin) then
						Kmin = kval(nchange)
					endif
					if(kval(nchange).gt.Kmax) then
						Kmax = kval(nchange)
					endif
				endif
				Kmean = Kmean + kval(nchange)
700			continue
			Kmean = Kmean/ncont
c
			call hessian_elastic_full(natom_select,
     1			corresp_select,mass_sqrt,coord_select,hessian)
c
			call dcopy(9*natom_select*natom_select,hessian,
     1			1,eigenvect,1)
			call dsyev('V','U',3*natom_select,eigenvect,
     1			3*natom_select,eigenval,work1,lwork,info)
c
			call rescale_eigvect(natom_select,
     1			corresp_select,mass_sqrt,eigenvect,nmodes)
c
			call score_allostery(natom_select,nactive,
     1			nligand,diff_coord,nmode1,nmode2,eigenval,
     2			eigenvect,work2,score)
c
			do 800 nchange = 1,ncont
c
				iat = listcontact(1,nchange)
				jat = listcontact(2,nchange)
c
				call dscore_allostery(iat,jat,
     1				natom_select,coord_select,nactive,
     2				nligand,diff_coord,nmodes,nmode1,
     3				nmode2,eigenval,eigenvect,work2,work3,
     4				work1,work4,work5,score,dscore)
c
				deriv_kval(nchange) = dscore
c
800			continue
c
			icycle = icycle + 1
c
			write(6,6) icycle,score,Kmin,Kmax,Kmean
			write(lun_log,6) icycle,score,Kmin,Kmax,Kmean
c
			if(icycle.eq.ncycle) goto 900
			if(score.le.tol) goto 900
c
			goto 600
c
		elseif(task(1:5).eq.'NEW_X') then
c
			goto 600
c
		endif
c
900	continue
c
	write(6,7)
	write(6,*) ' '
	write(lun_log,7)
	write(lun_log,*) ' '
c
	call free(ptr_xval)
	call free(ptr_lower)
	call free(ptr_upper)
	call free(ptr_nbound)
	call free(ptr_wa)
	call free(ptr_iwa)
c
c			
c ==========================================================================================
c =     11. Compute conformational entropy for final elastic network =======================
c
c	Calculation is based on:
c
c	I. Anddricioaei and M. Karplus, "On the calculation of emtropy from covariance 
c	matrices of the atomic fluctuations", J. Chem. Phys., 115, 6289-6292 (2001).
c
c	Constants:
c
c	Temp is the temperature in Kelvin
c	hbar is h/(2pi), where h is the Plank constant; h = 1.0546 10^(-34) J.s
c	kb is the Boltzmann constant: kb = 1.380662 10^{-23} J/K
c	t0 is the time unit in the AKMA system: t0 = 4.88882 10^(-14) s
c ==========================================================================================
c
	Temp = 300
	hbar = 1.0546d-34
	kb   = 1.3807d-23
	t0   = 4.88882d-14
c
	fac = hbar/(Temp*kb*t0)
	fac2 = 108.59d0
c
	entropy = 0
	write(6,2)
	write(lun_log,2)
	do 1000 i = 1,nmode2
		x = dsqrt(abs(eigenval(i)))
		A = x*fac
		B = x*fac2
		expA = exp(A)
		expmA = 1.d0/expA
		if(i.gt.6) then
			ent = A/(expA-1) - log(1-expmA)
		else
			ent = 0.d0
		endif
		entropy = entropy + ent
		write(6,3) i,abs(eigenval(i)),B,ent,entropy
		write(lun_log,3) i,abs(eigenval(i)),B,ent,entropy
1000	continue
	write(6,4)
	write(lun_log,4)
	write(lun_log,*) ' '
c
c ==========================================================================================
c	Compute overlap between modes and atomic displacements between two structures
c ==========================================================================================
c
	fname=fname_res(1:idx)//'_overlap.dat'
c
	lun_out = 11
	open(unit=lun_out,file=fname,status='unknown')
c
	call compute_overlap(nmode1,nmode2,nmodes,natom_select,
     1	corresp_select,eigenvect,coord_target,mass_atom,overlaps,
     2	besto,besti)
c
	call write_overlap(lun_out,nmode1,nmode2,overlaps)
c
	close(unit=lun_out)
c
	over5 = 0.d0
	over10 = 0.d0
	over = 0.d0
	do 1050 i = nmode1,nmode2
		val = overlaps(i)*overlaps(i)
		if(i-nmode1.lt.5) then
			over5 = over5 + val
		endif
		if(i-nmode1.lt.10) then
			over10 = over10 + val
		endif
		over = over + val
1050	continue
	over5  = sqrt(over5)
	over10 = sqrt(over10)
	over   = sqrt(over)
c
	write(6,*) 'After refinement:'
	write(6,*) 'Best overlap start/target      : ',besto
	write(6,*) 'Mode # for best overlap        : ',besti
	write(6,*) 'Cumul. overlap 5 modes         : ',over5
	write(6,*) 'Cumul. overlap 10 modes        : ',over10
	write(6,*) 'Cumul. overlap up to nkeep     : ',over
	write(6,*) ' '
	write(lun_log,*) 'After refinement:'
	write(lun_log,*) 'Best overlap start/target      : ',besto
	write(lun_log,*) 'Mode # for best overlap        : ',besti
	write(lun_log,*) 'Cumul. overlap 5 modes         : ',over5
	write(lun_log,*) 'Cumul. overlap 10 modes        : ',over10
	write(lun_log,*) 'Cumul. overlap up to nkeep     : ',over
	write(lun_log,*) ' '
c
	call free(ptr_overlaps)
c
c===================================================================================
c	Write final elastic network
c===================================================================================
c
	lun_eln = 14
	fname=fname_res(1:idx)//'.eln'
	open(unit=lun_eln,file=fname,status='unknown')
	call write_eln(lun_eln,natom_select,corresp_select,mass_sqrt,
     1		coord,bfact,fname_rpdb,flag_network,cutoffE)
	close(lun_eln)
c
	fname=fname_res(1:idx)//'.pml'
	open(unit=lun_eln,file=fname,status='unknown')
	call write_pml(lun_eln,natom_select,coord,corresp_select)
	close(lun_eln)
c
	fname=fname_res(1:idx)//'5.pml'
	open(unit=lun_eln,file=fname,status='unknown')
	call write_pml10(lun_eln,natom_select,coord,corresp_select)
	close(lun_eln)
c
	fname=fname_res(1:idx)//'.pdb'
	open(unit=lun_eln,file=fname,status='unknown')
	call write_pdb(lun_eln,natom,natom_select,coord,label,
     1		corresp2,bfact,occup)
	close(lun_eln)
c
	close(lun_log)
c
c===================================================================================
c	Cleanup arrays
c===================================================================================
c
	call free(ptr_corresp_select)
	call free(ptr_corresp2)
	call free(ptr_mass_atom)
	call free(ptr_eigenvect)
	call free(ptr_eigenval)
	call free(ptr_bfact)
	call free(ptr_coord_target)
	call free(ptr_coord_select)
	call free(ptr_diff_coord)
	call free(ptr_hessian)
	call free(ptr_mass_sqrt)
	call free(ptr_work1)
	call free(ptr_work2)
	call free(ptr_work3)
	call free(ptr_work4)
	call free(ptr_work5)
c
	stop
	end
c ===========================================================================================
c ===========================================================================================
c	read_flags.f
c ===========================================================================================
c ===========================================================================================
c
c	This subroutine "reads" the program inputs (given on the command line)
c
c ===========================================================================================
c ===========================================================================================
c
	subroutine read_flags(file_pdb,Kconst,cutoff,flag_network,
     1  nkeep,ncycle,flag_atom,flag_mass,fname_res,fname_target,
     2	fname_ligand,fname_active)
c
c ===========================================================================================
c ===========================================================================================
c
	integer	i,narg
	integer	nkeep,ncycle,flag_network,flag_atom,flag_mass
c
	real*8	Kconst,cutoff
c
	character	file_pdb*100,fname_res*100,fname_target*100
	character	fname_ligand*100,fname_active*100
	character	OPTION*1
	character	List(100)*100
c
	logical		Input,Skip
c
1	format(2x,90('='))
2	format(2x,'=',88x,'=')
3	format(2x,'=',40x,'AlloModes',39x,'=')
4	format(2x,'=',5x,'AlloModes finds allosteric paths using ',
     1  ' the concept of normal mode analysis    ',4x,'=')
5	format(2x,'=',5x,'Usage:',77x,'=',/,2x,'=',88x,'=',/,
     1  2x,'=',5x,'AlloModes.exe -i <pdb file> -o <res. file>',
     2  ' -t <target> -l <ligand> -s <active>',5x,'=',/,2x,'=',5x,
     3  ' [ -e network> -k <Kconst> -c <cutoff> -a <flag_atom>',
     4  ' -m <flag_mass> -n <nmodes> ',2x,'='/,2x,'=',5x,
     5  ' -r <ncycles> ]',68x,'=')
6	format(2x,'=',5x,'with:',78x,'=',/,
     1  2x,'=',10x,'-i <pdb file> ......path name to input PDB file',
     2  31x,'=',/,
     3  2x,'=',10x,'-o <res. file> .....path name to result file',
     4  34x,'=',/,
     &  2x,'=',10x,'-t <target file>... path name to PDB file for',
     &  ' target struct.               ',3x,'=',/,
     &  2x,'=',10x,'-l <ligand> ........PDB file for ligand pocket',
     &  32x,'=',/,
     &  2x,'=',10x,'-s <active> ........PDB file for active site',
     &  34x,'=',/,
     7  2x,'=',5x,'Optional arguments:',64x,'=',/,
     &  2x,'=',10x,'-e <network>........cutoff-based (1) or Delaunay',
     &  ' (2) network (default 1)',6x,'=',/,
     &  2x,'=',10x,'-k <Kconst> ........Force constant for springs ',
     &  '(elastic pot.) (default 1.0)',3x,'=',/,
     &  2x,'=',10x,'-c <cutoff> ........cutoff for elastic network',
     &  ' (only if e = 1) (default 10.0)',1x,'=',/,
     &  2x,'=',10x,'-a <flag_atom>..... select CA (0), backbone (1)',
     &  ' or all heavy atoms (2)',8x,'=',/,
     &  2x,'=',29x,' for elastic network (default 0)',27x,'=',/,
     &  2x,'=',10x,'-m <flag_mass>..... mass-weighted (1) or ',
     &  ' not (0) (default 1)',17x,'=',/,
     &  2x,'=',10x,'-n <nmodes>........ number of modes to ',
     &  ' keep (default 100)',20x,'=',/,
     &  2x,'=',10x,'-r <ncycles>....... number of cycles to ',
     &  ' refine kvals (default 100)',11x,'=')
7	format('Allowed options: i,o,n,e,k,c,a,m.',
     1  ' You used: ',a1)
c
	narg = iargc()
c
	if(narg.eq.0) then
		write(6,*) ' '
		write(6,1)
		write(6,1)
		write(6,2)
		write(6,3)
		write(6,2)
		write(6,4)
		write(6,2)
		write(6,5)
		write(6,2)
		write(6,6)
		write(6,2)
		write(6,1)
		write(6,1)
		write(6,*) ' '
		stop
	endif
c
	do 100 i = 1,narg
		call getarg(i,List(i))
100	continue
c
	Input = .false.
c
	Kconst       = 1.d0
	cutoff       = 10.d0
	flag_network = 1
	flag_atom    = 0
	fname_target = ' '
	flag_mass    = 1
	nkeep        = 26
	ncycle       = 100
c
	Skip  = .false.
c
	do 200 i = 1,narg
c
		if(Skip) then
			Skip = .false.
			goto 200
		endif
c
		if(List(i)(1:1).eq.'-') then
c
			OPTION=List(i)(2:2)
			if(OPTION.eq.'i') then
				file_pdb = List(i+1)
				Skip = .true.
				Input = .true.
			elseif(OPTION.eq.'o') then
				fname_res = List(i+1)
				Skip = .true.
			elseif(OPTION.eq.'k') then
				read(List(i+1),*) Kconst
				Skip = .true.
			elseif(OPTION.eq.'c') then
				read(List(i+1),*) cutoff
				Skip = .true.
			elseif(OPTION.eq.'a') then
				read(List(i+1),*) flag_atom
				Skip = .true.
			elseif(OPTION.eq.'e') then
				read(List(i+1),*) flag_network
				Skip = .true.
			elseif(OPTION.eq.'m') then
				read(List(i+1),*) flag_mass
				Skip = .true.
			elseif(OPTION.eq.'t') then
				fname_target = List(i+1)
				Skip = .true.
			elseif(OPTION.eq.'l') then
				fname_ligand = List(i+1)
				Skip = .true.
			elseif(OPTION.eq.'s') then
				fname_active = List(i+1)
				Skip = .true.
			elseif(OPTION.eq.'n') then
				read(List(i+1),*) nkeep
				Skip = .true.
			elseif(OPTION.eq.'r') then
				read(List(i+1),*) ncycle
				Skip = .true.
			else
				write(6,*) ' '
				write(6,7) OPTION
				write(6,*) ' '
				write(6,1)
				write(6,1)
				write(6,2)
				write(6,3)
				write(6,2)
				write(6,4)
				write(6,2)
				write(6,5)
				write(6,6)
				write(6,2)
				write(6,1)
				write(6,1)
				stop
			endif
c
		endif
c
200	continue
c
	if(.not.Input) then
		write(6,*) ' '
		write(6,*) 'You did not provide an input file'
		write(6,*) ' '
		stop
	endif
c
	return
	end
