c ==========================================================================================
c ==========================================================================================
c	Readpdb.f		Version 1 5/5/1992	Patrice Koehl
c
c	This program reads the present state of the molecule in
c	PDB format
c	At this stage, only one monomer of the molecule is considered
c
c ==========================================================================================
c ==========================================================================================
c	Readca.f		Version 1 5/5/1992	Patrice Koehl
c
c	This program reads the present state of the molecule in
c	PDB format (only CA atoms!)
c ==========================================================================================
c ==========================================================================================
c
	subroutine readca(lun_pdb,ierror)
c
	include 'toolkit.h'
c
	real*8	x,y,z
	real*8	coord(*)
c
	integer	nseq,nat,i,natom,ierr,nchain
	integer	nres,lun_pdb,ierror
	integer	itype(*),listatom(*)
	integer	chain(*),check(*)
c
	character name*4,name1*4,record*80
	character ch,ch1
	character ntest*5,idum2*5
	character atom(*)*4,seq(*)*4
	character label(*)*30
c
	pointer	(ptr_itype,itype)
	pointer	(ptr_listatom,listatom)
	pointer	(ptr_chain,chain)
	pointer	(ptr_check,check)
	pointer	(ptr_coord,coord)
	pointer	(ptr_seq,seq)
	pointer	(ptr_atom,atom)
	pointer	(ptr_atom_label,label)
c
c ==========================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /xyz/	 ptr_coord,ptr_check,ptr_listatom,natom
	common /names/	 ptr_atom
	common /tags/	 ptr_atom_label
c ==========================================================================================
c
1	format(12x,a4,1x,4x,a1,a5,3x,3f8.3)
2	format(a)
c
	ierror = 0
c
	do 100 i = 1,natom
		check(i) = 0
100	continue
c
	ch1 = 'x'
	nres = 1
	ntest = '     '
c
200     read(lun_pdb,2,end=300) record
        if(record(1:6).ne.'ATOM  ') goto 200
        read(record,1) name1,ch,idum2,x,y,z
	if(idum2.ne.ntest) then
c		if(ntest.ne.-10.and.idum2-ntest.ne.1) then
c			write(6,*) 'Warning : non consecutive residues !!'
c			write(6,*) 'Residues : ',ntest,idum2
c		endif
		ntest = idum2
		if(nres.eq.1) ch1 = ch
		nres = nres + 1
		if(ch.ne.ch1) then
			ch1 = ch
			chain(nres) = 1
			chain(nres-1) = 2
		else
			chain(nres) = 0
		endif
	endif
        if(name1(1:1).eq.' ') then
                name(1:4) = name1(2:4)//' '
        else
                name = name1
        endif
	if(name.ne.'CA  ') goto 200
	call atompos(nres,name,nat,ierr)
c	if(ierr.eq.0.or.ierr.eq.-1) goto 100
	if(check(nat).ne.0) goto 200
	check(nat) = check(nat) + 1
	coord(3*nat-2) = x
	coord(3*nat-1) = y
	coord(3*nat)   = z
	label(nat) = record(1:30)
	goto 200
300	continue
c
	if(nres.ne.nseq-1) then
		write(6,*) 'nres,nseq :',nres,nseq
		write(6,*) 'inconsistence in residue number !!'
		stop
	endif
c
	chain(nseq) = 2
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Readpdb_new.f		Version 1 5/5/1992	Patrice Koehl
c
c	This program reads the present state of the molecule in
c	PDB format
c	At this stage, only one monomer of the molecule is considered
c ==========================================================================================
c ==========================================================================================
c
	subroutine readpdb_new(lun_pdb,ierror,bfact)
c
	include 'toolkit.h'
c
	real*8	x,y,z,b
	real*8	coord(*),bfact(*)
c
	integer	nseq,nat,i,natom,ierr,nchain,nch
	integer	nres,lun_pdb,ierror,nblank
	integer	itype(*),listatom(*)
	integer	check(*)
	integer	chain(*)
	integer	atom_info(2,*)
c
	character name*4,name1*5,record*80
	character ch,ch1
	character ntest*5,idum2*5,resname*5
	character atom(*)*4,seq(*)*4
	character label(*)*30
c
	pointer	(ptr_itype,itype)
	pointer	(ptr_listatom,listatom)
	pointer	(ptr_chain,chain)
	pointer	(ptr_check,check)
	pointer	(ptr_coord,coord)
	pointer	(ptr_seq,seq)
	pointer	(ptr_atom,atom)
	pointer	(ptr_atom_label,label)
	pointer	(ptr_atom_info,atom_info)
c
c ==========================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /xyz/     ptr_coord,ptr_check,ptr_listatom,natom
	common /names/   ptr_atom
	common /tags/	 ptr_atom_label
	common /infos/	 ptr_atom_info
c ==========================================================================================
c
1	format(12x,a5,a5,a5,3x,3f8.3,6x,f6.2)
2	format(a)
c
	ierror = 0
c
	do 100 i = 1,natom
		check(i) = 0
100	continue
c
	ch1 = 'x'
	nres = 1
	chain(nres) = 1
	ntest = '     '
	nch = 1
c
200     read(lun_pdb,2,end=300) record
	if(record(1:6).eq.'ENDMDL') goto 300
        if(record(1:6).ne.'ATOM  ') goto 200
        read(record,1) name1,resname,idum2,x,y,z,b
	nblank = 0
	do 150 i = 1,4
		if(resname(i:i).eq.' ') nblank = nblank + 1
150	continue
	if(nblank.ge.2) goto 200
	read(resname(5:5),'(a1)') ch
	if(idum2.ne.ntest) then
c		if(ntest.ne.-10.and.idum2-ntest.ne.1) then
c			write(6,*) 'Warning : non consecutive residues !!'
c			write(6,*) 'Residues : ',ntest,idum2
c		endif
		ntest = idum2
		if(nres.eq.1) ch1 = ch
		if(ch.ne.ch1) then
			ch1 = ch
			nch = nch + 1
			nres = nres + 1
			chain(nres) = 2
			nres = nres + 1
			chain(nres) = 1
			nres = nres + 1
			chain(nres) = 0
		else
			nres = nres + 1
			chain(nres) = 0
		endif
	endif
	if(name1.eq.' H   ') name1=' HN  '
        if(name1(1:1).eq.' ') then
                name(1:4) = name1(2:5)
        else
		if(name1(1:1).ne.'H') then
			name=name1(2:4)//name1(1:1)
		else
                	name = name1(1:4)
		endif
        endif
	call atompos(nres,name,nat,ierr)
c	if(ierr.eq.0.or.ierr.eq.-1) goto 100
c	if(ierr.eq.0) then
c		write(6,*) name,'  ',name1,' could not be identified'
c	endif
	if(check(nat).ne.0) goto 200
	check(nat) = check(nat) + 1
	coord(3*nat-2) = x
	coord(3*nat-1) = y
	coord(3*nat)   = z
	label(nat) = record(1:30)
	atom_info(1,nat) = nch
	atom_info(2,nat) = nres
	bfact(nat) = b
	goto 200
300	continue
c
	if(nres.ne.nseq-1) then
		write(6,*) 'nres,nseq-1 :',nres,nseq-1
		write(6,*) 'inconsistence in residue number !!'
		stop
	endif
c
	chain(nseq) = 2
c
	nat = 0
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	extractseq.f		Version 1 24/7/1995	Patrice Koehl
c
c	this subroutine scans the pdb file to get the sequence of the
c	protein
c	Two 'pseudo' residues are added :
c	Nter and Cter (for HT1 and HT2, and for OT2, respectively)
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine extractseq(lun_pdb)
c
	include 'toolkit.h'
c
	integer	nseq,nchain,i,residue,lun_pdb,nblank
	integer	itype(*),chain(*)
c
	character resname*4,record*80
	character ch*1,ch1*1
	character itest*5,iread*5
	character seq(*)*4
c
	pointer	(ptr_itype,itype)
	pointer	(ptr_chain,chain)
	pointer	(ptr_seq,seq)
c
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
c
1	format(17x,a4,a1,a5)
2	format(a)
c
	ch1 = 'x'
	nseq  = 1
	seq(1) = 'NTER'
	itest = '     '
100     read(lun_pdb,2,end=200) record
	if(record(1:6).eq.'ENDMDL') goto 200
        if(record(1:6).ne.'ATOM  ') goto 100
        read(record,1) resname,ch,iread
c
c ==========================================================================================
c	Amino acids are usually stored as 3 letter code...if name has less than three letters,
c	discard
c ==========================================================================================
c
	nblank = 0
	do 150 i = 1,4
		if(resname(i:i).eq.' ') nblank = nblank + 1
150	continue
	if(nblank.ge.2) goto 100
c
	if(iread.ne.itest) then
		itest = iread
		if(nseq.eq.1) ch1 = ch
		if(ch.ne.ch1) then
			nseq = nseq + 1
			seq(nseq) = 'CTER'
			nseq = nseq + 1
			seq(nseq) = 'NTER'
			ch1 = ch
		endif
		nseq = nseq + 1
		seq(nseq) = resname
	endif
	goto 100
200	continue
	nseq = nseq + 1
	seq(nseq) = 'CTER'
c
        do 300 i = 1,nseq
                itype(i) = residue(seq(i))
300     continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Defchain.f		Version 1 5/5/1992	Patrice Koehl
c
c	This program counts the number of chain in the protein
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine defchain(lun_pdb,nchain,chname,nres)
c
	integer	nchain,lun_pdb
	integer	idum2,nres,ntest
c
	character record*80,chname(10)*1,ch1
c
1	format(a)
2       format(22x,i4)
c
	ch1 	= 'x'
	nchain 	= 0
	nres	= 0
	ntest 	= -10
c
100     read(lun_pdb,1,end=200) record
	if(record(1:6).eq.'ENDMDL') goto 200
        if(record(1:6).ne.'ATOM  ') goto 100
	if(record(27:27).ne.' ') goto 100
	if(record(22:22).ne.ch1) then
		nchain 			= nchain + 1
		chname(nchain) 		= record(22:22)
		ch1 			= chname(nchain)
	endif
	read(record,2) idum2
	if(idum2.ne.ntest) then
c                if(ntest.ne.-10.and.idum2-ntest.ne.1) then
c			write(6,*) 'Warning : non consecutive residues !!'
c			write(6,*) 'Residues : ',ntest,idum2
c                endif
                ntest = idum2
		nres  = nres +1
	endif
	goto 100
200	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c
c	xmass.f
c
c	Assigns a mass to each atom
c
c ==========================================================================================
c ==========================================================================================
c
	function xmass(name)
c
	real*8 xmass
c
	character*1 name
c
	if(name.eq.'H') then
		xmass = 1.008d0
	elseif(name.eq.'C') then
		xmass = 12.012d0
	elseif(name.eq.'N') then
		xmass = 14.007d0
	elseif(name.eq.'O') then
		xmass = 15.999d0
	elseif(name.eq.'S') then
		xmass = 32.060d0
	elseif(name.eq.'F') then
		xmass = 55.847d0
	elseif(name.eq.'P') then
		xmass = 30.974d0
	else
c		write(6,*) 'Warning: problem with atom :',name
		xmass = 0
	endif
c
	return
	end
