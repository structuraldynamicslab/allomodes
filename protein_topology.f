c ==========================================================================================
c ==========================================================================================
c	Topology.f	Version 1 12/7/1995		Patrice Koehl
c
c	Version 2.0 : 5/26/06	Includes definition of bonds, angles, dihed and improper
c
c	This file contains a blockdata which initialise the name of the 
c	standard amino acids that can be found in a protein, the number of
c	atoms and dihedrals which define each of these amino acids, as well
c	as the name, category and charge of each of these atoms.
c	Only the 20 standard amino acids are considered at this stage
c	Atom category are defined according to CHARMM19 : only polar
c	hydrogens are included
c
c	Atom types, according to CHARMM (used for VdW interactions) :
c
c	Hydrogens :
c			1.  H  (H which can H-bond to neutral atom)
c			2.  HC (H which can H-bond to charged atom)
c			3.  HA (aliphatic hydrogen)
c			4.  HT (TIPS3P water hydrogen)
c			5.  LP (ST2 lone pair)
c	Carbons
c			6.  CT   (aliphatic carbon)
c			7.  C    (carbonyl carbon)
c			8.  CH1E (extended carbon, with 1 hydrogen)
c			9.  CH2E (extended carbon, with 2 hydrogens)
c			10. CH3E (extended carbon, with 3 hydrogens)
c			11. CR1E (ring carbons)
c			12. CM   (C in carbonmonoxide)
c	Nitrogens
c			13. N	(peptide N with no hydrogen atoms attached)
c			14. NR  (N in aromatic ring with no hydrogen atoms)
c			15. NP	(Pyrole N)
c			16. NH1 (Peptide N bound to one hydrogen)
c			17. NH2 (Peptide N bound to two hydrogen atoms)
c			18. NH3 (Peptide N bound to three hydrogen atoms)
c			19. NC2 (Guanidinium N bound to 2 hydrogens)
c
c	Oxygens
c			20. O   (carbonyl O)
c			21. OC  (carboxy O)
c			22. OH1 (hydroxy O)
c			23. OH2 (ST2 water oxygen)
c			24. OM  (O in carbonmonoxide) 
c			25. OT  (TIPS3P water oxygen)
c			26. OS  (esther oxygen)
c
c	Sulfurs
c			27. S   (sulfur)
c			28. SH1E (extended atom S with one hydrogen)
c	
c	Iron
c			29. FE  (iron)
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine topology_new
c
	include 'toolkit.h'
c
	integer		i,ntype
	integer		idihed(nresdef),iatom(nresdef)
	integer		atomtype(nresdef,20)
c
	real*8		atomcharge(nresdef,20)
c
	character	nameatom(nresdef,20)*4,nameres(nresdef)*4
	character	atom_name(natom_type)*4
c
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
	common /name/	 ntype,nameres
	common /atomtyp/ atom_name
c
	data (atom_name(i),i=1,29) /'H   ','HC  ','HA  ','HT  ','LP  ',
     1	'CT  ','C   ','CH1E','CH2E','CH3E','CR1E','CM  ','N   ','NR  ',
     2	'NP  ','NH1 ','NH2 ','NH3 ','NC2 ','O   ','OC  ','OH1 ','OH2 ',
     3  'OM  ','OT  ','OS  ','S   ','SH1E','FE  '/
c
	data ntype/24/
c
	data (nameres(i),i=1,24)/'GLY','ALA','VAL','ILE','LEU','PHE',
     1		'PRO','MET','TRP','CYS','SER','THR','ASN','GLN','TYR',
     2		'HIS','ASP','GLU','LYS','ARG','HSD','HSP','NTER','CTER'/
c
	data (idihed(i),i=1,24)/0,0,1,2,2,2,1,3,2,1,2,2,3,4,2,2,2,3,
     1				5,7,2,2,0,0/
c
	data (iatom(i),i=1,24)/6,6,8,9,9,12,8,9,16,7,8,9,11,12,14,12,9,
     1				10,13,17,12,13,2,1/
c
c ==========================================================================================
c	1. Glycine : no sidechain, only 5 atoms to define the backbone
c	For sake of simplicity, a fictitious atoms is added in place of
c	CB, such that all residues (except NTER and CTER), have 6 atoms
c	in their backbone
c ==========================================================================================
c
	data (nameatom(1,i),i=1,6) /'N','HN','CA','X','C','O'/
	data (atomtype(1,i),i=1,6) /16,1,9,0,7,20/
	data (atomcharge(1,i),i=1,6)/-0.35,0.25,0.10,0,0.55,-0.55/
c
c ==========================================================================================
c	2. Alanine : only CB to define the sidechain
c ==========================================================================================
c
	data (nameatom(2,i),i=1,6) /'N','HN','CA','CB','C','O'/
	data (atomtype(2,i),i=1,6) /16,1,8,10,7,20/
	data (atomcharge(2,i),i=1,6)/-0.35,0.25,0.10,0.,0.55,-0.55/
c
c ==========================================================================================
c	3. Valine : 
c ==========================================================================================
c
	data (nameatom(3,i),i=1,8) /'N','HN','CA','CB','C','O',
     1		'CG1','CG2'/
	data (atomtype(3,i),i=1,8) /16,1,8,8,7,20,10,10/
	data (atomcharge(3,i),i=1,8) /-0.35,0.25,0.10,0,0.55,-0.55,
     1	0.,0./
c
c ==========================================================================================
c	4. Isoleucine :
c ==========================================================================================
c
	data (nameatom(4,i),i=1,9) /'N','HN','CA','CB','C','O','CG1',
     1	'CG2','CD1'/
	data (atomtype(4,i),i=1,9) /16,1,8,8,7,20,9,10,10/
	data (atomcharge(4,i),i=1,9) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0./
c
c ==========================================================================================
c	5. Leucine :
c ==========================================================================================
c
	data (nameatom(5,i),i=1,9) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2'/
	data (atomtype(5,i),i=1,9) /16,1,8,9,7,20,8,10,10/
	data (atomcharge(5,i),i=1,9) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0./
c
c ==========================================================================================
c	6. Phenylalanine :
c ==========================================================================================
c
	data (nameatom(6,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2','CE1','CE2','CZ'/
	data (atomtype(6,i),i=1,12) /16,1,8,9,7,20,7,11,11,11,11,11/
	data (atomcharge(6,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,
     1	-0.55,0.,0.,0.,0.,0.,0./
c
c ==========================================================================================
c	7. Proline :
c ==========================================================================================
c
	data (nameatom(7,i),i=1,8) /'N','X','CA','CB','C','O','CG','CD'/
	data (atomtype(7,i),i=1,8) /13,0,8,9,7,20,9,9/
	data (atomcharge(7,i),i=1,8) /-0.2,0.,0.10,0.,0.55,-0.55,
     1	0.,0.10/
c
c ==========================================================================================
c	8. Methionine
c ==========================================================================================
c
	data (nameatom(8,i),i=1,9) /'N','HN','CA','CB','C','O','CG',
     1	'SD','CE'/
	data (atomtype(8,i),i=1,9) /16,1,8,9,7,20,9,27,10/
	data (atomcharge(8,i),i=1,9) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.06,-0.12,0.06/
c
c ==========================================================================================
c	9. Tryptophane
c ==========================================================================================
c
	data (nameatom(9,i),i=1,16) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2','HE1'/
	data (atomtype(9,i),i=1,16) /16,1,8,9,7,20,7,11,7,16,7,
     1	11,11,11,11,1/
	data (atomcharge(9,i),i=1,16) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	-0.03,0.06,0.10,-0.36,-0.04,-0.03,0.,0.,0.,0.3/
c
c ==========================================================================================
c	10. Cysteine
c ==========================================================================================
c
	data (nameatom(10,i),i=1,7) /'N','HN','CA','CB','C','O','SG'/
	data (atomtype(10,i),i=1,7) /16,1,8,9,7,20,28/
	data (atomcharge(10,i),i=1,7) /-0.35,0.25,0.10,0.19,0.55,
     1	-0.55,-0.19/
c
c ==========================================================================================
c	11. Serine
c ==========================================================================================
c
	data (nameatom(11,i),i=1,8) /'N','HN','CA','CB','C','O','OG',
     1	'HG1'/
	data (atomtype(11,i),i=1,8) /16,1,8,9,7,20,22,1/
	data (atomcharge(11,i),i=1,8) /-0.35,0.25,0.10,0.25,0.55,-0.55,
     1	-0.65,0.40/
c
c ==========================================================================================
c	12. Threonin
c ==========================================================================================
c
	data (nameatom(12,i),i=1,9) /'N','HN','CA','CB','C','O','CG2',
     1	'OG1','HG1'/
	data (atomtype(12,i),i=1,9) /16,1,8,8,7,20,10,22,1/
	data (atomcharge(12,i),i=1,9) /-0.35,0.25,0.10,0.25,0.55,-0.55,
     1	0.,-0.65,0.40/
c
c ==========================================================================================
c	13. Asparagine
c ==========================================================================================
c
	data (nameatom(13,i),i=1,11) /'N','HN','CA','CB','C','O','CG',
     1	'OD1','ND2','HD21','HD22'/
	data (atomtype(13,i),i=1,11) /16,1,8,9,7,20,7,20,17,1,1/
	data (atomcharge(13,i),i=1,11) /-0.35,0.25,0.10,0,0.55,-0.55,
     1	0.55,-0.55,-0.60,0.30,0.30/
c
c ==========================================================================================
c	14. Glutamine
c ==========================================================================================
c
	data (nameatom(14,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'CD','OE1','NE2','HE21','HE22'/
	data (atomtype(14,i),i=1,12) /16,1,8,9,7,20,9,7,20,17,1,1/
	data (atomcharge(14,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.55,-0.55,-0.6,0.3,0.3/
c
c ==========================================================================================
c	15. Tyrosine
c ==========================================================================================
c
	data (nameatom(15,i),i=1,14) /'N','HN','CA','CB','C','O','CG',
     1	'CD1','CD2','CE1','CE2','CZ','OH','HH'/
	data (atomtype(15,i),i=1,14) /16,1,8,9,7,20,7,11,11,11,11,7,
     1	22,1/
	data (atomcharge(15,i),i=1,14)/-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0.,0.,0.,0.25,-0.65,0.40/
c
c ==========================================================================================
c	16. Histidine (ND1 protonated)
c ==========================================================================================
c
	data (nameatom(16,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'ND1','CD2','CE1','NE2','HD1'/
	data (atomtype(16,i),i=1,12) /16,1,8,9,7,20,7,16,11,11,14,1/
	data (atomcharge(16,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.10,-0.40,0.10,0.30,-0.40,0.30/
c
c ==========================================================================================
c	17. Aspartic Acid
c ==========================================================================================
c
	data (nameatom(17,i),i=1,9) /'N','HN','CA','CB','C','O','CG',
     1	'OD1','OD2'/
	data (atomtype(17,i),i=1,9) /16,1,8,9,7,20,7,21,21/
	data (atomcharge(17,i),i=1,9) /-0.35,0.25,0.10,-0.16,0.55,-0.55,
     1	0.36,-0.6,-0.6/
c
c ==========================================================================================
c	18. Glutamic Acid
c ==========================================================================================
c
	data (nameatom(18,i),i=1,10) /'N','HN','CA','CB','C','O','CG',
     1	'CD','OE1','OE2'/
	data (atomtype(18,i),i=1,10) /16,1,8,9,7,20,9,7,21,21/
	data (atomcharge(18,i),i=1,10) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	-0.16,0.36,-0.6,-0.6/
c
c ==========================================================================================
c	19. Lysine
c ==========================================================================================
c
	data (nameatom(19,i),i=1,13) /'N','HN','CA','CB','C','O','CG',
     1  'CD','CE','NZ','HZ1','HZ2','HZ3'/
	data (atomtype(19,i),i=1,13) /16,1,8,9,7,20,9,9,9,18,2,2,2/
	data (atomcharge(19,i),i=1,13) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.,0.,0.25,-0.30,0.35,0.35,0.35/
c
c ==========================================================================================
c	20. Arginine
c ==========================================================================================
c
	data (nameatom(20,i),i=1,17) /'N','HN','CA','CB','C','O','CG',
     1	'CD','NE','CZ','HE','NH1','NH2','HH11','HH12','HH21','HH22'/
	data (atomtype(20,i),i=1,17) /16,1,8,9,7,20,9,9,16,7,1,19,19,
     1	2,2,2,2/
	data (atomcharge(20,i),i=1,17) /-0.35,0.25,0.10,0.,0.55,
     1	-0.55,0.,0.10,-0.4,0.50,0.30,-0.45,-0.45,0.35,0.35,0.35,0.35/
c
c ==========================================================================================
c	21. Histidine (HSD : NE2 protonated)
c ==========================================================================================
c
	data (nameatom(21,i),i=1,12) /'N','HN','CA','CB','C','O','CG',
     1	'ND1','CD2','CE1','NE2','HE2'/
	data (atomtype(21,i),i=1,12) /16,1,8,9,7,20,7,14,11,11,16,1/
	data (atomcharge(21,i),i=1,12) /-0.35,0.25,0.10,0.,0.55,-0.55,
     1	0.10,-0.40,0.10,0.30,-0.40,0.30/
c
c ==========================================================================================
c	22. Histidine (HSC : protonated ND1 and NE2)
c ==========================================================================================
c
	data (nameatom(22,i),i=1,13) /'N','HN','CA','CB','C','O','CG',
     1	'ND1','CD2','CE1','NE2','HD1','HE2'/
	data (atomtype(22,i),i=1,13) /16,1,8,9,7,20,7,14,11,11,16,1,1/
	data (atomcharge(22,i),i=1,13) /-0.35,0.25,0.10,0.10,0.55,-0.55,
     1	0.15,-0.30,0.20,0.45,-0.30,0.35,0.35/
c
c ==========================================================================================
c	23. Nter
c ==========================================================================================
c
	data (nameatom(23,i),i=1,2) /'HT1','HT2'/
	data (atomtype(23,i),i=1,2) /2,2/
	data (atomcharge(23,i),i=1,2) /0.35,0.35/
c
c ==========================================================================================
c	24. Cter
c ==========================================================================================
c
	data (nameatom(24,i),i=1,1) /'OT2'/
	data (atomtype(24,i),i=1,1) /21/
	data (atomcharge(24,i),i=1,1) /-0.57/
c
c ==========================================================================================
c	return to main program
c ==========================================================================================
c
	end
c ==========================================================================================
c ==========================================================================================
c	count_atom_dihed.f		Version 1 21/7/1995	Patrice Koehl
c
c	This subroutine counts the number of atoms and number of dihedral angles
c	for the molecule considered, based on its sequence
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine count_atom_dihed(natom,ndihed)
c
c ==========================================================================================
c	toolkit.h contains the dimension of the arrays
c ==========================================================================================
c
	include 'toolkit.h'
c
	integer		i,nseq,ndihed,natom,nchain
	integer		itype(*)
	integer		idihed(nresdef),iatom(nresdef)
	integer		atomtype(nresdef,20)
	integer		chain(*)
c
	real*8		atomcharge(nresdef,20)
c
	character 	seq(*)*4
	character	nameatom(nresdef,20)*4
c
	pointer	(ptr_itype,itype)
	pointer	(ptr_chain,chain)
	pointer	(ptr_seq,seq)
c
c ==========================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
c ==========================================================================================
c
	natom  = 0
	ndihed = 0
c
c ==========================================================================================
c	sum contribution of each residue
c ==========================================================================================
c
	do 100 i = 1,nseq
c
		natom = natom +  iatom(itype(i))
		if(itype(i).ne.23.and.itype(i).ne.24) then
			ndihed = ndihed + 3+idihed(itype(i))
		else
			ndihed = ndihed + idihed(itype(i))
		endif
c
100	continue
c
c ==========================================================================================
c	Add interchain degrees of freedom
c ==========================================================================================
c
	ndihed = ndihed + 6*(nchain-1)
c
c ==========================================================================================
c	Remove last omega of each chain
c ==========================================================================================
c
	ndihed = ndihed - nchain
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	Prepare_protein.f		Version 1 21/7/1995	Patrice Koehl
c
c	This subroutine prepares the molecule : 
c		- gives name of the atoms
c		- gives number of dihedral angles, number of atoms
c		  (total, and per residue)
c		- build up exclusion lists
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine prepare_protein
c
c ==========================================================================================
c	toolkit.h contains the dimension of the arrays
c ==========================================================================================
c
	include 'toolkit.h'
c
	integer		i,j,nseq,ndihed,natom,nchain
	integer		itype(*),listdihed(*)
	integer		idihed(nresdef),iatom(nresdef)
	integer		atomtype(nresdef,20),listatom(*)
	integer		chain(*),occup(*)
c
	real*8		atomcharge(nresdef,20)
	real*8		dihed(*),coord(*)
c
	character 	atom(*)*4,seq(*)*4
	character	nameatom(nresdef,20)*4
c
	pointer	(ptr_itype,itype)
	pointer	(ptr_listdihed,listdihed)
	pointer	(ptr_listatom,listatom)
	pointer	(ptr_chain,chain)
	pointer	(ptr_occup,occup)
	pointer (ptr_dihed,dihed)
	pointer	(ptr_coord,coord)
	pointer	(ptr_atom,atom)
	pointer	(ptr_seq,seq)
c
c ==========================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /angle/   ptr_dihed,ptr_listdihed,ndihed
	common /xyz/	 ptr_coord,ptr_occup,ptr_listatom,natom
	common /names/	 ptr_atom
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
c ==========================================================================================
c
	ndihed = 0
	natom  = 0
c
c ==========================================================================================
	do 200 i = 1,nseq
c
c 	====================================================================================
c		1. Store all info for each residue :
c 	====================================================================================
c
		do 100 j = 1,iatom(itype(i))
			atom(natom+j)   = nameatom(itype(i),j)
100		continue
c
		listatom(i)  = iatom(itype(i))
		if(itype(i).ne.23.and.itype(i).ne.24) then
			listdihed(i) = 3+idihed(itype(i))
		else
			listdihed(i) = idihed(itype(i))
		endif
c
		ndihed       = ndihed + listdihed(i)
		natom        = natom + listatom(i)
c
200	continue
c
	ndihed = ndihed + 6*(nchain-1) - nchain
c
	do 300 i = 2,nseq-1
		if(itype(i+1).eq.24) listdihed(i)=listdihed(i)-1
300	continue
c
c	write(6,*) ' '
c	write(6,*) ' Number of residue              : ',nseq-2
c	write(6,*) ' Number of atoms (with polar H) : ',natom
c	write(6,*) ' '
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Topology_dihed.f		Version 1 9/26/2006	Patrice Koehl
c
c	THis blockdata defines the dihedral angle "structure" of each amino acid
c	sidechain.
c
c	It defines:
c
c		- bondtopo:	for each dihedral angle, gives the four atoms
c				involved
c		- bondup:	for each dihedral angle, number of dihedral angles
c				that branches from it
c		- bondatm	number of atoms that depends specifically on one given
c				dihedral angle
c		- bondlist	list of atoms that depends on one given dihedral angle
c ==========================================================================================
c ==========================================================================================
c
	subroutine topology_dihed
c
	include 'toolkit.h'
c
	integer i
	integer	bondtopo(4,nangside,nresdef)
	integer	bondup(nangside,nresdef)
	integer	bondatm(nangside,nresdef)
	integer	bondlist(20,nangside,nresdef)
c
	common /topo_dihed/ bondtopo,bondup,bondatm,bondlist
c
c ==========================================================================================
c	1. Gly and 2. Alanine
c ==========================================================================================
c
c	No sidechain dihedral angle
c
c ==========================================================================================
c	3. Val
c ==========================================================================================
c
	data (bondtopo(i,1,3),i=1,4)/1,3,4,7/
	data (bondup(i,3),i=1,1)/0/
	data (bondatm(i,3),i=1,1)/2/
	data (bondlist(i,1,3),i=1,2)/7,8/
c
c ==========================================================================================
c	4. Ile
c ==========================================================================================
c
	data (bondtopo(i,1,4),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,4),i=1,4)/3,4,7,9/
	data (bondup(i,4),i=1,2)/1,0/
	data (bondatm(i,4),i=1,2)/2,1/
	data (bondlist(i,1,4),i=1,2)/7,8/
	data (bondlist(i,2,4),i=1,1)/9/
c
c ==========================================================================================
c	5. Leu
c ==========================================================================================
c
	data (bondtopo(i,1,5),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,5),i=1,4)/3,4,7,8/
	data (bondup(i,5),i=1,2)/1,0/
	data (bondatm(i,5),i=1,2)/1,2/
	data (bondlist(i,1,5),i=1,1)/7/
	data (bondlist(i,2,5),i=1,2)/8,9/
c
c ==========================================================================================
c	6. Phe
c ==========================================================================================
c
	data (bondtopo(i,1,6),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,6),i=1,4)/3,4,7,8/
	data (bondup(i,6),i=1,2)/1,0/
	data (bondatm(i,6),i=1,2)/1,5/
	data (bondlist(i,1,6),i=1,1)/7/
	data (bondlist(i,2,6),i=1,5)/8,9,10,11,12/
c
c ==========================================================================================
c	7. Pro
c ==========================================================================================
c
	data (bondtopo(i,1,7),i=1,4)/1,3,4,7/
	data (bondup(i,7),i=1,1)/0/
	data (bondatm(i,7),i=1,1)/2/
	data (bondlist(i,1,7),i=1,2)/7,8/
c
c ==========================================================================================
c	8. Met
c ==========================================================================================
c
	data (bondtopo(i,1,8),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,8),i=1,4)/3,4,7,8/
	data (bondtopo(i,3,8),i=1,4)/4,7,8,9/
	data (bondup(i,8),i=1,3)/1,1,0/
	data (bondatm(i,8),i=1,3)/1,1,1/
	data (bondlist(i,1,8),i=1,1)/7/
	data (bondlist(i,2,8),i=1,1)/8/
	data (bondlist(i,3,8),i=1,1)/9/
c
c ==========================================================================================
c	9. Trp
c ==========================================================================================
c
	data (bondtopo(i,1,9),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,9),i=1,4)/3,4,7,8/
	data (bondup(i,9),i=1,2)/1,0/
	data (bondatm(i,9),i=1,2)/1,9/
	data (bondlist(i,1,9),i=1,1)/7/
	data (bondlist(i,2,9),i=1,9)/8,9,10,11,12,13,14,15,16/
c
c ==========================================================================================
c	10 Cysteine
c ==========================================================================================
c
	data (bondtopo(i,1,10),i=1,4)/1,3,4,7/
	data (bondup(i,10),i=1,1)/0/
	data (bondatm(i,10),i=1,1)/1/
	data (bondlist(i,1,10),i=1,1)/7/
c
c ==========================================================================================
c	11 Ser
c ==========================================================================================
c
	data (bondtopo(i,1,11),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,11),i=1,4)/3,4,7,8/
	data (bondup(i,11),i=1,2)/1,0/
	data (bondatm(i,11),i=1,2)/1,1/
	data (bondlist(i,1,11),i=1,1)/7/
	data (bondlist(i,2,11),i=1,1)/8/
c
c ==========================================================================================
c	12. Thr
c ==========================================================================================
c
	data (bondtopo(i,1,12),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,12),i=1,4)/3,4,8,9/
	data (bondup(i,12),i=1,2)/1,0/
	data (bondatm(i,12),i=1,2)/2,1/
	data (bondlist(i,1,12),i=1,2)/7,8/
	data (bondlist(i,2,12),i=1,1)/9/
c
c ==========================================================================================
c	13. Asn
c ==========================================================================================
c
	data (bondtopo(i,1,13),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,13),i=1,4)/3,4,7,8/
	data (bondtopo(i,3,13),i=1,4)/4,7,9,10/
	data (bondup(i,13),i=1,3)/1,1,0/
	data (bondatm(i,13),i=1,3)/1,2,2/
	data (bondlist(i,1,13),i=1,1)/7/
	data (bondlist(i,2,13),i=1,2)/8,9/
	data (bondlist(i,3,13),i=1,2)/10,11/
c
c ==========================================================================================
c	14. Gln
c ==========================================================================================
c
	data (bondtopo(i,1,14),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,14),i=1,4)/3,4,7,8/
	data (bondtopo(i,3,14),i=1,4)/4,7,8,9/
	data (bondtopo(i,4,14),i=1,4)/7,8,10,11/
	data (bondup(i,14),i=1,4)/1,1,1,0/
	data (bondatm(i,14),i=1,4)/1,1,2,2/
	data (bondlist(i,1,14),i=1,1)/7/
	data (bondlist(i,2,14),i=1,1)/8/
	data (bondlist(i,3,14),i=1,2)/9,10/
	data (bondlist(i,4,14),i=1,2)/11,12/
c
c ==========================================================================================
c	15. Tyr
c ==========================================================================================
c
	data (bondtopo(i,1,15),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,15),i=1,4)/3,4,7,8/
	data (bondup(i,15),i=1,2)/1,0/
	data (bondatm(i,15),i=1,2)/1,7/
	data (bondlist(i,1,15),i=1,1)/7/
	data (bondlist(i,2,15),i=1,7)/8,9,10,11,12,13,14/
c
c ==========================================================================================
c	16. His
c ==========================================================================================
c
	data (bondtopo(i,1,16),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,16),i=1,4)/3,4,7,8/
	data (bondup(i,16),i=1,2)/1,0/
	data (bondatm(i,16),i=1,2)/1,5/
	data (bondlist(i,1,16),i=1,1)/7/
	data (bondlist(i,2,16),i=1,5)/8,9,10,11,12/
c
c ==========================================================================================
c	17. Asp
c ==========================================================================================
c
	data (bondtopo(i,1,17),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,17),i=1,4)/3,4,7,8/
	data (bondup(i,17),i=1,2)/1,0/
	data (bondatm(i,17),i=1,2)/1,2/
	data (bondlist(i,1,17),i=1,1)/7/
	data (bondlist(i,2,17),i=1,2)/8,9/
c
c ==========================================================================================
c	18. Glu
c ==========================================================================================
c
	data (bondtopo(i,1,18),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,18),i=1,4)/3,4,7,8/
	data (bondtopo(i,3,18),i=1,4)/4,7,8,9/
	data (bondup(i,18),i=1,3)/1,1,0/
	data (bondatm(i,18),i=1,3)/1,1,2/
	data (bondlist(i,1,18),i=1,1)/7/
	data (bondlist(i,2,18),i=1,1)/8/
	data (bondlist(i,3,18),i=1,2)/9,10/
c
c ==========================================================================================
c	19. Lys
c ==========================================================================================
c
	data (bondtopo(i,1,19),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,19),i=1,4)/3,4,7,8/
	data (bondtopo(i,3,19),i=1,4)/4,7,8,9/
	data (bondtopo(i,4,19),i=1,4)/7,8,9,10/
	data (bondtopo(i,5,19),i=1,4)/8,9,10,11/
	data (bondup(i,19),i=1,5)/1,1,1,1,0/
	data (bondatm(i,19),i=1,5)/1,1,1,1,3/
	data (bondlist(i,1,19),i=1,1)/7/
	data (bondlist(i,2,19),i=1,1)/8/
	data (bondlist(i,3,19),i=1,1)/9/
	data (bondlist(i,4,19),i=1,1)/10/
	data (bondlist(i,5,19),i=1,3)/11,12,13/
c
c ==========================================================================================
c	20. Arg
c ==========================================================================================
c
	data (bondtopo(i,1,20),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,20),i=1,4)/3,4,7,8/
	data (bondtopo(i,3,20),i=1,4)/4,7,8,9/
	data (bondtopo(i,4,20),i=1,4)/7,8,9,10/
	data (bondtopo(i,5,20),i=1,4)/8,9,10,12/
	data (bondtopo(i,6,20),i=1,4)/9,10,12,14/
	data (bondtopo(i,7,20),i=1,4)/9,10,13,16/
	data (bondup(i,20),i=1,7)/1,1,1,1,1,1,0/
	data (bondatm(i,20),i=1,7)/1,1,1,2,2,2,2/
	data (bondlist(i,1,20),i=1,1)/7/
	data (bondlist(i,2,20),i=1,1)/8/
	data (bondlist(i,3,20),i=1,1)/9/
	data (bondlist(i,4,20),i=1,2)/10,11/
	data (bondlist(i,5,20),i=1,2)/12,13/
	data (bondlist(i,6,20),i=1,2)/14,15/
	data (bondlist(i,7,20),i=1,2)/16,17/
c
c ==========================================================================================
c	21. Hse
c ==========================================================================================
c
	data (bondtopo(i,1,21),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,21),i=1,4)/3,4,7,8/
	data (bondup(i,21),i=1,2)/1,0/
	data (bondatm(i,21),i=1,2)/1,5/
	data (bondlist(i,1,21),i=1,1)/7/
	data (bondlist(i,2,21),i=1,5)/8,9,10,11,12/
c
c ==========================================================================================
c	22. Hsp
c ==========================================================================================
c
	data (bondtopo(i,1,22),i=1,4)/1,3,4,7/
	data (bondtopo(i,2,22),i=1,4)/3,4,7,8/
	data (bondup(i,22),i=1,2)/1,0/
	data (bondatm(i,22),i=1,2)/1,6/
	data (bondlist(i,1,22),i=1,1)/7/
	data (bondlist(i,2,22),i=1,6)/8,9,10,11,12,13/
c
c ==========================================================================================
c	23. Nter and 24. Cter
c ==========================================================================================
c
c	No sidechain
c
	end
c ==========================================================================================
c ==========================================================================================
c	Side_dihed_info.f	Version 1 9/26/3006		Patrice Koehl
c
c	This subroutine prepares the side chains for the symbolic representa-
c	tion of the protein :
c	it defines the units V(a) : sets of atoms whose relative positions
c	within each of them remain fixed for any conformational changes of
c	the molecule. Each atom of V(a) depend only on the dihedral angle
c	D(a-1), when all previous dihedral angles are kept fixed.
c
c	Input of the subroutine :
c				- itype	: residue type considered (from 1 to 20)
c				- nsub	: max number of possible subsets
c					  within a side chain
c				- nat	: max number of atom in one subset
c
c	Information from common blocks
c
c				- bondtopo	: array which contains the atom
c					  involved for each dihedral angle
c					  for each residue
c				- bondup	: for each dihedral angle, indicate
c					  # of following dihedral angles
c				- bondatm	: # of atoms in each subset of atoms
c				- bondlist : list of atoms in a subset
c
c	Output of the subroutine :
c				- nside : number of subunits V(a) in the
c				           side chain of the residue
c				- nlist : number of atoms in each subset
c					   (array)
c				- list	: list of atoms in each subset
c					  (array)
c				- bond	: atoms that define the dihedral
c					  angles in the sidechain (array)
c				- p	: number of dihedral angle that 
c					  branch from a given dihedral
c					  angle
c ==========================================================================================
c ==========================================================================================
c
	subroutine side_dihed_info(itype,nsub,nat,nside,nlist,
     1		list,bond,p,ps1)
c
	include 'toolkit.h'
c
	integer	i,itype,nsub,nat
	integer	nside,ps1,j
	integer	nlist(nsub),list(nat,nsub)
	integer	bond(4,nsub),p(nsub)
	integer	bondtopo(4,nangside,nresdef)
	integer	bondup(nangside,nresdef)
	integer	bondatm(nangside,nresdef)
	integer iatom(nresdef),idihed(nresdef)
	integer	bondlist(20,nangside,nresdef)
	integer	atomtype(nresdef,20)
c
	real*8	atomcharge(nresdef,20)
c
	character nameatom(nresdef,20)*4
c
c ==========================================================================================
	common /topo_dihed/ bondtopo,bondup,bondatm,bondlist
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
c ==========================================================================================
c
	nside = idihed(itype)
c
	ps1 = bondup(1,itype)
c
	do 300 i = 1,nside
		nlist(i) = bondatm(i,itype)
		do 100 j = 1,nlist(i)
			list(j,i) = bondlist(j,i,itype)
100		continue
		if(i.ge.2) then
			p(i-1) = bondup(i,itype)
			do 200 j = 1,4
				bond(j,i-1) = bondtopo(j,i,itype)
200			continue
		endif
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Prepare_dihed.f		Version 1 9/26/2006	Patrice Koehl
c
c	This subroutine creates a symbolic representation of the molecule :
c	For each dihedral angle D(a), we define an atom set V(a+1) which
c	contains the atoms whose coordinates only depend on D(a), when all 
c	previous D(b), b from 1 to a-1, are kept fixed.
c	Also, for each dihedral angle D(a), we define 5 numbers :
c	bond(j,a) for j = 1 to 4 gives the atoms that define a ; also
c	p(a) defines the number of dihedral angles that branches from D(a).
c	Finally, for each atom, we indicate the subset V(a) to which it belongs
c
c	Two parameters are needed for prepare_dihed:
c			-natotv		: maximum number of atom in a subset
c			- nangside	: maximum number of subsets in a side
c					  chain
c	Two new common blocks are defined :
c			/defda/		: contains information on each dihedral
c					  angle : list of atoms that define the
c					  rotatable bond, and number of angles
c				 	  that branches from it
c			/defva/		: for each atom, gives the subset V(a)
c					  that contains it. (array listva)
c ==========================================================================================
c ==========================================================================================
c
	subroutine prepare_dihed
c
c ==========================================================================================
c ==========================================================================================
c
	include 'toolkit.h'
c
	integer		i,j,k,l
	integer		nside,nat,nv,nref,nchain,ndi
	integer		ps1
	integer		nseq,natom,nangtot
c
	integer		itype(*)
	integer		listatom(*),chain(*)
	integer		atomtype(nresdef,20)
	integer		check(*)
	integer		bond(4,*)
	integer		p(*),nlist(*),list(natotv,*)
	integer		bonds(4,nangside),ps(nangside)
	integer		nlists(nangside),lists(natotv,nangside)
	integer		listva(*),typevar(*)
	integer		idihed(nresdef),iatom(nresdef)
c
	real*8		coord(*)
	real*8		atomcharge(nresdef,20)
c
	character	nameatom(nresdef,20)*4,seq(*)*4
	character	atom(*)*4,dihed_names(*)*6
c
	pointer	(ptr_itype,itype)
	pointer (ptr_chain,chain)
	pointer	(ptr_coord,coord)
	pointer (ptr_check,check)
	pointer	(ptr_listatom,listatom)
	pointer	(ptr_bond,bond)
	pointer	(ptr_p,p)
	pointer	(ptr_listva,listva)
	pointer	(ptr_nlist,nlist)
	pointer	(ptr_list,list)
	pointer	(ptr_seq,seq)
	pointer	(ptr_atom,atom)
	pointer	(ptr_typevar,typevar)
	pointer	(ptr_dihed_names,dihed_names)
c
c ==========================================================================================
	common  /default/	nameatom,atomtype,atomcharge,idihed,
     1				iatom
	common /protein/ 	ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /xyz/     	ptr_coord,ptr_check,ptr_listatom,natom
	common /defda/   	ptr_typevar,ptr_bond,ptr_p,nref,nv
	common /defva/	 	ptr_listva
	common /names/	 	ptr_atom
	common /angle_names/	ptr_dihed_names
c ==========================================================================================
c
	nangtot = 0
	do 100 i = 1,nseq
		nangtot = nangtot + 3 + idihed(itype(i))
100	continue
	nangtot = nangtot + 6*(nchain-1)
c
	ptr_typevar 	= malloc(4*nangtot)
	ptr_bond 	= malloc(4*4*nangtot)
	ptr_p    	= malloc(4*nangtot)
	ptr_listva 	= malloc(4*natom)
	ptr_nlist  	= malloc(4*nangtot)
	ptr_list   	= malloc(4*natotv*nangtot)
	ptr_dihed_names = malloc(6*nangtot)
c
	nat = 0
	nv  = 0
	ndi = 0
c
	do 900 i = 1,nseq
c
c ==========================================================================================
c	For each residue : define first V(a) ( HN and Calpha) and
c	following D(a) : PHI. In the first residue, N is included
c	in this subset. (this V(a) is in fact defined according to the
c	preceding OMEGA (remember that V(a) follows D(a-1)))
c ==========================================================================================
c
		if(itype(i).eq.23.or.itype(i).eq.24) then
			nat = nat + listatom(i)
			goto 900
		endif
c
		nv = nv + 1
		if(itype(i-1).eq.23) then
			nlist(nv) = listatom(i-1)+3
			do 200 j = 1,listatom(i-1)+3
				list(j,nv) = nat - listatom(i-1) + j
200			continue
			bond(1,nv) = nat-listatom(i-1)+1
		else
			nlist(nv)  = 2
			list(1,nv) = 2 + nat
			list(2,nv) = 3 + nat
			bond(1,nv) = 5 + nat - listatom(i-1)
		endif
		bond(2,nv) = 1 + nat
		bond(3,nv) = 3 + nat
		bond(4,nv) = 5 + nat
c
		ndi = ndi + 1
		dihed_names(ndi) = 'PHI   '
c
c ==========================================================================================
c	Now define V(a) corresponding to PHI : it contains CO
c	and Cbeta. D(a) will either be PSI in the absence of side chain
c	or CHI1 in the presence of a side chain. In the later case, p(PHI)
c	is set to 2 since the molecule branches at this point.
c ==========================================================================================
c
		nv = nv + 1
		nlist(nv)  = 2
		list(1,nv) = 4 + nat
		list(2,nv) = 5 + nat
c
c ==========================================================================================
c	At this stage proline is kept planar so that its Cgamma is included
c	in this subset :
c ==========================================================================================
c
		if(itype(i).eq.7) then
			nlist(nv) = 4
			do 300 j = 3,4
				list(j,nv) = j+4+nat
300			continue
		endif
c
		call side_dihed_info(itype(i),nangside,natotv,nside,
     1			nlists,lists,bonds,ps,ps1)
c
		if(nside.eq.0) then
			p(nv-1)    = 1
			bond(1,nv) = 1 +nat
			bond(2,nv) = 3+nat
			bond(3,nv) = 5+nat
			bond(4,nv) = 1 + nat + listatom(i)
			p(nv)      = 1
c
			ndi = ndi + 1
			dihed_names(ndi) = 'PSI   '
c
c ==========================================================================================
c	If a side chain exits, it is incorporated at this level :
c ==========================================================================================
c
		else
			p(nv-1)    = 2
			bond(1,nv) = 1+nat
			bond(2,nv) = 3+nat
			bond(3,nv) = 4+nat
			bond(4,nv) = 7+nat
			p(nv) = ps1
			ndi = ndi + 1
			dihed_names(ndi) = 'CHI1  '
			do 600 k = 1,nside-1
				nv = nv + 1
				nlist(nv) = nlists(k)
				do 400 l = 1,nlists(k)
					list(l,nv) = lists(l,k)+nat
400				continue
				do 500 l = 1,4
					bond(l,nv) = bonds(l,k)+nat
500				continue
				p(nv) = ps(k)
				ndi = ndi + 1
				dihed_names(ndi) = 'CHI   '
				write(dihed_names(ndi)(4:4),'(i1)') k+1
600			continue
c
			nv = nv + 1
			nlist(nv) = nlists(nside)
			do 700 l = 1,nlists(nside)
				list(l,nv) = lists(l,nside) +nat
700			continue
			bond(1,nv) = 1+nat
			bond(2,nv) = 3+nat
			bond(3,nv) = 5+nat
			bond(4,nv) = 1+nat+listatom(i)
			p(nv)      = 1
c
			ndi = ndi + 1
			dihed_names(ndi) = 'PSI   '
		endif
c
c ==========================================================================================
c	Finally, if we have not reached the last residue, we define
c	V(a) for PSI, which contains O and N from the next residue.
c	Then D(a) (i.e. OMEGA) is defined, before cycling to the next
c	residue. For the last residue, p(PSI) is set to zero and
c	 a V(end) is defined out of the loop.
c ==========================================================================================
c
		if(itype(i+1).eq.24) then
			p(nv) = 0
			nv = nv + 1
			nlist(nv) = 1+listatom(i+1)
			list(1,nv) = 6+nat
			bond(1,nv) = 0
			bond(2,nv) = 0
			bond(3,nv) = 0
			bond(4,nv) = 0
			do 800 l = 1,listatom(i+1)
				list(l+1,nv) = nat+listatom(i)+l
800			continue
			ndi = ndi + 1
			dihed_names(ndi) = 'TRANSX'
			if(i.ne.nseq-1) then
				p(nv-1) = 1
				p(nv) = 1
				do 850 l = 1,5
					nv = nv + 1
					nlist(nv) = 0
					bond(1,nv) = 0
					bond(2,nv) = 0
					bond(3,nv) = 0
					bond(4,nv) = 0
					p(nv) = 1
850				continue
				ndi = ndi + 1
				dihed_names(ndi) = 'TRANSY'
				ndi = ndi + 1
				dihed_names(ndi) = 'TRANSZ'
				ndi = ndi + 1
				dihed_names(ndi) = 'EULER1'
				ndi = ndi + 1
				dihed_names(ndi) = 'EULER2'
				ndi = ndi + 1
				dihed_names(ndi) = 'EULER3'
			endif
			nat = nat + listatom(i)
			goto 900
		endif
c
		nv         = nv + 1
		nlist(nv)  = 2
		list(1,nv) = 6 + nat
		bond(1,nv) = 3+nat
		bond(2,nv) = 5+nat
c
		nat        = nat + listatom(i)
		list(2,nv) = 1+nat
		bond(3,nv) = 1+nat
		bond(4,nv) = 3+nat
		p(nv)      = 1
		ndi = ndi + 1
		dihed_names(ndi) = 'OMEGA '
c
900	continue
c
	if(ndi.ne.nv) then
		write(6,*) 'Problem!!!'
		write(6,*) 'ndi,nv:',ndi,nv
		stop
	endif
c
	nref = nv - 1
c
c ==========================================================================================
c	now define for each atom the subset V(a) to which it belongs
c ==========================================================================================
c
	do 1100 i = 1,nv
		do 1000 j = 1,nlist(i)
			listva(list(j,i)) = i
1000		continue
1100	continue
c
	typevar(1) = 1
	ps1 = 0
	do 1200 i = 2,nv-1
		if(nlist(i).eq.0.or.nlist(i+1).eq.0) then
			ps1 = ps1 + 1
			typevar(i) = 1 + ps1
		else
			ps1 = 0
			typevar(i) = 1
		endif
1200	continue
c
	call free(ptr_nlist)
	call free(ptr_list)
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Buildmap.f		Version 1 8/13/1990		Patrice Koehl
c
c	This subroutine builds up a map which will be needed for derivative
c	calculations.
c	It builds a matrix whose rows are the dihedral angle and column the
c	atom subsets V(a). A 1 on position (i,j) means that V(j) is on the
c	movable side of D(i). ( The molecule is constructed from residu 1
c	such that anything on the left of a dihedral angle D(a) is non
c	movable for this dihedral angle)
c	This program is based on a multistack algorithm : it builds up 
c	the molecule from its end part, taking in account its branching
c	by creating a new side chain at each time it is needed ;
c	
c ==========================================================================================
c ==========================================================================================
c
	subroutine buildmap(nangtot,list,ma)
c
c ==========================================================================================
c ==========================================================================================
c
	include 'toolkit.h'
c
	integer maxchain
	parameter (maxchain=5)
c
	integer		nangtot
c
	integer*1	ma(nangtot,*)
c
	integer	i,j,k,iang
	integer	nref,nv,bond(4,*),p(*)
	integer	nval(maxchain),list(maxchain,*)
	integer	typevar(*)
c
	pointer	(ptr_bond,bond)
	pointer	(ptr_p,p)
	pointer	(ptr_typevar,typevar)
c
c ==========================================================================================
	common /defda/ ptr_typevar,ptr_bond,ptr_p,nref,nv
c ==========================================================================================
c
	do 200 i = 1,nref
		do 100 j = 1,nv
			ma(i,j) = 0
100		continue
200	continue
c
	do 300 i = 1,maxchain
		nval(i) = 0
300	continue
c
	iang = nref
c
	do 1300 i = nv-1,1,-1
c
c		if(p(i).eq.-1) goto 1300
c
		if(p(i).eq.0) then
			do 400 j = 1,maxchain
				if(nval(j).eq.0) goto 500
400			continue
500			nval(j) = 1
			list(j,1) = i+1
			ma(iang,i+1) = 1
		elseif(p(i).eq.1) then
			do 600 j = maxchain,1,-1
				if(nval(j).ne.0) goto 700
600			continue
700			nval(j) = nval(j) + 1
			list(j,nval(j)) = i+1
			do 800 k = 1,nval(j)
				ma(iang,list(j,k)) = 1
800			continue
		elseif(p(i).eq.2) then
			do 900 j = maxchain,1,-1
				if(nval(j).ne.0) goto 1000
900			continue
1000			do 1100 k = 1,nval(j)
				nval(j-1) = nval(j-1) + 1
				list(j-1,nval(j-1)) = list(j,k)
1100			continue
			nval(j-1) = nval(j-1) + 1
			list(j-1,nval(j-1)) = i+1	
			do 1200 k = 1,nval(j-1)
				ma(iang,list(j-1,k)) = 1
1200			continue
			nval(j) = 0
		endif
c
		iang = iang - 1
c
1300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Atompos.f		Version 1 24/7/1995	Patrice Koehl
c
c	This subroutine, knowing the name and the residu number of an
c	atom, gives its position in the array of atoms build for the 
c	calculation
c	Input :
c		- ires	: residu which contains the atom
c		- name	: name of the atom
c	Output :
c		- ipos	: atom number
c ==========================================================================================
c ==========================================================================================
c
	subroutine atompos(ires,name,ipos,iadd)
c
	include 'toolkit.h'
c
	integer		i,j,nseq,nchain,natom,ires,ipos,iadd
	integer		itype(*),listatom(*)
	integer		iatom(nresdef),idihed(nresdef)
	integer		atomtype(nresdef,20)
	integer		chain(*),occup(*)
c
	real*8		coord(*)
	real*8		atomcharge(nresdef,20)
c
	character 	name*4,seq(*)*4
	character	nameatom(nresdef,20)*4
c
	pointer	(ptr_itype,itype)
	pointer	(ptr_chain,chain)
	pointer	(ptr_seq,seq)
	pointer	(ptr_coord,coord)
	pointer	(ptr_occup,occup)
	pointer	(ptr_listatom,listatom)
c
c ==========================================================================================
	common /protein/ ptr_itype,ptr_chain,ptr_seq,nseq,nchain
	common /xyz/     ptr_coord,ptr_occup,ptr_listatom,natom
	common /default/ nameatom,atomtype,atomcharge,idihed,iatom
c ==========================================================================================
c
	ipos = 0
	do 10 i = 1,ires - 1
		ipos = ipos + listatom(i)
10	continue
c
        iadd = 0
c
	if(name.eq.'H   ') name = 'HN'
c
	if(name.eq.'HT3') name = 'HN'
c
	if(name.eq.'HT1') then
		iadd = -1
		goto 40
	endif
	if(name.eq.'HT2') then
		iadd = 0
		goto 40
	endif
c
	if(name.eq.'OT1 ') name = 'O'
c
	if(name.eq.'OT2') then
		iadd = 1 + listatom(ires)
		goto 40
	endif
c
	if(name.eq.'OXT') then
		iadd = 1 + listatom(ires)
		goto 40
	endif
c
	if(itype(ires).eq.4.and.name.eq.'CD') name = 'CD1'
c
	if(itype(ires).eq.11.and.name.eq.'HG') name = 'HG1'
c
	j = itype(ires)
	do 20 i = 1,iatom(j)
		if(name.eq.nameatom(j,i)) then
			iadd = i
			goto 40
		endif
20	continue
c
40	continue
c
	ipos = ipos + iadd
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Residue.f		Version 1 24/7/1992	Patrice koehl
c
c       This functions gives a number (from 1 to ntype) to each amino acid
c
c ==========================================================================================
c ==========================================================================================
c
        integer function residue(aa)
c
	include 'toolkit.h'
c
        character*4     aa
        character*4     nameres(nresdef)
c
        integer i,ntype
c
        common /name/   ntype,nameres
c
	if(aa.eq.'NTER') then
		residue = 23
		return
	endif
c
	if(aa.eq.'CTER') then
		residue = 24
		return
	endif
c
        do 100 i = 1,ntype
                if(aa.eq.nameres(i)) then
                        residue = i
                        return
                endif
100     continue
c
c	By default, set residue as Gly
c
	residue = 1
c
        return
	end
