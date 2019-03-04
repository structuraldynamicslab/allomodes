c ==========================================================================================
c ==========================================================================================
C	toolkit.h	Version 1 1/16/2002	Patrice Koehl
c
c	This include file contains all the dimensions of the
c	arrays needed for all routines in Toolkit
c ==========================================================================================
c ==========================================================================================
c
	integer	nresdef
	integer	natom_type,natrestot
	integer	natsidetot,nangside
	integer	natotv
c
c ==========================================================================================
c	Definitions:
c
c		nresdef		: maximum number of types of residues
c		natom_type	: maximum number of atom types
c		natrestot	: maximum number of atoms in one residue
c       	nangside  	: maximum number of angles in a side chain
c		natsidetot	: maximum number of atoms in a side chain
c		natotv		: maximum number of atoms in a fixed subset of atoms
c
c ==========================================================================================
c
	parameter	(nresdef	= 25)
	parameter	(natom_type	= 30)
	parameter	(natrestot	= 20)
	parameter       (natsidetot	= 14)
	parameter       (nangside  	= 10)
	parameter	(natotv		= 20)
c
