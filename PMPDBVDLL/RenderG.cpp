//
//  RenderG.c
//  PDBV
//
//  Created by Christian Iseli on 23.09.2024.
//

#include "pch.h"
#include <stdio.h>
#include "Globals.h"
#include "RenderG.h"

/* ----GLOBALS------------------------------------------------------------------------- */

short theLayer;
unsigned short theGroup;
unsigned short theAtom;
unsigned short theAtomB;
unsigned short AtomSymbol;
GLYp *GroupsHdl;
ATOM_TYPEp *AtomNamesHdl;
Boolean isCA;
