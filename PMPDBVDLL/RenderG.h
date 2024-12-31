//
//  RenderG.h
//  PDBV
//
//  Created by Christian Iseli on 23.09.2024.
//

#include "PMPDBVDLL.h"

#pragma once
extern "C" {

#ifndef RenderG_h
#define RenderG_h

extern short theLayer;
extern unsigned short theGroup;
extern unsigned short theAtom;
extern unsigned short theAtomB;
extern unsigned short AtomSymbol;
extern GLYp *GroupsHdl;
extern ATOM_TYPEp *AtomNamesHdl;
extern Boolean isCA;

#endif /* RenderG_h */
}
