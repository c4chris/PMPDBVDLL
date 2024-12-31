//
//  SharedData.h
//  PDBV
//
//  Created by Christian Iseli on 21.09.2024.
//
#pragma once

#include "PMPDBVDLL.h"

extern "C" {

#ifndef Globals_h
#define Globals_h

#include <stdio.h>

#define kmaxPathChar 1024
#define kmaxNameChar 32
#define maxLayerNb 128

#define kNoCoord 9e99
#define kMaxCharForResNum 4       /* assume that PDB file have residues from 0 and 9999 */

typedef double DBL;
typedef bool Boolean;
#define false 0
#define FALSE 0
#define true 1
#define TRUE 1

/* ----TYPES--------------------------------------------------------------------------- */

//struct Space_Point
//{
//  DBL x, y, z;
//};

//struct FSpace_Point
//{
//  float x, y, z;
//};

//typedef struct Space_Point Point3D;
//typedef Point3D *Point3D_Ptr;
//typedef struct FSpace_Point FPoint3D;

#include "grptypes.h"

typedef struct PDBmodel_struct PDBmodel;
struct PDBmodel_struct
{
  GLYp *GroupsHdl;
  ATOM_TYPEp *AtomNamesHdl;
  char LayerName[kmaxNameChar];
  char originalFilename[kmaxPathChar];
  unsigned long flags;
  unsigned short nbAtomGroups;
  unsigned short nbChains;
};

/* ----------------- Protein Model Part ----------------- */

extern short nbLayer;                         /* nb of molecules loaded */
extern PDBmodel PDB[maxLayerNb];

extern aaDef aaName[];
int myLog(const char*, ...);

#endif /* Globals_h */
}
