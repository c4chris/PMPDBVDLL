//
//  pdb.h
//  PDBV
//
//  Created by Christian Iseli on 22.09.2024.
//
#pragma once

#include "PMPDBVDLL.h"

extern "C" {

#ifndef pdb_h
#define pdb_h

#include <stdio.h>

/* ----PROTOTYPES---------------------------------------------------------------------- */

void Fill_Header(GLYp p, unsigned long flags, unsigned long atomNumStart, char chain, char *name, char *nb, char kind);
int AllocateMemoryForNewLayer(short layer);
short pdb_Input(short whichlayer, FILE *outfile, const char *inName, long filesize, long *openFilePos);

/* buffer to hold atoms coord and type */
/* during pdb file input. */
/* Once all atoms of a group are read, */
/* they are sorted according to the predifined */
/* structures. */

struct aa_struct
{
  Point3D atm_coord;
  char atm_name[kAAnameMaxLength + 1];
  unsigned short Bfactor;
  unsigned short atm_occupancy;
  char altResNum;
  char chain;
  char res_insertion; /* __ts */
};

#endif /* pdb_h */
}
