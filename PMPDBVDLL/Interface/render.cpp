//
//  render.c
//  PDBV
//
//  Created by Christian Iseli on 23.09.2024.
//

#include "pch.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Memory.h>

#include "PMPDBVDLL.h"
#include "Globals.h"
#include "utils.h"
#include "render.h"
#include "RenderG.h"

static void Bond(Point3D A, Point3D B);
static void Atom(Point3D A);
static Boolean checkIfPresent(Point3D A);

extern void GL_atm(FPoint3D A);
extern void GL_bnd(FPoint3D A, FPoint3D B);

static void Bond(Point3D A, Point3D B)
{
  FPoint3D FA = {(float)A.x, (float)A.y, (float)A.z};
  FPoint3D FB = {(float)B.x, (float)B.y, (float)B.z};
  GL_bnd(FA, FB);
} /* Bond */

/* ------------------------------------------------------------------------------------ */
static void Atom(Point3D A)
{
  FPoint3D FA = {(float)A.x, (float)A.y, (float)A.z};
  GL_atm(FA);
} /* Atom */

/* ------------------------------------------------------------------------------------ */

short Draw_BackBone(GLYp p)
{
    AtomSymbol = 7; /*AtomKind = 'N';*/
    theAtom = 0;
    Atom((*p).N);
    AtomSymbol = 6; /*AtomKind = 'C';*/
    /* CA-C */
    theAtom = 2;
    Atom((*p).C);
    isCA = true;
    theAtom = 1;
    Atom((*p).CA);
    isCA = false;
    AtomSymbol = 8; /*AtomKind = 'O';*/
    theAtom = 3;
    Atom((*p).O);
    return 0;
} /* Draw_BackBone */

/* ------------------------------------------------------------------------------------ */

void Draw_ALA(ALAp p)
{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);

} /* Draw_ALA */

/* ------------------------------------------------------------------------------------ */

void Draw_LYS(LYSp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD);
  theAtom = 7;
  Atom((*p).CE);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 8;
  Atom((*p).NZ);

} /* Draw_LYS */

/* ------------------------------------------------------------------------------------ */

void Draw_ARG(ARGp p)

{
  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD);
  theAtom = 8;
  Atom((*p).CZ);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 7;
  Atom((*p).NE);
  theAtom = 9;
  Atom((*p).NH1);
  theAtom = 10;
  Atom((*p).NH2);

} /* Draw_ARG */

/* ------------------------------------------------------------------------------------ */

void Draw_SER(SERp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 5;
  Atom((*p).OG);

} /* Draw_SER */

/* ------------------------------------------------------------------------------------ */

void Draw_THR(THRp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 6;
  Atom((*p).CG2);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 5;
  Atom((*p).OG1);

} /* Draw_THR */

/* ------------------------------------------------------------------------------------ */

void Draw_LEU(LEUp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD1);
  theAtom = 7;
  Atom((*p).CD2);

} /* Draw_LEU */

/* ------------------------------------------------------------------------------------ */

void Draw_ILE(ILEp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG1);
  theAtom = 6;
  Atom((*p).CG2);
  theAtom = 7;
  Atom((*p).CD1);

} /* Draw_ILE */

/* ------------------------------------------------------------------------------------ */

void Draw_VAL(VALp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG1);
  theAtom = 6;
  Atom((*p).CG2);

} /* Draw_VAL */

/* ------------------------------------------------------------------------------------ */

void Draw_PRO(PROp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  //Atom((*p).N); /* FIXME - CI - not sure what to do here... */

} /* Draw_PRO */

/* ------------------------------------------------------------------------------------ */

void Draw_CYS(CYSp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  AtomSymbol = 16; /*AtomKind == 'S'*/
  theAtom = 5;
  Atom((*p).SG);

} /* Draw_CYS */

/* ------------------------------------------------------------------------------------ */

void Draw_MET(METp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 7;
  Atom((*p).CE);
  AtomSymbol = 16; /*AtomKind == 'S'*/
  theAtom = 6;
  Atom((*p).SD);

} /* Draw_MET */

/* ------------------------------------------------------------------------------------ */

void Draw_ASP(ASPp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 6;
  Atom((*p).OD1);
  theAtom = 7;
  Atom((*p).OD2);

} /* Draw_ASP */

/* ------------------------------------------------------------------------------------ */

void Draw_ASN(ASNp p)

{

  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 6;
  Atom((*p).OD1);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 7;
  Atom((*p).ND2);

} /* Draw_ASN */

/* ------------------------------------------------------------------------------------ */

void Draw_GLU(GLUp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 7;
  Atom((*p).OE1);
  theAtom = 8;
  Atom((*p).OE2);

} /* Draw_GLU */

/* ------------------------------------------------------------------------------------ */

void Draw_GLN(GLNp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 7;
  Atom((*p).OE1);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 8;
  Atom((*p).NE2);

} /* Draw_GLN */

/* ------------------------------------------------------------------------------------ */

void Draw_PHE(PHEp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD1);
  theAtom = 7;
  Atom((*p).CD2);
  theAtom = 8;
  Atom((*p).CE1);
  theAtom = 9;
  Atom((*p).CE2);
  theAtom = 10;
  Atom((*p).CZ);

} /* Draw_PHE */

/* ------------------------------------------------------------------------------------ */

void Draw_HIS(HISp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 7;
  Atom((*p).CD2);
  theAtom = 8;
  Atom((*p).CE1);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 9;
  Atom((*p).NE2);
  theAtom = 6;
  Atom((*p).ND1);

} /* Draw_HIS */

/* ------------------------------------------------------------------------------------ */

void Draw_TRP(TRPp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD1);
  theAtom = 7;
  Atom((*p).CD2);
  theAtom = 9;
  Atom((*p).CE2);
  theAtom = 10;
  Atom((*p).CE3);
  theAtom = 11;
  Atom((*p).CZ2);
  theAtom = 12;
  Atom((*p).CZ3);
  theAtom = 13;
  Atom((*p).CH2);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 8;
  Atom((*p).NE1);

} /* Draw_TRP */

/* ------------------------------------------------------------------------------------ */

void Draw_TYR(TYRp p)

{
  if (Draw_BackBone((GLYp)p))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 4;
  Atom((*p).CB);
  theAtom = 5;
  Atom((*p).CG);
  theAtom = 6;
  Atom((*p).CD1);
  theAtom = 7;
  Atom((*p).CD2);
  theAtom = 8;
  Atom((*p).CE1);
  theAtom = 9;
  Atom((*p).CE2);
  theAtom = 10;
  Atom((*p).CZ);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 11;
  Atom((*p).OH);

} /* Draw_TYR */

/* ------------------------------------------------------------------------------------ */

void Draw_GLY(GLYp p) { Draw_BackBone(p); } /* Draw_GLY */

/* ------------------------------------------------------------------------------------ */

short Draw_NT_BackBone(nt_Ap p, unsigned short rnaOpos)

{
  /* ------ trace "backbone" -------- */

  /* ring */
  AtomSymbol = 6; /*AtomKind = 'C';*/
  theAtom = 10;
  Atom((*p).C1_);
  theAtom = 9;
  Atom((*p).C2_);
  theAtom = 7;
  Atom((*p).C3_);
  theAtom = 5;
  Atom((*p).C4_);
  theAtom = 4;
  Atom((*p).C5_);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  if (rnaOpos)
  {
    theAtom = rnaOpos;
    Atom((*(HETATMp)p).coord[rnaOpos]);
  }
  theAtom = 1;
  Atom((*p).O1P);
  theAtom = 2;
  Atom((*p).O2P);
  theAtom = 3;
  Atom((*p).O5_);
  theAtom = 6;
  Atom((*p).O4_);
  theAtom = 8;
  Atom((*p).O3_);
  AtomSymbol = 16; /*AtomKind == 'P'*/
  theAtom = 0;
  Atom((*p).P);

  return (0);

} /* Draw_NT_BackBone */

/* ------------------------------------------------------------------------------------ */

void Draw_A(nt_Ap p)

{
  unsigned short rnaOpos;

  if ((*p).header.nbAtoms != kNTA_nbatms)
    rnaOpos = 0;
  else
    rnaOpos = kNTA_nbatms - 1; /* pos of extra O present in RNA */

  if (Draw_NT_BackBone((nt_Ap)p, rnaOpos))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  //theAtom = 10; // already in backbone
  //Atom((*p).C1_);
  theAtom = 12;
  Atom((*p).C8);
  theAtom = 14;
  Atom((*p).C5);
  theAtom = 15;
  Atom((*p).C6);
  theAtom = 18;
  Atom((*p).C2);
  theAtom = 20;
  Atom((*p).C4);
  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 11;
  Atom((*p).N9);
  theAtom = 13;
  Atom((*p).N7);
  theAtom = 16;
  Atom((*p).N6);
  theAtom = 17;
  Atom((*p).N1);
  theAtom = 19;
  Atom((*p).N3);

} /* Draw_A */

/* ------------------------------------------------------------------------------------ */

void Draw_TorU(nt_Tp p, Boolean isT)

{
  unsigned short rnaOpos;

  if (isT)
    rnaOpos = 0;
  else
    rnaOpos = kNTU_nbatms - 1; /* pos of extra O present in RNA */

  if (Draw_NT_BackBone((nt_Ap)p, rnaOpos))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  //theAtom = 10; // already in backbone
  //Atom((*p).C1_);
  theAtom = 17;
  Atom((*p).C5);
  theAtom = 18;
  if (isT)
  {
    Atom((*p).C5M);
    theAtom = 19;
    Atom((*p).C6);
  }
  else
    Atom((*(nt_U *)p).C6);
  theAtom = 12;
  Atom((*p).C2);
  theAtom = 15;
  Atom((*p).C4);

  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 14;
  Atom((*p).N3);
  theAtom = 11;
  Atom((*p).N1);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 13;
  Atom((*p).O2);
  theAtom = 16;
  Atom((*p).O4);

} /* Draw_TorU */

/* ------------------------------------------------------------------------------------ */

void Draw_C(nt_Cp p)

{
  unsigned short rnaOpos;

  if ((*p).header.nbAtoms != kNTC_nbatms)
    rnaOpos = 0;
  else
    rnaOpos = kNTC_nbatms - 1; /* pos of extra O present in RNA */

  if (Draw_NT_BackBone((nt_Ap)p, rnaOpos))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  //theAtom = 10; // already in backbone
  //Atom((*p).C1_);
  theAtom = 12;
  Atom((*p).C2);
  theAtom = 16;
  Atom((*p).C4);
  theAtom = 17;
  Atom((*p).C5);
  theAtom = 18;
  Atom((*p).C6);

  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 11;
  Atom((*p).N1);
  theAtom = 14;
  Atom((*p).N3);
  theAtom = 15;
  Atom((*p).N4);
  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 13;
  Atom((*p).O2);

} /* Draw_C */

/* ------------------------------------------------------------------------------------ */

void Draw_G(nt_Gp p)

{
  unsigned short rnaOpos;

  if ((*p).header.nbAtoms != kNTG_nbatms)
    rnaOpos = 0;
  else
    rnaOpos = kNTG_nbatms - 1; /* pos of extra O present in RNA */

  if (Draw_NT_BackBone((nt_Ap)p, rnaOpos))
    return;

  AtomSymbol = 6; /*AtomKind = 'C';*/
  //theAtom = 10; // already in backbone
  //Atom((*p).C1_);
  theAtom = 18;
  Atom((*p).C2);
  theAtom = 21;
  Atom((*p).C4);
  theAtom = 14;
  Atom((*p).C5);
  theAtom = 15;
  Atom((*p).C6);
  theAtom = 12;
  Atom((*p).C8);

  AtomSymbol = 7; /*AtomKind = 'N';*/
  theAtom = 17;
  Atom((*p).N1);
  theAtom = 19;
  Atom((*p).N2);
  theAtom = 20;
  Atom((*p).N3);
  theAtom = 13;
  Atom((*p).N7);
  theAtom = 11;
  Atom((*p).N9);

  AtomSymbol = 8; /*AtomKind = 'O';*/
  theAtom = 16;
  Atom((*p).O6);

} /* Draw_G */

/* ------------------------------------------------------------------------------------ */
