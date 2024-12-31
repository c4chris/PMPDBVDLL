//
//  render.h
//  PDBV
//
//  Created by Christian Iseli on 23.09.2024.
//
#pragma once

#include "PMPDBVDLL.h"

extern "C" {

#ifndef render_h
#define render_h

#include <stdio.h>

short Draw_BackBone(GLYp p);
void Draw_ALA(ALAp p);
void Draw_GLY(GLYp p);
void Draw_ALA(ALAp p);
void Draw_LYS(LYSp p);
void Draw_ARG(ARGp p);
void Draw_SER(SERp p);
void Draw_THR(THRp p);
void Draw_LEU(LEUp p);
void Draw_ILE(ILEp p);
void Draw_VAL(VALp p);
void Draw_PRO(PROp p);
void Draw_CYS(CYSp p);
void Draw_MET(METp p);
void Draw_ASP(ASPp p);
void Draw_ASN(ASNp p);
void Draw_GLU(GLUp p);
void Draw_GLN(GLNp p);
void Draw_PHE(PHEp p);
void Draw_HIS(HISp p);
void Draw_TRP(TRPp p);
void Draw_TYR(TYRp p);
void Draw_GLY(GLYp p);

short Draw_NT_BackBone(nt_Ap p, unsigned short rnaOpos);
void Draw_A(nt_Ap p);
void Draw_TorU(nt_Tp p, Boolean isT);
void Draw_C(nt_Cp p);
void Draw_G(nt_Gp p);
void Draw_U(nt_Up p);

#endif /* render_h */
}
