//
//  SwiftCallback.c
//  PDBV
//
//  Created by Christian Iseli on 24.09.2024.
//

#include "pch.h"
#include "PMPDBVDLL.h"
#include "Globals.h"
#include "RenderG.h"
#include "render.h"
#include "SwiftCallback.h"

static SwiftCallbacks sCallbacks;
static void GL_resize(void);
static void Render_Molec(short layer);

extern void SwiftCBSetup(const SwiftCallbacks *callbacks) { sCallbacks = *callbacks; }

void GL_atm(FPoint3D A)
{
  GLYp AGptr = GroupsHdl[theGroup];
  sCallbacks.GL_atm(A, theLayer, theGroup, theAtom, AGptr->header.chain, 0);
}

void GL_bnd(FPoint3D A, FPoint3D B)
{
  GLYp AGptr = GroupsHdl[theGroup];
  sCallbacks.GL_bnd(A, B, theLayer, theGroup, theAtom, theAtomB, AGptr->header.chain, 0);
}

void GL_init_Atom_Colors (void) {
  sCallbacks.GL_Atom_Color("O", 1.0, 0.0, 0.0, true);
  sCallbacks.GL_Atom_Color("N", 0.24, 0.32, 1.0, false);
  sCallbacks.GL_Atom_Color("C", 1.0, 1.0, 1.0, false);
  sCallbacks.GL_Atom_Color("S", 1.0, 0.94, 0.0, false);
  sCallbacks.GL_Atom_Color("P", 1.0, 0.66, 0.09, false);
  sCallbacks.GL_Atom_Color("H", 0.3, 0.8, 1.0, false);
  sCallbacks.GL_Atom_Color("X", 0.625, 0.625, 0.625, false);
  sCallbacks.GL_Atom_Color("StrongHBond", 0.0, 1.0, 0.0, false);
  sCallbacks.GL_Atom_Color("WeakHBond", 0.5625, 0.6875, 0.5625, false);
  sCallbacks.GL_Atom_Color("Clash", 1.0, 0.43, 1.0, false);
  sCallbacks.GL_Atom_Color("SSbnd", 1.0, 0.94, 0.0, false);
} /* GL_init_Atom_Colors */

void UpdateMoleculeGLstatus(short layer)
{
  //    myLog("Forcing Rebuilding of glMolecule for layer %d\n",layer);
  if (sCallbacks.GL_ForceRebuild != NULL)
      sCallbacks.GL_ForceRebuild(layer);
  if (sCallbacks.GL_All_Invisible != NULL)
    sCallbacks.GL_All_Invisible(layer);
  GLRender();
    
}
/* ------------------------------------------------------------------------------------ */

void GLRender(void)
{
  GL_resize();
}

static void GL_resize(void)
{
  short layer;

  for (layer = 0; layer <= nbLayer; layer++)
  {
    Render_Molec(layer);
  }
} /* GL_resize */
/* ------------------------------------------------------------------------------------ */

static void Render_Molec(short layer)
{
  register unsigned short i;

  theLayer = layer;

  GroupsHdl = PDB[layer].GroupsHdl;
  AtomNamesHdl = PDB[layer].AtomNamesHdl;

  for (i = 0; i < PDB[layer].nbAtomGroups; i++)
  {

    theGroup = i;
    GLYp AGptr = GroupsHdl[i];

    switch ((*(GLYp)AGptr).header.kind)
    {
    case kLEU:
      Draw_LEU((LEUp)AGptr);
      break;
    case kALA:
      Draw_ALA((ALAp)AGptr);
      break;
    case kGLY:
      Draw_GLY((GLYp)AGptr);
      break;
    case kSER:
      Draw_SER((SERp)AGptr);
      break;
    case kVAL:
      Draw_VAL((VALp)AGptr);
      break;
    case kGLU:
      Draw_GLU((GLUp)AGptr);
      break;
    case kLYS:
      Draw_LYS((LYSp)AGptr);
      break;
    case kTHR:
      Draw_THR((THRp)AGptr);
      break;
    case kASP:
      Draw_ASP((ASPp)AGptr);
      break;
    case kILE:
      Draw_ILE((ILEp)AGptr);
      break;
    case kPRO:
      Draw_PRO((PROp)AGptr);
      break;
    case kARG:
      Draw_ARG((ARGp)AGptr);
      break;
    case kASN:
      Draw_ASN((ASNp)AGptr);
      break;
    case kGLN:
      Draw_GLN((GLNp)AGptr);
      break;
    case kPHE:
      Draw_PHE((PHEp)AGptr);
      break;
    case kTYR:
      Draw_TYR((TYRp)AGptr);
      break;
    case kHIS:
      Draw_HIS((HISp)AGptr);
      break;
    case kMET:
      Draw_MET((METp)AGptr);
      break;
    case kCYS:
      Draw_CYS((CYSp)AGptr);
      break;
    case kTRP:
      Draw_TRP((TRPp)AGptr);
      break;
    case kG:
      Draw_G((nt_Gp)AGptr);
      break;
    case kA:
      Draw_A((nt_Ap)AGptr);
      break;
    case kT:
      Draw_TorU((nt_Tp)AGptr, true);
      break;
    case kC:
      Draw_C((nt_Cp)AGptr);
      break;
    case kU:
      Draw_TorU((nt_Tp)AGptr, false);
      break;
    default:
      break;
    } /*switch*/
  } /* for every acide amine */
} /* Render_Molec */
