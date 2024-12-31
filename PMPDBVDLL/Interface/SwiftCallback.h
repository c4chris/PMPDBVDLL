//
//  SwiftCallback.h
//  PDBV
//
//  Created by Christian Iseli on 24.09.2024.
//
#pragma once

#include "PMPDBVDLL.h"

extern "C" {

#ifndef SwiftCallback_h
#define SwiftCallback_h

#include <stdio.h>
#ifdef _MSC_VER
#define _Nonnull
#endif

//struct SwiftCallbacks
//{
//  void (*_Nonnull GL_atm)(FPoint3D A, short layer, unsigned short group, short atom, char chain, long color);
//  void (*_Nonnull GL_bnd)(FPoint3D A, FPoint3D B, short layer, unsigned short group, short atomA, short atomB, char chain, long color);
//  void (*_Nonnull GL_Atom_Color)(const char * _Nonnull name, double red, double green, double blue, Boolean isMetallic);
//  void (*_Nonnull GL_All_Invisible)(short layer);
//  void (*_Nonnull GL_ForceRebuild)(short layer);
//  void (*_Nonnull GL_New_Shape)(void);
//  void (*_Nonnull GL_Add_Vertex)(Point3D V, Point3D N, unsigned int color);
//  void (*_Nonnull GL_Draw_Shape)(short layer, char type);
//};
//typedef struct SwiftCallbacks SwiftCallbacks;

//void SwiftCBSetup(const SwiftCallbacks *_Nonnull callbacks);
void UpdateMoleculeGLstatus(short layer);
void GLRender(void);

#endif /* SwiftCallback_h */
}
