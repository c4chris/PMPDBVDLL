//
//  init.c
//  PDBV
//
//  Created by Christian Iseli on 22.09.2024.
//

#include "pch.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "Globals.h"
#include "init.h"
#include "utils.h"
#include "pdb.h"
#include "SwiftCallback.h"

int myLog(const char* fmt, ...)
{
  va_list va;
  va_start(va, fmt);
  int res = -1;
  if (printfDelegate != NULL) {
    res = printfDelegate(fmt, va);
  } else {
    res = vfprintf(stdout, fmt, va);
  }
  return res;
}

static void ResetLayerInfos(short layer)
{

  PDB[layer].GroupsHdl = NULL;
  PDB[layer].AtomNamesHdl = NULL;
  PDB[layer].flags = 0;
  PDB[layer].nbAtomGroups = 0;
  PDB[layer].nbChains = 0;
  UpdateMoleculeGLstatus(layer);

} /* ResetLayerInfos */

short InitGlobals(void)
{

  short layer;

  myLog("This is InitGlobals() starting\n");
  myLog("bundleRsrcDir is %s\n", bundleRsrcDir);
  myLog("downloadDir is %s\n", downloadDir);
  myLog("tempDir is %s\n", tempDir);
  
  nbLayer = -1;

  for (layer = 0; layer < maxLayerNb; layer++)
  {
    ResetLayerInfos(layer);
  }

  myLog("This is InitGlobals() ending\n");
  return (1);

} /* InitGlobals */

/* ------------------------------------------------------------------------------------ */
void ShowMessage(const char *s, Boolean abortSPDBV)
{
  myLog("Message: %s\n", s);
  if (abortSPDBV)
    exit(999);
} /* ShowMessage */

/* ------------------------------------------------------------------------------------ */

void ShowMessageCode(const char *s, long e, Boolean abortSPDBV)
{
  myLog("MessageCode: %s (%ld)\n", s, e);
  if (abortSPDBV)
    exit(999);
} /* ShowMessageCode */

/* ------------------------------------------------------------------------------------ */
