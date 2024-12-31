//
//  init.h
//  PDBV
//
//  Created by Christian Iseli on 22.09.2024.
//
#pragma once

#include "PMPDBVDLL.h"

extern "C" {

#ifndef init_h
#define init_h

#include <stdio.h>

/* ----PROTOTYPES---------------------------------------------------------------------- */

void ShowMessage(const char *s, Boolean abortSPDBV);
void ShowMessageCode(const char *s, long e, Boolean abortSPDBV);
//short InitGlobals(void);

#endif /* init_h */
}
