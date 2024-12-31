//
//  utils.h
//  PDBV
//
//  Created by Christian Iseli on 22.09.2024.
//
#pragma once

#include "PMPDBVDLL.h"

extern "C" {

#ifndef utils_h
#define utils_h

#include <stdio.h>

//extern const char* bundleRsrcDir;
//extern const char* downloadDir;
//extern const char* tempDir;

/* ----PROTOTYPES---------------------------------------------------------------------- */

void *safe_realloc(void **x, int l);

#define kEOF 1
#define kBufSize 1024
#define kOccupancyDotPos 57
#define kBfactorDotPos 63

short myfgets_(char *linbuf, int nChar, FILE *infile);

#endif /* utils_h */
}
