//
//  utils.c
//  PDBV
//
//  Created by Christian Iseli on 22.09.2024.
//

#include "pch.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Globals.h"
#include "utils.h"

/* ------------------------------------------------------------------------------------ */

void *safe_realloc(void **x, int l)
{
  void *tmp = *x;

  /*myLog("REALLOC_IN %d\n", *x); */
  tmp = realloc(tmp, l);
  if (tmp != NULL)
  {
    *x = tmp;
  }
  /*myLog("REALLOC_out %d\n", *x); */
  return (tmp);

} /* safe_realloc */

/* ------------------------------------------------------------------------------------ */

/* this routine should be used for column formatted files, e.g. PDB;
    it initializes nChar characters in linbuf and reads max. this
    amount into the buffer .... */

short myfgets_(char *linbuf, int nChar, FILE *infile)
{
  unsigned short i;
  char c;

  /* Initialize the buffer with \0 */
  memset(linbuf, '\0', nChar);

  i = 0;
  do
  {
    if (i == kBufSize)
      return (kEOF);
    linbuf[i++] = c = fgetc(infile);
    if (feof(infile))
      return (kEOF);
  } while ((c != '\n') && (c != '\r') && (i < nChar));

  linbuf[i - 1] = 0;

  c = fgetc(infile);

  if ((c != '\n') && (c != '\r'))
    ungetc(c, infile);

  return (0);

} /*  myfgets_ (char *linbuf, int nChar, FILE *infile) */

const char *bundleRsrcDir = "/tmp";
const char *downloadDir = "/tmp/download";
const char *tempDir = "/tmp/temp";

/* ------------------------------------------------------------------------------------ */
