//
//  pdb.c
//  PDBV
//
//  Created by Christian Iseli on 22.09.2024.
//

#include "pch.h"
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Globals.h"

#include "pdb.h"
#include "init.h"
#include "utils.h"

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

/* ------------------------------------------------------------------------------------ */

#define kOverflow -1
#define kMaxAtomBuffer 75
#define kReadGroupChunk 250L /* allocate (& increase) memory for groups by 100 atoms groups chunks steps */

static short layer;

/* info about group currently being parsed */
static char currentName[8];
static char currentLabel[8];
static char currentChain;
static char currentAltResNum, justReadAltResNum;
static char current_label_alt_id; /* __ts */
static char accepted_label_alt_id;
static unsigned long currentAtomNumberStart;
static GENERIC_aa generic_aa[kMaxAtomBuffer]; /* structure used to load all atoms of the group under process, before dispatching into the proper aa type */

static char linbuf[kBufSize + 2];
static FILE *infile;
static GLYp *GroupsHdl;
static ATOM_TYPEp *AtomNamesHdl;

static unsigned long attributes;

static Point3D R[4]; /* matrix */
static char justReadChain;
static unsigned long isWaterFlag;
static long posWhenCloseFile;
static Boolean hasSkippedANISOU;
static int uniqueID = 1;

static short processFile(void);
static short GetAtom(Point3D *coordPtr, ATOM_TYPEp typePtr, const char query[5], const char atm_type[2], short n);
static Boolean GetNextAtom(PDB_line *pdb);
static short Read_aa(unsigned short i);
static short Read_NT_Backbone(unsigned short i, long AGsize, long ATsize, GLYp *p, ATOM_TYPEp *AtTypePtr);
static short Read_A(nt_Ap p, ATOM_TYPEp *AtTypePtr, char dnaflag);
static short Read_G(nt_Gp p, ATOM_TYPEp *AtTypePtr, char dnaflag);
static short Read_T(nt_Tp p, ATOM_TYPEp *AtTypePtr, char dnaflag);
static short Read_C(nt_Cp p, ATOM_TYPEp *AtTypePtr, char dnaflag);
static short Read_U(nt_Up p, ATOM_TYPEp *AtTypePtr, char dnaflag);

/* ------------------------------------------------------------------------------------ */

int AllocateMemoryForNewLayer(short layerL)
{
  /* allocate memory for Groups and Atom Names */
  PDB[layerL].GroupsHdl = (GLYp *) malloc((long)kReadGroupChunk * sizeof(void *));
  if (PDB[layerL].GroupsHdl == NULL)
  {
    ShowMessage("Not enough memory for a new Layer", FALSE);
    return (1);
  }
  PDB[layerL].AtomNamesHdl = (ATOM_TYPEp *) malloc((long)kReadGroupChunk * sizeof(void *));
  if (PDB[layerL].AtomNamesHdl == NULL)
  {
    ShowMessage("Not enough memory for a new Layer", FALSE);
    free(PDB[layerL].GroupsHdl);
    return (1);
  }
  return (0);

} /* AllocateMemoryForNewLayer */
/* ------------------------------------------------------------------------------------ */

void doPDBinput(const char *inName, const char *layerName)
{
  long pos = 0;
  nbLayer += 1;
  //strcpy(PDB[nbLayer].originalFilename, inName);
  strcpy_s(PDB[nbLayer].originalFilename, kmaxPathChar, inName);
  //strcpy(PDB[nbLayer].LayerName, layerName);
  strcpy_s(PDB[nbLayer].LayerName, kmaxNameChar, layerName);
  if (pdb_Input(nbLayer, NULL, inName, 0, &pos) != 0)
  {
    myLog("Failed to loading PDB %s\n", inName);
    return;
  }
  myLog("Done loading PDB %s, layer = %u, nbGroups = %u, LayerName = %s\n", inName, nbLayer, PDB[nbLayer].nbAtomGroups, PDB[nbLayer].LayerName);

}

short pdb_Input(short whichlayer, FILE *outfile, const char *inName, long filesize, long *openFilePos)
{
  char s[256];
  short err;

  layer = whichlayer; /* fill static variable with proper value */

  /* open file */

  //if (!(infile = fopen(inName, "rb")))
  if (fopen_s(&infile, inName, "rb") != 0)
  {
    //strcpy(s, "Cannot open file: ");
    //strcat(s, inName);
    strcpy_s(s, 256, "Cannot open file: ");
    strcat_s(s, 256, inName);
    ShowMessage(s, FALSE);

    *openFilePos = 0L;
    return (1);
  }

  if (AllocateMemoryForNewLayer(layer))
    return (1);
  GroupsHdl = PDB[layer].GroupsHdl;
  AtomNamesHdl = PDB[layer].AtomNamesHdl;

  err = processFile();
  fclose(infile);
  if (err == kOverflow)
  {
    //strcpy(s, "Sorry, I cannot load ");
    //strcat(s, PDB[layer].LayerName);
    //strcat(s, ": out of memory");
    strcpy_s(s, 256, "Sorry, I cannot load ");
    strcat_s(s, 256, PDB[layer].LayerName);
    strcat_s(s, 256, ": out of memory");
    ShowMessage(s, FALSE);
    return (1);
  }

  if (PDB[layer].nbAtomGroups == 0)
  {
    //strcpy(s, PDB[layer].LayerName);
    //strcat(s, ": File Ignored (either it is not a valid PDB file, or it contains only a Carbon Alpha trace). Would you like to Open the file as text ?");
    strcpy_s(s, 256, PDB[layer].LayerName);
    strcat_s(s, 256, ": File Ignored (either it is not a valid PDB file, or it contains only a Carbon Alpha trace). Would you like to Open the file as text ?");
    return (1);
  }

  return (0);

} /*pdb_Input */

/* ------------------------------------------------------------------------------------ */

/*
ATOM      1  O5*   C A   1       6.236  17.562  -9.026  1.00 27.53      1D56  59
*/
static Boolean GetNextAtom(PDB_line *pdb)
{
  float Bfactor, occupancy;

  do
  {
    if (myfgets_(linbuf, kBufSize, infile) == kEOF)
      return (true);
    if (strncmp(linbuf, "//", 2) == 0)
      continue;
    if (strncmp(linbuf, "ANISOU", 6) == 0)
    {
      hasSkippedANISOU = true;
      continue;
    }
    if ((strncmp(linbuf, "ATOM", 4) != 0) && (strncmp(linbuf, "HET", 3) != 0))
      return (false);
    else
    {
      //sscanf(linbuf, "%s", &(*pdb).type[0]);
      //sscanf(&linbuf[7], "%d", &(*pdb).atomNumber);
      //sscanf(&linbuf[12], "%4c", &(*pdb).atomName[0]);
      //sscanf(&linbuf[17], "%3c", &(*pdb).labelName[0]);
#if (kMaxCharForResNum == 3)
      //sscanf(&linbuf[23], "%3c", &(*pdb).labelNumber[0]); /* in fact should be 22, and read 4 chars */
#else
      //sscanf(&linbuf[22], "%4c", &(*pdb).labelNumber[0]); /* in fact should be 22, and read 4 chars */
#endif
      //sscanf(&linbuf[30], "%8lf%8lf%8lf", &(*pdb).coord.x, &(*pdb).coord.y, &(*pdb).coord.z);
      sscanf_s(linbuf, "%s", &(*pdb).type[0], 8);
      sscanf_s(&linbuf[7], "%d", &(*pdb).atomNumber);
      sscanf_s(&linbuf[12], "%4c", &(*pdb).atomName[0], 8);
      sscanf_s(&linbuf[17], "%3c", &(*pdb).labelName[0], 8);
#if (kMaxCha_srForResNum == 3)
      sscanf_s(&linbuf[23], "%3c", &(*pdb).labelNumber[0]); /* in fact should be 22, and read 4 chars */
#else
      sscanf_s(&linbuf[22], "%4c", &(*pdb).labelNumber[0], 8); /* in fact should be 22, and read 4 chars */
#endif
      sscanf_s(&linbuf[30], "%8lf%8lf%8lf", &(*pdb).coord.x, &(*pdb).coord.y, &(*pdb).coord.z);
      if (pdb->coord.x > 9999.98)
        pdb->coord.x = kNoCoord;
      if (pdb->coord.y > 9999.98)
        pdb->coord.y = kNoCoord;
      if (pdb->coord.z > 9999.98)
        pdb->coord.z = kNoCoord;
#ifdef INSERTION_CODE
      (*pdb).res_insertion = linbuf[26]; /* __ts */
#endif
      current_label_alt_id = linbuf[16]; /* __ts then ng */
      if ((current_label_alt_id != ' ') && (accepted_label_alt_id == ' '))
        accepted_label_alt_id = current_label_alt_id; /*  remember first alternate id  we see and ignore all other later */

      if (linbuf[kOccupancyDotPos] == '.')
      {
        char s[8];
        //strncpy(s, &linbuf[kOccupancyDotPos - 2], 5);
        strncpy_s(s, 8, &linbuf[kOccupancyDotPos - 2], 5);
        s[5] = 0;
        //sscanf(s, "%f", &occupancy);
        sscanf_s(s, "%f", &occupancy);
        (*pdb).atm_occupancy = (unsigned short)(100 * (occupancy + 0.001));
      }
      else
        (*pdb).atm_occupancy = 100;

      if (linbuf[kBfactorDotPos] == '.')
      { /* There is a problem with scanf, if there is no space between the columns
            therefor I copy the columns first to the stack .... __TS */
        char Bfac_column[12];
        //strncpy(Bfac_column, &linbuf[kBfactorDotPos - 3], 6);
        //sscanf(Bfac_column, "%f", &Bfactor);
        strncpy_s(Bfac_column, 12, &linbuf[kBfactorDotPos - 3], 6);
        sscanf_s(Bfac_column, "%f", &Bfactor);
        if (Bfactor <= 0)
        {
          Bfactor = 0;
        }

        (*pdb).Bfactor = (unsigned short)(100 * (Bfactor + 0.001));
      }
      else
        (*pdb).Bfactor = 9999;

      justReadChain = linbuf[21]; /* awful, must clean this */
      justReadAltResNum = linbuf[26];
      return (false);
    }
  } while (true);

} /*GetNextAtom */
/* ------------------------------------------------------------------------------------ */

static short processFile(void)
{
  PDB_line pdb;
  unsigned short i;
  short k;
  GLYp AGptr;
  ATOM_TYPEp AtomTypePtr;
  char s[100];
  short err;

  accepted_label_alt_id = ' ';
  do
  {
    if (myfgets_(linbuf, kBufSize, infile) == kEOF)
      return (true);
    pdb.type[0] = 0;
    //sscanf(linbuf, "%s", &pdb.type[0]);      /* quick scan to retrieve type of infos stored on this file line  */
    sscanf_s(linbuf, "%s", &pdb.type[0], 8);      /* quick scan to retrieve type of infos stored on this file line  */
  } while (strncmp(pdb.type, "//", 2) == 0); /* skip groups that have been commented */

  /* --------------- Let's input the file ---------------- */
  do
  {

    i = PDB[layer].nbAtomGroups;

    /* Increase group Handle if necessary */
    if ((i % kReadGroupChunk) == (kReadGroupChunk - 1))
    {
      if ((safe_realloc((void **)&PDB[layer].GroupsHdl, (i + kReadGroupChunk) * sizeof(void *))) == 0)
        return (kOverflow);
      if ((safe_realloc((void **)&PDB[layer].AtomNamesHdl, (i + kReadGroupChunk) * sizeof(void *))) == 0)
        return (kOverflow);
      GroupsHdl = PDB[layer].GroupsHdl;
      AtomNamesHdl = PDB[layer].AtomNamesHdl;
    }

    pdb.type[0] = 0;
    //sscanf(linbuf, "%s", &pdb.type[0]); /* quick scan to retrieve type of infos stored on this file line  */
    sscanf_s(linbuf, "%s", &pdb.type[0], 8); /* quick scan to retrieve type of infos stored on this file line  */

    /* treat ATOM type line info */
    if (strncmp(pdb.type, "ATOM", 4) == 0)
    {
      int index;
      float Bfactor, occupancy;

      pdb.labelName[0] = 0;
      //sscanf(linbuf, "%s", &pdb.type[0]);
      //sscanf(&linbuf[6], "%d", &pdb.atomNumber);
      //sscanf(&linbuf[12], "%4c", &pdb.atomName[0]);
      sscanf_s(linbuf, "%s", &pdb.type[0], 8);
      sscanf_s(&linbuf[6], "%d", &pdb.atomNumber);
      sscanf_s(&linbuf[12], "%4c", &pdb.atomName[0], 8);
      pdb.atomName[4] = 0;
      /*            sscanf(linbuf,"%s %d%1c%4c",&pdb.type,&pdb.atomNumber,&trash,&pdb.atomName);*/
      //sscanf(&linbuf[17], "%3c", &pdb.labelName[0]);
      sscanf_s(&linbuf[17], "%3c", &pdb.labelName[0], 8);
      pdb.labelName[3] = 0;
#if (kMaxCharForResNum == 3)
      sscanf(&linbuf[23], "%3c", &pdb.labelNumber[0]); /* in fact should be 22, and read 4 chars */
#else
      //sscanf(&linbuf[22], "%4c", &pdb.labelNumber[0]); /* in fact should be 22, and read 4 chars */
      sscanf_s(&linbuf[22], "%4c", &pdb.labelNumber[0], 8); /* in fact should be 22, and read 4 chars */
#endif
      pdb.labelNumber[kMaxCharForResNum] = 0;
      //sscanf(&linbuf[30], "%8lf%8lf%8lf", &pdb.coord.x, &pdb.coord.y, &pdb.coord.z);
      sscanf_s(&linbuf[30], "%8lf%8lf%8lf", &pdb.coord.x, &pdb.coord.y, &pdb.coord.z);

      if (pdb.coord.x > 9999.98)
        pdb.coord.x = kNoCoord;
      if (pdb.coord.y > 9999.98)
        pdb.coord.y = kNoCoord;
      if (pdb.coord.z > 9999.98)
        pdb.coord.z = kNoCoord;
#ifdef INSERTION_CODE
      pdb.res_insertion = linbuf[26]; /* __ts */
#endif
      current_label_alt_id = linbuf[16];
      if ((current_label_alt_id != ' ') && (accepted_label_alt_id == ' '))
        accepted_label_alt_id = current_label_alt_id; /*  remember first alternate id  we see and ignore all other later */

      if (linbuf[kOccupancyDotPos] == '.')
      {
        char s2[8];
        //strncpy(s2, &linbuf[kOccupancyDotPos - 2], 5);
        strncpy_s(s2, 8, &linbuf[kOccupancyDotPos - 2], 5);
        s2[5] = 0;
        //sscanf(s2, "%f", &occupancy);
        sscanf_s(s2, "%f", &occupancy);
        pdb.atm_occupancy = (unsigned short)(100 * (occupancy + 0.001));
      }
      else
        pdb.atm_occupancy = 100;
      if (linbuf[kBfactorDotPos] == '.')
      { /* There is a problem with scanf, if there is no space between the columns
            therefor I copy the columns first to the stack .... __TS */
        char Bfac_column[12];
        //strncpy(Bfac_column, &linbuf[kBfactorDotPos - 3], 6);
        //sscanf(Bfac_column, "%f", &Bfactor);
        strncpy_s(Bfac_column, 12, &linbuf[kBfactorDotPos - 3], 6);
        sscanf_s(Bfac_column, "%f", &Bfactor);
        if (Bfactor <= 0)
        {
          Bfactor = 0;
        }

        pdb.Bfactor = (unsigned short)(100 * (Bfactor + 0.001));
      }
      else
      {
        pdb.Bfactor = 9999;
      }
      /* check if aa or NT */
      index = -1;
      do
      {
        if (index == 20)
          break;
      } while (strncmp(aaName[++index].aa, pdb.labelName, 3) != 0);
      if ((index == 20) && (strncmp(pdb.labelName, " DT", 3) != 0) && (strncmp(pdb.labelName, " DC", 3) != 0) && (strncmp(pdb.labelName, " DG", 3) != 0) &&
          (strncmp(pdb.labelName, " DA", 3) != 0) && (strncmp(pdb.labelName, "  T", 3) != 0) && (strncmp(pdb.labelName, "  C", 3) != 0) &&
          (strncmp(pdb.labelName, "  G", 3) != 0) && (strncmp(pdb.labelName, "  A", 3) != 0) && (strncmp(pdb.labelName, "  U", 3) != 0))
      {
	/* assume that it is a modified amino-acid, and replace it by an ALA */
	ShowMessageCode("Unknown Amino-acid replaced by GLY at pos.", PDB[layer].nbAtomGroups, false);
	i = 65000;
      }

      /* remember name of current ATM group in order to stop as soon as it changes */
      //strncpy(currentName, pdb.labelName, 3);
      strncpy_s(currentName, 8, pdb.labelName, 3);
      currentName[3] = 0;
      //strncpy(currentLabel, pdb.labelNumber, kMaxCharForResNum);
      strncpy_s(currentLabel, 8, pdb.labelNumber, kMaxCharForResNum);
      currentLabel[kMaxCharForResNum] = 0;
      currentChain = linbuf[21];
      currentAltResNum = linbuf[26];

      current_label_alt_id = linbuf[16];
      if ((current_label_alt_id != ' ') && (accepted_label_alt_id == ' '))
        accepted_label_alt_id = current_label_alt_id; /*  remember first alternate id  we see and ignore all other later */
      currentAtomNumberStart = (unsigned long)pdb.atomNumber;

      /* read all infos into the "generic" aa structure until a new group is reached */
      for (k = 0; k < kMaxAtomBuffer; k++)
      {
        generic_aa[k].atm_name[0] = 0;
        generic_aa[k].atm_name[1] = 0;
        generic_aa[k].atm_name[2] = 0;
        generic_aa[k].atm_name[3] = 0;
      }
      k = 0;
      do
      {
        if ((current_label_alt_id == ' ') /*|| (pdb.res_insertion == accepted_label_alt_id)*/ || (current_label_alt_id == accepted_label_alt_id))
        {
          generic_aa[k].chain = currentChain;
          generic_aa[k].altResNum = currentAltResNum;

#ifdef INSERTION_CODE
          generic_aa[k].res_insertion = pdb.res_insertion;
#endif
          generic_aa[k].Bfactor = pdb.Bfactor;
          generic_aa[k].atm_occupancy = pdb.atm_occupancy;
          generic_aa[k].atm_coord = pdb.coord;
          //strncpy(generic_aa[k].atm_name, pdb.atomName, 4); /* store name of atom  (for example CG1) */
          strncpy_s(generic_aa[k].atm_name, 12, pdb.atomName, 4); /* store name of atom  (for example CG1) */

          generic_aa[k++].atm_name[4] = 0;

          pdb.labelName[0] = 0; /* reset label name (to avoid remanence problems) */
        }
        else
        {
          current_label_alt_id = ' '; /* avoid remanence */
        }
        /* now read next line */
        if (GetNextAtom(&pdb)) /* if eof, break */
          break;

      } while ((currentChain == justReadChain) && (strncmp(currentName, pdb.labelName, 3) == 0) &&
               (strncmp(currentLabel, pdb.labelNumber, kMaxCharForResNum) == 0) && (currentAltResNum == justReadAltResNum));

      if (i == 65000) /* (cf above) assume it is a GLY; (often the case with UNK as in pdb file 1155C) */
      {
        if (k < 3)
        {
          myLog("***  EXPDB file %s needs some clean up ***\n", PDB[layer].LayerName);
          exit(1);
        }
        i = PDB[layer].nbAtomGroups;
        //strcpy(pdb.labelName, "GLY");
        //strcpy(currentName, "GLY");
        strcpy_s(pdb.labelName, 8, "GLY");
        strcpy_s(currentName, 8, "GLY");
      }
      if (k > 0) /* else, we probably skipped completely an alternate residue as no atom has ever  been loaded */
        err = Read_aa(i);
      else
        err = 0;
      if (err == kOverflow)
        return (kOverflow);
      else if (err == 99) /* not an amino acid */
      {

        /* check for NT */
        if (currentName[2] == 'T')
        {
          if (Read_NT_Backbone(i, (long)sizeof(nt_T), (long)kNTT_nbatms * sizeof(ATOM_TYPE), &AGptr, &AtomTypePtr) == kOverflow)
            return (kOverflow);
          Read_T((nt_Tp)AGptr, &AtomTypePtr, currentName[1]);
        }
        else if (currentName[2] == 'C')
        {
          if (Read_NT_Backbone(i, (long)sizeof(nt_C), (long)kNTC_nbatms * sizeof(ATOM_TYPE), &AGptr, &AtomTypePtr) == kOverflow)
            return (kOverflow);
          Read_C((nt_Cp)AGptr, &AtomTypePtr, currentName[1]);
        }
        else if (currentName[2] == 'G')
        {
          if (Read_NT_Backbone(i, (long)sizeof(nt_G), (long)kNTG_nbatms * sizeof(ATOM_TYPE), &AGptr, &AtomTypePtr) == kOverflow)
            return (kOverflow);
          Read_G((nt_Gp)AGptr, &AtomTypePtr, currentName[1]);
        }
        else if (currentName[2] == 'A')
        {
          if (Read_NT_Backbone(i, (long)sizeof(nt_A), (long)kNTA_nbatms * sizeof(ATOM_TYPE), &AGptr, &AtomTypePtr) == kOverflow)
            return (kOverflow);
          Read_A((nt_Ap)AGptr, &AtomTypePtr, currentName[1]);
        }
        else if (currentName[2] == 'U')
        {
          if (Read_NT_Backbone(i, (long)sizeof(nt_U), (long)kNTU_nbatms * sizeof(ATOM_TYPE), &AGptr, &AtomTypePtr) == kOverflow)
            return (kOverflow);
          Read_U((nt_Up)AGptr, &AtomTypePtr, currentName[1]);
        }
        else
        {
          //strcpy(s, "Unknown group: ");
          //strncat(s, currentName, 3);
          strcpy_s(s, 100, "Unknown group: ");
          strncat_s(s, 100, currentName, 3);
          ShowMessageCode(s, pdb.atomNumber, FALSE); /* unknown struct */
        }
      }
    } /* "ATOM" */
    else if (strncmp(pdb.type, "END", 3) == 0) /* should differ from ENDMDL */
      return (kEOF);                           /*end */

    else /* meaningless line, so get next one */
    {
      if (myfgets_(linbuf, kBufSize, infile) == kEOF)
        return (kEOF);
    }

  } while (TRUE);

} /* process file */

/* ------------------------------------------------------------------------------------ */
void Fill_Header(GLYp p, unsigned long flags, unsigned long atomNumStart, char chain, char *name, char *nb, char kind)
{
  unsigned short /*firstLetter,*/ qqq;

  for (qqq = 0; qqq <= kAAnameMaxLength; qqq++)
    (*p).header.aa[qqq] = 0;

  /* I prefer to have the label number left justified */
  //strcpy((*p).header.nb, "    ");
  strcpy_s((*p).header.nb, 8, "    ");
#if (kMaxCharForResNum == 3)
  sprintf((*p).header.nb, "%-3ld", atol(nb));
#else
  //sprintf((*p).header.nb, "%-4ld", atol(nb));
  sprintf_s((*p).header.nb, 8, "%-4ld", atol(nb));
#endif

  (*p).header.nb[kMaxCharForResNum] = 0;

  //strncpy((*p).header.aa, name, 3);
  strncpy_s((*p).header.aa, 12, name, 3);
  (*p).header.aa[3] = 0;
  (*p).header.chain = chain;
  (*p).header.attributes = flags;
  (*p).header.startingAtomNumber = atomNumStart;
  (*p).header.kind = kind;
} /*Fill_Header */
/* ------------------------------------------------------------------------------------ */
/* this function retrieve the coord from an atom according to the type query
    in order to "fill" the structure with the correct atom. Just in case the pdb
    file is not ordered in the same way as ours structures */

static short GetAtom(Point3D *coordPtr, ATOM_TYPEp typePtr, const char query[5], const char atm_type[2], short n /*, unsigned short *Bfactor*/)
{
  short k = -1;

  //strncpy((*typePtr).atm_name, query, 5);
  strncpy_s((*typePtr).atm_name, 6, query, 5);
  (*typePtr).atm_flags = 0;
  do
  {
    if (++k == kMaxAtomBuffer)
    {
      /*ShowMessage("\pAn Expected Atom is missing",FALSE);     */
      (*coordPtr).x = kNoCoord;
      (*coordPtr).y = kNoCoord;
      (*coordPtr).z = kNoCoord;
      /*if (*Bfactor != kDoNotProcessBfactor)
       *Bfactor = 9999;*/
      (*typePtr).atm_Bfactor = 9999;
      (*typePtr).atm_occupancy = 100;
#ifdef INSERTION_CODE
      (*typePtr).res_insertion = ' '; /* Prevent empty res_insertions for missing atoms 2.7.02 __TS */
#endif
      return (1);
    }
  } while ((strncmp(generic_aa[k].atm_name, query, n) != 0));

  if (generic_aa[k].atm_coord.x <= 9999.98)
    (*coordPtr) = generic_aa[k].atm_coord;
  else
  {
    (*coordPtr).x = kNoCoord;
    (*coordPtr).y = kNoCoord;
    (*coordPtr).z = kNoCoord;
  }

  (*typePtr).atm_Bfactor = generic_aa[k].Bfactor;
  (*typePtr).atm_occupancy = generic_aa[k].atm_occupancy;
#ifdef INSERTION_CODE
  (*typePtr).res_insertion = generic_aa[k].res_insertion;
#endif
  /*if (*Bfactor != kDoNotProcessBfactor)
   *Bfactor = Max(*Bfactor,generic_aa_Bfactor[k]>*Bfactor);*/

  return (0);

} /*GetAtom */

/* ------------------------------------------------------------------------------------ */
static short Read_aa(unsigned short i)
{
  short err;
  short index;
  short c, j, x;
  char atmname[5];
  char atmtype[2];
  short nbatms;
  ATOM_TYPEp AtTypePtr;

  short k = 0;
  GLYp p;
  long l;

  index = -1;
  do
  {
    if (index == 20)
      return (99);
  } while (strncmp(aaName[++index].aa, currentName, 3) != 0);
  nbatms = 4 + aaName[index].nbSideAtm;

  /* first of all: allocate memory */
  l = sizeof(AtomHeader);
  if (((l * 8) % 64) != 0)
  {
    ShowMessageCode("HEADER not on 64 bits boundaries =", l, true);
  }
  l += (long)nbatms * sizeof(Point3D);
  l += (long)aaName[index].nbPolarH * sizeof(Point3D); /* just for testing purposes */
  GroupsHdl[i] = (GLYp) malloc(l);
  if (GroupsHdl[i] == NULL)
    return (kOverflow);

  p = GroupsHdl[i];

  /* just for testing purposes */
  AtomNamesHdl[i] = (ATOM_TYPEp) calloc((long)(aaName[index].nbPolarH + nbatms), sizeof(ATOM_TYPE));
  if (AtomNamesHdl[i] == NULL)
    return (kOverflow);

  AtTypePtr = AtomNamesHdl[i];
  /* then fill header and retrieve backbone atoms infos */

  err = 0;

  Fill_Header((GLYp)p, attributes, currentAtomNumberStart, currentChain, currentName, currentLabel, aaName[index].kind);

  err += GetAtom(&(*p).N, (AtTypePtr)++, " N  ", "N ", 3);
  if (nbatms > 4) /* not GLY */
    err += GetAtom(&(*p).CA, (AtTypePtr)++, " CA ", "CH", 3);
  else
    err += GetAtom(&(*p).CA, (AtTypePtr)++, " CA ", "C2", 3);

  err += GetAtom(&(*p).C, (AtTypePtr)++, " C  ", "C ", 3);
  err += GetAtom(&(*p).O, (AtTypePtr)++, " O  ", "O ", 3);

  if (err > 0) /* discard the group */
  {
    free(AtomNamesHdl[i]);
    free(GroupsHdl[i]);
    GroupsHdl[i] = NULL;
    AtomNamesHdl[i] = NULL;
    return (0);
  }

  PDB[layer].nbAtomGroups++;

  (*p).header.nbAtoms = nbatms;
  (*p).header.kind = aaName[index].kind;

  if (nbatms == 4) /* GLY */
    goto bail;

  /* now load sidechain */
  err = 0;

  c = 0;
  j = 4;
  do
  {
    atmname[0] = aaName[index].sidechainTypes[c++];
    atmname[1] = aaName[index].sidechainTypes[c++];
    atmname[2] = aaName[index].sidechainTypes[c++];
    atmname[3] = aaName[index].sidechainTypes[c++];
    atmname[4] = '\0';
    atmtype[0] = aaName[index].sidechainTypes[c++];
    atmtype[1] = aaName[index].sidechainTypes[c++];
    if (atmname[3] == ' ')
      x = 3;
    else
    {
      if (atmname[1] <= 'Z')
        x = 4;
      else
      {
        atmname[1] -= 'a';
        atmname[1] += 'A';
        x = 3;
      }
    }
    err += GetAtom(&(*(HETATMp)p).coord[j++], (AtTypePtr)++, atmname, atmtype, x);
  } while (j < nbatms);

bail:

  return (0);

} /*    Read_aa */

/* ------------------------------------------------------------------------------------ */

static short Read_NT_Backbone(unsigned short i, long AGsize, long ATsize, GLYp *p, ATOM_TYPEp *AtTypePtr)
{
  short err = 0;
  short k = 0;

  /* first of all: allocate memory */

  GroupsHdl[i] = (GLYp) malloc(AGsize);
  if (GroupsHdl[i] == NULL)
    return (kOverflow);

  *p = GroupsHdl[i];

  AtomNamesHdl[i] = (ATOM_TYPEp) calloc(1, ATsize);
  if (AtomNamesHdl[i] == NULL)
    return (kOverflow);

  *AtTypePtr = AtomNamesHdl[i];

  /* then fill header and retrieve backbone atoms infos */
  Fill_Header((GLYp)(*p), 0, currentAtomNumberStart, currentChain, currentName, currentLabel, 'n');

  err += GetAtom(&(*(nt_Ap)(*p)).P, (*AtTypePtr)++, " P  ", "P ", 3);

  if (GetAtom(&(*(nt_Ap)(*p)).O1P, (*AtTypePtr), " OP1", "O2", 4))        /* test new nomenclature */
    err += GetAtom(&(*(nt_Ap)(*p)).O1P, (*AtTypePtr)++, " O1P", "O2", 4); /* else try old one */
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).O2P, (*AtTypePtr), " OP2", "O2", 4)) /* test new nomenclature */
    err += GetAtom(&(*(nt_Ap)(*p)).O2P, (*AtTypePtr)++, " O2P", "O2", 4);
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).O5_, (*AtTypePtr), " O5'", "OS", 4))        /* in fact: OS if connected, Otherwise OH */
    err += GetAtom(&(*(nt_Ap)(*p)).O5_, (*AtTypePtr)++, " O5*", "OS", 4); /* in fact: OS if connected, Otherwise OH */
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).C5_, (*AtTypePtr), " C5'", "C2", 4))
    err += GetAtom(&(*(nt_Ap)(*p)).C5_, (*AtTypePtr)++, " C5*", "C2", 4);
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).C4_, (*AtTypePtr), " C4'", "CH", 4))
    err += GetAtom(&(*(nt_Ap)(*p)).C4_, (*AtTypePtr)++, " C4*", "CH", 4);
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).O4_, (*AtTypePtr), " O4'", "OS", 4))
    err += GetAtom(&(*(nt_Ap)(*p)).O4_, (*AtTypePtr)++, " O4*", "OS", 4);
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).C3_, (*AtTypePtr), " C3'", "CH", 4))
    err += GetAtom(&(*(nt_Ap)(*p)).C3_, (*AtTypePtr)++, " C3*", "CH", 4);
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).O3_, (*AtTypePtr), " O3'", "OS", 4))
    err += GetAtom(&(*(nt_Ap)(*p)).O3_, (*AtTypePtr)++, " O3*", "OS", 4);
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).C2_, (*AtTypePtr), " C2'", "C2", 4))
    err += GetAtom(&(*(nt_Ap)(*p)).C2_, (*AtTypePtr)++, " C2*", "C2", 4);
  else
    (*AtTypePtr)++;

  if (GetAtom(&(*(nt_Ap)(*p)).C1_, (*AtTypePtr), " C1'", "CH", 4))
    err += GetAtom(&(*(nt_Ap)(*p)).C1_, (*AtTypePtr)++, " C1*", "CH", 4);
  else
    (*AtTypePtr)++;

  if (err > 0) /* discard the group */
  {
    free(AtomNamesHdl[i]);
    free(GroupsHdl[i]);
    GroupsHdl[i] = NULL;
    AtomNamesHdl[i] = NULL;
    return (0);
  }

  PDB[layer].nbAtomGroups++;
  return (0);

} /*Read_NT_Backbone */

/* ------------------------------------------------------------------------------------ */

static short Read_A(nt_Ap p, ATOM_TYPEp *AtTypePtr, char dnaflag)
{
  short err = 0;

  (*p).header.kind = kA;
  /*    (*p).header.oneLetterCode = 'a';*/
  (*p).header.nbAtoms = kNTA_nbatms;
  err += GetAtom(&(*p).N9, (*AtTypePtr)++, " N9 ", "N*", 4);
  err += GetAtom(&(*p).C8, (*AtTypePtr)++, " C8 ", "CE", 4);
  err += GetAtom(&(*p).N7, (*AtTypePtr)++, " N7 ", "NB", 4);
  err += GetAtom(&(*p).C5, (*AtTypePtr)++, " C5 ", "CB", 4);
  err += GetAtom(&(*p).C6, (*AtTypePtr)++, " C6 ", "CA", 4);
  err += GetAtom(&(*p).N6, (*AtTypePtr)++, " N6 ", "N2", 4);
  err += GetAtom(&(*p).N1, (*AtTypePtr)++, " N1 ", "NC", 4);
  err += GetAtom(&(*p).C2, (*AtTypePtr)++, " C2 ", "CI", 4);
  err += GetAtom(&(*p).N3, (*AtTypePtr)++, " N3 ", "NC", 4);
  err += GetAtom(&(*p).C4, (*AtTypePtr)++, " C4 ", "CB", 4);
  if (dnaflag != 'D') /* attempt to load RNA */
  {
    if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2'", "OH", 4))   /* if no O2' (new nomenclature) we might deal with an old RNA */
      if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2*", "OH", 4)) /* if no O2* we definitely deal with DNA */
        (*p).header.nbAtoms--;
  }
  else
    (*p).header.nbAtoms--;

  return (err);

} /*Read_A */

/* ------------------------------------------------------------------------------------ */

static short Read_G(nt_Gp p, ATOM_TYPEp *AtTypePtr, char dnaflag)
{
  short err = 0;

  (*p).header.kind = kG;
  /*    (*p).header.oneLetterCode = 'g';*/
  (*p).header.nbAtoms = kNTG_nbatms;
  err += GetAtom(&(*p).N9, (*AtTypePtr)++, " N9 ", "N*", 4);
  err += GetAtom(&(*p).C8, (*AtTypePtr)++, " C8 ", "CE", 4);
  err += GetAtom(&(*p).N7, (*AtTypePtr)++, " N7 ", "NB", 4);
  err += GetAtom(&(*p).C5, (*AtTypePtr)++, " C5 ", "CB", 4);
  err += GetAtom(&(*p).C6, (*AtTypePtr)++, " C6 ", "CA", 4);
  err += GetAtom(&(*p).O6, (*AtTypePtr)++, " O6 ", "O ", 4);
  err += GetAtom(&(*p).N1, (*AtTypePtr)++, " N1 ", "NA", 4);
  err += GetAtom(&(*p).C2, (*AtTypePtr)++, " C2 ", "CA", 4);
  err += GetAtom(&(*p).N2, (*AtTypePtr)++, " N2 ", "N2", 4);
  err += GetAtom(&(*p).N3, (*AtTypePtr)++, " N3 ", "NC", 4);
  err += GetAtom(&(*p).C4, (*AtTypePtr)++, " C4 ", "CB", 4);

  if (dnaflag != 'D') /* attempt to load RNA */
  {
    if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2'", "OH", 4))   /* if no O2' (new nomenclature) we might deal with an old RNA */
      if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2*", "OH", 4)) /* if no O2* we definitely deal with DNA */
        (*p).header.nbAtoms--;
  }
  else
    (*p).header.nbAtoms--;

  return (err);

} /*Read_G */

/* ------------------------------------------------------------------------------------ */

static short Read_T(nt_Tp p, ATOM_TYPEp *AtTypePtr, char dnaflag)
{
  short err = 0;

  (*p).header.kind = kT;
  /*    (*p).header.oneLetterCode = 't';*/
  (*p).header.nbAtoms = kNTT_nbatms /*-1*/;
  err += GetAtom(&(*p).N1, (*AtTypePtr)++, " N1 ", "N*", 4);
  err += GetAtom(&(*p).C2, (*AtTypePtr)++, " C2 ", "C ", 4);
  err += GetAtom(&(*p).O2, (*AtTypePtr)++, " O2 ", "O ", 4);
  err += GetAtom(&(*p).N3, (*AtTypePtr)++, " N3 ", "NA", 4);
  err += GetAtom(&(*p).C4, (*AtTypePtr)++, " C4 ", "C ", 4);
  err += GetAtom(&(*p).O4, (*AtTypePtr)++, " O4 ", "O ", 4);
  err += GetAtom(&(*p).C5, (*AtTypePtr)++, " C5 ", "CM", 4);
  if (GetAtom(&(*p).C5M, (*AtTypePtr)++, " C7 ", "C3", 4)) /* new nomenclature */
    err += GetAtom(&(*p).C5M, (*AtTypePtr)++, " C5M", "C3", 4);
  err += GetAtom(&(*p).C6, (*AtTypePtr)++, " C6 ", "CJ", 4);

  if (dnaflag != 'D') /* attempt to load RNA */
  {
    if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2'", "OH", 4))   /* if no O2' (new nomenclature) we might deal with an old RNA */
      if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2*", "OH", 4)) /* if no O2* we definitely deal with DNA */
        (*p).header.nbAtoms--;
  }
  else
    (*p).header.nbAtoms--;

  return (err);

} /*Read_T */

/* ------------------------------------------------------------------------------------ */

static short Read_C(nt_Cp p, ATOM_TYPEp *AtTypePtr, char dnaflag)
{
  short err = 0;

  (*p).header.kind = kC;
  /*    (*p).header.oneLetterCode = 'c';*/
  (*p).header.nbAtoms = kNTC_nbatms;
  err += GetAtom(&(*p).N1, (*AtTypePtr)++, " N1 ", "N*", 4);
  err += GetAtom(&(*p).C2, (*AtTypePtr)++, " C2 ", "C ", 4);
  err += GetAtom(&(*p).O2, (*AtTypePtr)++, " O2 ", "O ", 4);
  err += GetAtom(&(*p).N3, (*AtTypePtr)++, " N3 ", "NC", 4);
  err += GetAtom(&(*p).N4, (*AtTypePtr)++, " N4 ", "N2", 4);
  err += GetAtom(&(*p).C4, (*AtTypePtr)++, " C4 ", "CA", 4);
  err += GetAtom(&(*p).C5, (*AtTypePtr)++, " C5 ", "CJ", 4);
  err += GetAtom(&(*p).C6, (*AtTypePtr)++, " C6 ", "CJ", 4);
  if (dnaflag != 'D') /* attempt to load RNA */
  {
    if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2'", "OH", 4))   /* if no O2' (new nomenclature) we might deal with an old RNA */
      if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2*", "OH", 4)) /* if no O2* we definitely deal with DNA */
        (*p).header.nbAtoms--;
  }
  else
    (*p).header.nbAtoms--;

  return (err);

} /*Read_C */

/* ------------------------------------------------------------------------------------ */

static short Read_U(nt_Up p, ATOM_TYPEp *AtTypePtr, char dnaflag)
{
  short err = 0;

  (*p).header.kind = kU;
  /*    (*p).header.oneLetterCode = 'u';*/
  (*p).header.nbAtoms = kNTU_nbatms;
  err += GetAtom(&(*p).N1, (*AtTypePtr)++, " N1 ", "N*", 4);
  err += GetAtom(&(*p).C2, (*AtTypePtr)++, " C2 ", "C ", 4);
  err += GetAtom(&(*p).O2, (*AtTypePtr)++, " O2 ", "O ", 4);
  err += GetAtom(&(*p).N3, (*AtTypePtr)++, " N3 ", "NA", 4);
  err += GetAtom(&(*p).C4, (*AtTypePtr)++, " C4 ", "C ", 4);
  err += GetAtom(&(*p).O4, (*AtTypePtr)++, " O4 ", "O ", 4);
  err += GetAtom(&(*p).C5, (*AtTypePtr)++, " C5 ", "CM", 4);
  err += GetAtom(&(*p).C6, (*AtTypePtr)++, " C6 ", "CJ", 4);

  if (dnaflag != 'D') /* attempt to load RNA, should always be the case... */
  {
    if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2'", "OH", 4))   /* if no O2' (new nomenclature) we might deal with an old RNA */
      if (GetAtom(&(*p).O2_, (*AtTypePtr), " O2*", "OH", 4)) /* if no O2* we definitely deal with DNA */
        (*p).header.nbAtoms--;
  }
  else
    (*p).header.nbAtoms--;

  return (err);

} /*Read_U */

/* ------------------------------------------------------------------------------------ */
