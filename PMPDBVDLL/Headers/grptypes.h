//
//  grptypes.h
//  PDBV
//
//  Created by Christian Iseli on 22.09.2024.
//
#pragma once

#include "PMPDBVDLL.h"

extern "C" {

#ifndef grptypes_h
#define grptypes_h

enum {
  kALA = 'A',
  kCYS = 'C',
  kASP = 'D',
  kGLU = 'E',
  kPHE = 'F',
  kGLY = 'G',
  kHIS = 'H',
  kILE = 'I',
  kLYS = 'K',
  kLEU = 'L',
  kMET = 'M',
  kASN = 'N',
  kPRO = 'P',
  kGLN = 'Q',
  kARG = 'R',
  kSER = 'S',
  kTHR = 'T',
  kVAL = 'V',
  kTRP = 'W',
  kHETaa = 'X',
  kTYR = 'Y',
  kG = 'g',
  kA = 'a',
  kT = 't',
  kC = 'c',
  kU = 'u'
};

typedef struct aaDef_struct aaDef;
struct aaDef_struct
{
  char aa[4];
  char kind;
  short nbSideAtm;
  short nbPolarH;
  char sidechainTypes[104];
  Boolean sideCanHbond;
};

/* AA kind < NT kind < HETATM kind == 0xff */
#define kHETkind 0x7f
#define kAAkind kTYR /* last aa */
#define kNTkind kU   /* last nt */

/*------------------ PDB line ------------------- */
typedef struct PDB_line_struct PDB_line;
struct PDB_line_struct
{
  Point3D coord;
  char type[8];
  char atomName[8];
  char labelName[8];
  char labelNumber[8];
  unsigned short Bfactor; /* multiplied by 100 */
  unsigned short atm_occupancy;
  int atomNumber;
  char res_insertion; /* __ts */
};

/*--------------------- HEADER ----------------------- */
/* WARNING: must be aligned on 64bits boundaries on PC */
/*---------------------------------------------------- */
#define kAAnameMaxLength 11

typedef struct AtomHeader_struct AtomHeader;
struct AtomHeader_struct
{
  unsigned long startingAtomNumber;
  unsigned long attributes;
  unsigned short nbAtoms;
  unsigned short label_seq_id;
  char label_entity_id[12];
  char label_asym_id[12];
  char aa[12];
  char nb[8];
  char chain;
  unsigned char kind;
  char pad[4];
};

typedef struct ATOM_TYPE_struct ATOM_TYPE;
struct ATOM_TYPE_struct
{
  char atm_name[6];
  unsigned short atm_Bfactor;
  unsigned short atm_occupancy;
  unsigned short atm_flags;
  char res_insertion; /* __ts */
};
typedef ATOM_TYPE *ATOM_TYPEp;
typedef ATOM_TYPEp *ATOM_TYPEh;

/* -------------- Amino Acids ---------------------- */

struct ALA_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D H;
};

struct CYS_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D SG;
  Point3D H;
  Point3D HG;
};

struct ASP_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D OD1;
  Point3D OD2;
  Point3D H;
};

struct GLU_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD;
  Point3D OE1;
  Point3D OE2;
  Point3D H;
};

struct PHE_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD1;
  Point3D CD2;
  Point3D CE1;
  Point3D CE2;
  Point3D CZ;
  Point3D H;
  Point3D HD1;
  Point3D HD2;
  Point3D HE1;
  Point3D HE2;
  Point3D HZ;
};

struct GLY_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D H;
};

struct HIS_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D ND1;
  Point3D CD2;
  Point3D CE1;
  Point3D NE2;
  Point3D H;
  Point3D HD1;
};

struct ILE_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG1;
  Point3D CG2;
  Point3D CD1;
  Point3D H;
};

struct LYS_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD;
  Point3D CE;
  Point3D NZ;
  Point3D H;
  Point3D HZ1;
  Point3D HZ2;
  Point3D HZ3;
};

struct LEU_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD1;
  Point3D CD2;
  Point3D H;
};

struct MET_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D SD;
  Point3D CE;
  Point3D H;
};

struct ASN_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D OD1;
  Point3D ND2;
  Point3D H;
  Point3D HD21;
  Point3D HD22;
};

struct PRO_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD;
};

struct GLN_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD;
  Point3D OE1;
  Point3D NE2;
  Point3D H;
  Point3D HE21;
  Point3D HE22;
};

struct ARG_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD;
  Point3D NE;
  Point3D CZ;
  Point3D NH1;
  Point3D NH2;
  Point3D H;
  Point3D HE;
  Point3D HH11;
  Point3D HH12;
  Point3D HH21;
  Point3D HH22;
};

struct SER_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D OG;
  Point3D H;
  Point3D HG;
};

struct THR_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D OG1;
  Point3D CG2;
  Point3D H;
  Point3D HG1;
};

struct VAL_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG1;
  Point3D CG2;
  Point3D H;
};

struct TRP_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD1;
  Point3D CD2;
  Point3D NE1;
  Point3D CE2;
  Point3D CE3;
  Point3D CZ2;
  Point3D CZ3;
  Point3D CH2;
  Point3D H;
  Point3D HD1;
  Point3D HE1;
  Point3D HE3;
  Point3D HZ2;
  Point3D HZ3;
  Point3D HH2;
};

struct TYR_struct
{
  AtomHeader header;
  Point3D N;
  Point3D CA;
  Point3D C;
  Point3D O;
  Point3D CB;
  Point3D CG;
  Point3D CD1;
  Point3D CD2;
  Point3D CE1;
  Point3D CE2;
  Point3D CZ;
  Point3D OH;
  Point3D H;
  Point3D HD1;
  Point3D HD2;
  Point3D HE1;
  Point3D HE2;
  Point3D HH;
};

#define kNTA_nbatms 22
struct ntA_struct
{
  AtomHeader header;
  Point3D P;
  Point3D O1P;
  Point3D O2P;
  Point3D O5_;
  Point3D C5_;
  Point3D C4_;
  Point3D O4_;
  Point3D C3_;
  Point3D O3_;
  Point3D C2_;
  Point3D C1_;
  Point3D N9;
  Point3D C8;
  Point3D N7;
  Point3D C5;
  Point3D C6;
  Point3D N6;
  Point3D N1;
  Point3D C2;
  Point3D N3;
  Point3D C4;
  Point3D O2_;
};

#define kNTG_nbatms 23
struct ntG_struct
{
  AtomHeader header;
  Point3D P;
  Point3D O1P;
  Point3D O2P;
  Point3D O5_;
  Point3D C5_;
  Point3D C4_;
  Point3D O4_;
  Point3D C3_;
  Point3D O3_;
  Point3D C2_;
  Point3D C1_;
  Point3D N9;
  Point3D C8;
  Point3D N7;
  Point3D C5;
  Point3D C6;
  Point3D O6;
  Point3D N1;
  Point3D C2;
  Point3D N2;
  Point3D N3;
  Point3D C4;
  Point3D O2_;
};

#define kNTT_nbatms 21
struct ntT_struct
{
  AtomHeader header;
  Point3D P;
  Point3D O1P;
  Point3D O2P;
  Point3D O5_;
  Point3D C5_;
  Point3D C4_;
  Point3D O4_;
  Point3D C3_;
  Point3D O3_;
  Point3D C2_;
  Point3D C1_;
  Point3D N1;
  Point3D C2;
  Point3D O2;
  Point3D N3;
  Point3D C4;
  Point3D O4;
  Point3D C5;
  Point3D C5M;
  Point3D C6;
  Point3D O2_;
};

#define kNTC_nbatms 20
struct ntC_struct
{
  AtomHeader header;
  Point3D P;
  Point3D O1P;
  Point3D O2P;
  Point3D O5_;
  Point3D C5_;
  Point3D C4_;
  Point3D O4_;
  Point3D C3_;
  Point3D O3_;
  Point3D C2_;
  Point3D C1_;
  Point3D N1;
  Point3D C2;
  Point3D O2;
  Point3D N3;
  Point3D N4;
  Point3D C4;
  Point3D C5;
  Point3D C6;
  Point3D O2_;
};

#define kNTU_nbatms 20
struct ntU_struct
{
  AtomHeader header;
  Point3D P;
  Point3D O1P;
  Point3D O2P;
  Point3D O5_;
  Point3D C5_;
  Point3D C4_;
  Point3D O4_;
  Point3D C3_;
  Point3D O3_;
  Point3D C2_;
  Point3D C1_;
  Point3D N1;
  Point3D C2;
  Point3D O2;
  Point3D N3;
  Point3D C4;
  Point3D O4;
  Point3D C5;
  Point3D C6;
  Point3D O2_;
};

struct HETATM_struct
{
  AtomHeader header;
  Point3D coord[16]; /* liste des points 3D de la structure */
};

typedef struct aa_struct GENERIC_aa;

typedef struct GLY_struct ANY_GROUP;
typedef struct ALA_struct ALA;
typedef struct CYS_struct CYS;
typedef struct ASP_struct ASP;
typedef struct GLU_struct GLU;
typedef struct PHE_struct PHE;
typedef struct GLY_struct GLY;
typedef struct HIS_struct HIS;
typedef struct ILE_struct ILE;
typedef struct LYS_struct LYS;
typedef struct LEU_struct LEU;
typedef struct MET_struct MET;
typedef struct ASN_struct ASN;
typedef struct PRO_struct PRO;
typedef struct GLN_struct GLN;
typedef struct ARG_struct ARG;
typedef struct SER_struct SER;
typedef struct THR_struct THR;
typedef struct VAL_struct VAL;
typedef struct TRP_struct TRP;
typedef struct TYR_struct TYR;
typedef struct ntG_struct nt_G;
typedef struct ntA_struct nt_A;
typedef struct ntT_struct nt_T;
typedef struct ntC_struct nt_C;
typedef struct ntU_struct nt_U;
typedef struct HETATM_struct HETATM;

typedef ALA *ALAp;
typedef CYS *CYSp;
typedef ASP *ASPp;
typedef GLU *GLUp;
typedef PHE *PHEp;
typedef GLY *GLYp;
typedef HIS *HISp;
typedef ILE *ILEp;
typedef LYS *LYSp;
typedef LEU *LEUp;
typedef MET *METp;
typedef ASN *ASNp;
typedef PRO *PROp;
typedef GLN *GLNp;
typedef ARG *ARGp;
typedef SER *SERp;
typedef THR *THRp;
typedef VAL *VALp;
typedef TRP *TRPp;
typedef TYR *TYRp;
typedef nt_G *nt_Gp;
typedef nt_A *nt_Ap;
typedef nt_T *nt_Tp;
typedef nt_C *nt_Cp;
typedef nt_U *nt_Up;
typedef HETATM *HETATMp;

#endif /* grptypes_h */
}
