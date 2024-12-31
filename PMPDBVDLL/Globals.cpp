//
//  SharedData.c
//  PDBV
//
//  Created by Christian Iseli on 21.09.2024.
//

#include "pch.h"
#include "Globals.h"

short nbLayer;
PDBmodel PDB[maxLayerNb];

aaDef aaName[] = {{"ALA", 'A', 1, 1, " CB C3 H  H ", FALSE},
                  {"ARG", 'R', 7, 6, " CB C2 CG C2 CD C2 NE N2 CZ CA NH1N2 NH2N2 H  H  HE H 1HH1H 2HH1H 1HH2H 2HH2H ", TRUE},
                  {"ASN", 'N', 4, 3, " CB C2 CG C  oD1O  nD2N  H  H 1HD2H 2HD2H ", TRUE},
                  {"ASP", 'D', 4, 1, " CB C2 CG C  OD1O2 OD2O2 H  H ", TRUE},
                  {"CYS", 'C', 2, 2, " CB C2 SG SH H  H  HG H ", TRUE},
                  {"GLN", 'Q', 5, 3, " CB C2 CG C2 CD C  oE1O  nE2N  H  H 1HE2H 2HE2H ", TRUE},
                  {"GLU", 'E', 5, 1, " CB C2 CG C2 CD C  OE1O2 OE2O2 H  H ", TRUE},
                  {"GLY", 'G', 0, 1, " H  H ", FALSE},
                  {"HIS", 'H', 6, 2, " CB C2 CG CC nD1NA cD2CF cE1CP nE2NB H  H  HD1H ", TRUE},
                  {"ILE", 'I', 4, 1, " CB CH CG1C2 CG2C3 cD1C3 H  H ", FALSE},
                  {"LEU", 'L', 4, 1, " CB C2 CG CH CD1C3 CD2C3 H  H ", FALSE},
                  {"LYS", 'K', 5, 4, " CB C2 CG C2 CD C2 CE C2 NZ N3 H  H  HZ1H  HZ2H  HZ3H ", TRUE},
                  {"MET", 'M', 4, 1, " CB C2 CG C2 SD S  CE C3 H  H ", TRUE},
                  {"PHE", 'F', 7, 6, " CB C2 CG CA CD1CD CD2CD CE1CD CE2CD CZ CD H  H  HD1H  HD2H  HE1H  HE2H  HZ H ", FALSE},
                  {"PRO", 'P', 3, 0, " CB C2 CG C2 CD C2", FALSE},
                  {"SER", 'S', 2, 2, " CB C2 OG OH H  H  HG H ", TRUE},
                  {"THR", 'T', 3, 2, " CB CH oG1OH cG2C3 H  H  HG1H ", TRUE},
                  {"TRP", 'W', 10, 7, " CB C2 CG C* CD1CG CD2CB nE1NA CE2CN CE3CD CZ2CD CZ3CD CH2CD H  H  HD1H  HE1H  HE3H  HZ2H  HZ3H  HH2H ", TRUE},
                  {"TYR", 'Y', 8, 6, " CB C2 CG CA CD1CD CD2CD CE1CD CE2CD CZ C  OH OH H  H  HD1H  HD2H  HE1H  HE2H  HH H ", TRUE},
                  {"VAL", 'V', 3, 1, " CB CH CG1C3 CG2C3 H  H ", FALSE}};

/* Pointer to the C++ function that handles printf calls in error/debug situations */
int (*printfDelegate)(const char* fmt, va_list va);
