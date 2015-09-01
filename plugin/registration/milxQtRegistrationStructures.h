// This file contains all the structures and constants shared between most files of the plugin
#ifndef MILXQTRegistrationStructures_H
#define MILXQTRegistrationStructures_H

#include <stdio.h>

// Type of the registration F3D or Aladin
typedef enum
{
    F3D,
    Aladin,
    None
} RegType;

// Parameters for F3D algorithm
typedef struct {
    char referenceName[FILENAME_MAX + 1];
    char floatingName[FILENAME_MAX + 1];
    char outputControlPointGridName[FILENAME_MAX + 1];
    char outputWarpedName[FILENAME_MAX + 1];
    int maxit;  /* maxit should be -1 by default*/
    float spacing[3]; /* spacing[0] = sx, spacing[1] = sy, spacing[2] = sz */
    unsigned int ln;
    unsigned int lp;
    bool nopy;
    bool useSym;
    bool cpp2Def;
    char defOutputName[FILENAME_MAX + 1];
} ParamsF3D;


// Parameters for aladin algorithm
typedef struct {
    char referenceName[FILENAME_MAX + 1];
    char floatingName[FILENAME_MAX + 1];
    char outputResultName[FILENAME_MAX + 1];
    bool rigOnly; /* false by default (rigid + Aladin) */
    bool aF3Direct; /* false by default (default is rigid then Aladin) */
    int maxit;  /* maxit should be -1 by default*/
    unsigned int ln;
    unsigned int lp;
    bool useSym;
    float percentBlock; /* Percentage of block to use, default 50 */
} ParamsAladin;


#endif // MILXQTRegistrationStructures_H