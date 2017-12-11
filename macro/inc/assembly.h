#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "macro.h"

#ifdef PETSC
#include "petscksp.h"
int assembly_b_petsc(void);
int assembly_AM_petsc(void);
int assembly_A_petsc(void);
#endif

#endif
