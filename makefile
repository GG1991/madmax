###################################################
# Makefile utilities
# make --quite : prints only "echos" instruction
#

PETSC = 1
SLEPC = 1
PARMETIS = 1

SRC_DIR= ./src
OBJ_DIR= ./obj
DEP_DIR= ./inc
SPU_DEP_DIR= ./inc

SPU_INC_DIR= ./inc
SPU_OBJ_DIR= ./obj
SPU_SRC_DIR= ./src

MAC_DIR= macro
MAC_OBJ_DIR= ${MAC_DIR}/obj
MAC_SRC_DIR= ${MAC_DIR}/src
MAC_INC_DIR= ${MAC_DIR}/inc

MIC_DIR= micro
MIC_OBJ_DIR= ${MIC_DIR}/obj
MIC_SRC_DIR= ${MIC_DIR}/src
MIC_INC_DIR= ${MIC_DIR}/inc

OPT = 0
ifeq ($(OPT),1)
CFLAGS=-O3
else
CFLAGS=-g -O0 -lm -std=gnu11
endif

DEPS_SPUTNIK = ${DEP_DIR}/list.h \
	       ${DEP_DIR}/function.h \
	       ${DEP_DIR}/trace.h \
	       ${DEP_DIR}/gmsh.h \
	       ${DEP_DIR}/vtk.h \
	       ${DEP_DIR}/material.h

DEPS_MAC = ${MAC_INC_DIR}/macro.h \
           ${DEP_DIR}/solvers.h \
           ${DEP_DIR}/util.h \
           ${DEP_DIR}/gmsh.h \
           ${DEP_DIR}/comm.h \
	   ${DEP_DIR}/mesh.h \
           ${DEP_DIR}/material.h \
           ${DEP_DIR}/vtk.h \
           ${DEP_DIR}/myio.h \
           ${DEP_DIR}/function.h \
           ${DEP_DIR}/trace.h

DEPS_MIC = ${MIC_INC_DIR}/micro.h \
	   ${MIC_INC_DIR}/micro_struct.h \
           ${DEP_DIR}/comm.h \
           ${DEP_DIR}/solvers.h \
           ${DEP_DIR}/material.h \
           ${DEP_DIR}/myio.h \
           ${DEP_DIR}/vtk.h \
           ${DEP_DIR}/geometry.h \
           ${DEP_DIR}/trace.h

MAC_OBJ  = ${MAC_OBJ_DIR}/main.o \
           $(MAC_OBJ_DIR)/assembly.o \
           $(MAC_OBJ_DIR)/alloc.o \
           $(MAC_OBJ_DIR)/init.o \
           $(MAC_OBJ_DIR)/finalize.o \
           $(MAC_OBJ_DIR)/output.o \
           $(MAC_OBJ_DIR)/boundary.o \
           $(MAC_OBJ_DIR)/comm_line.o \
           $(SPU_OBJ_DIR)/geometry.o \
           $(SPU_OBJ_DIR)/material.o \
           $(SPU_OBJ_DIR)/solvers.o \
           $(SPU_OBJ_DIR)/gmsh.o \
           $(SPU_OBJ_DIR)/comm.o \
           $(SPU_OBJ_DIR)/util.o \
           $(SPU_OBJ_DIR)/list.o \
           $(SPU_OBJ_DIR)/mesh.o \
           $(SPU_OBJ_DIR)/fem.o \
           $(SPU_OBJ_DIR)/function.o \
           $(SPU_OBJ_DIR)/myio.o \
           $(SPU_OBJ_DIR)/vtk.o

MIC_OBJ  = ${MIC_OBJ_DIR}/main.o \
           ${MIC_OBJ_DIR}/homogenize.o \
           ${MIC_OBJ_DIR}/micro_struct.o \
           ${MIC_OBJ_DIR}/init.o \
           ${MIC_OBJ_DIR}/finalize.o \
           ${MIC_OBJ_DIR}/alloc.o \
           ${MIC_OBJ_DIR}/assembly.o \
           ${MIC_OBJ_DIR}/output.o \
           ${MIC_OBJ_DIR}/comm_line.o \
           ${SPU_OBJ_DIR}/util.o \
           ${SPU_OBJ_DIR}/ell.o \
           ${SPU_OBJ_DIR}/mesh.o \
           $(SPU_OBJ_DIR)/comm.o \
           $(SPU_OBJ_DIR)/trace.o \
           $(SPU_OBJ_DIR)/material.o \
           $(SPU_OBJ_DIR)/vtk.o \
           $(SPU_OBJ_DIR)/solvers.o \
           $(SPU_OBJ_DIR)/myio.o \
           $(SPU_OBJ_DIR)/geometry.o \
           ${SPU_OBJ_DIR}/fem.o \
           $(SPU_OBJ_DIR)/list.o

##############################
# PARMETIS VARIABLES

PARMETIS_DIR = ${HOME}/libs/parmetis-4.0.3

PARMETIS_INC = -I${PARMETIS_DIR}/include -I${PARMETIS_DIR}/metis/include

PARMETIS_HEA = ${PARMETIS_DIR}/include/parmetis.h ${PARMETIS_DIR}/metis/include/metis.h

INC_FLAG = -I${DEP_DIR} ${PARMETIS_INC} -I${SPU_INC_DIR}  -I${MAC_INC_DIR} -I${MIC_INC_DIR} 

.PHONY: clean_ all

all: ${MAC_DIR}/macro ${MIC_DIR}/micro

MAC_LDFLAG = ${HOME}/libs/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a \
             ${HOME}/libs/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a       \
             ${HOME}/libs/slepc-3.7.4/arch-linux-opt/lib/libslepc.so                  \
	     -lgsl -lgslcblas -lm

${MAC_DIR}/macro: ${MAC_OBJ}
	gcc -o ${MAC_DIR}/macro $^ ${PETSC_KSP_LIB} ${MAC_LDFLAG} ${SLEPC_EPS_LIB}
	@echo "MACRO great :) !" 

MIC_LDFLAG = -lgsl -lgslcblas -lm

${MIC_DIR}/micro: ${MIC_OBJ}
	gcc -o ${MIC_DIR}/micro $^ ${PETSC_KSP_LIB} ${MIC_LDFLAG}
	@echo "MICRO great :) !" 

${SPU_OBJ_DIR}/%.o: ${SPU_SRC_DIR}/%.c ${DEPS_SPUTNIK} ${PARMETIS_HEA} ${SLEPC_EPS_LIB}
	${PETSC_COMPILE} -DPARMETIS -DPETSC -DMPI -c ${CFLAGS} ${INC_FLAG}  -o $@ $< 
	@echo ">>> "$@

${MAC_OBJ_DIR}/%.o: ${MAC_SRC_DIR}/%.c ${DEPS_MAC}
	${PETSC_COMPILE} -DPARMETIS -DSLEPC -DPETSC -c ${CFLAGS} ${INC_FLAG} -o $@ $< ${SLEPC_EPS_LIB}
	@echo ">>> "$@

${MIC_OBJ_DIR}/%.o: ${MIC_SRC_DIR}/%.c ${DEPS_MIC}
	${PETSC_COMPILE} -c -DPETSC ${CFLAGS} ${INC_FLAG} -o $@ $<
	@echo ">>> "$@

vars:
	@echo "DEPS     = " ${DEPS}
	@echo "INC_FLAG = " ${INC_FLAG}
	@echo "HOME     = " ${HOME}
	@echo "LDFLAG   = " ${LDFLAG}
	@echo "PETSC_COMPILE = " ${PETSC_COMPILE}

clean_:	    
	rm -f obj/* macro/obj/* micro/obj/* macro/macro micro/micro

#include ${PETSC_DIR}/lib/petsc/conf/variables	
#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
