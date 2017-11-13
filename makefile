###################################################
# Makefile utilities
# make --quite : prints only "echos" instruction
#

PWD:= $(shell pwd)

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
CFLAGS=-g -O0
endif

DEPS_SPUTNIK = ${DEP_DIR}/sputnik.h \
	       ${DEP_DIR}/list.h \
	       ${DEP_DIR}/boundary.h \
	       ${DEP_DIR}/fun.h \
	       ${DEP_DIR}/trace.h \
	       ${DEP_DIR}/gmsh.h \
	       ${DEP_DIR}/material.h

DEPS_MACMIC =  ${MAC_INC_DIR}/macro.h \
	       ${MIC_INC_DIR}/micro.h \
	       ${DEP_DIR}/sputnik.h

DEPS_MAC =  ${MIC_INC_DIR}/macro.h \
            ${DEP_DIR}/sputnik.h \
            ${DEP_DIR}/util.h \
            ${DEP_DIR}/gmsh.h \
            ${DEP_DIR}/comm.h \
            ${DEP_DIR}/material.h \
            ${DEP_DIR}/trace.h

DEPS_MIC =  ${MIC_INC_DIR}/micro.h \
            ${DEP_DIR}/comm.h \
            ${DEP_DIR}/material.h \
            ${DEP_DIR}/trace.h

DEP_DIRS = ${DEP_DIR} ${MAC_INC_DIR} ${MIC_INC_DIR}

SPU_OBJ  = $(SPU_OBJ_DIR)/mesh.o     \
           $(SPU_OBJ_DIR)/time.o     \
           $(SPU_OBJ_DIR)/util.o     \
           $(SPU_OBJ_DIR)/trace.o    \
           $(SPU_OBJ_DIR)/material.o \
           $(SPU_OBJ_DIR)/gmsh.o     \
           $(SPU_OBJ_DIR)/comm.o    

MAC_OBJ  = ${MAC_OBJ_DIR}/main.o     \
           $(SPU_OBJ_DIR)/material.o

MIC_OBJ  = ${MIC_OBJ_DIR}/main.o \
           ${MIC_OBJ_DIR}/homogenize.o \
           $(SPU_OBJ_DIR)/comm.o \
           $(SPU_OBJ_DIR)/trace.o \
           $(SPU_OBJ_DIR)/material.o \
           $(SPU_OBJ_DIR)/list.o

EXT_OBJ  = $(SPU_OBJ_DIR)/fem.o \
           $(SPU_OBJ_DIR)/list.o \
           $(SPU_OBJ_DIR)/fun.o 
 
EXT_DEP  = $(SPU_INC_DIR)/fem.h

DEPS+= $(EXT_DEP)

SPU_OBJ+= $(EXT_OBJ)

##############################
# PARMETIS VARIABLES

PARMETIS_DIR = ${HOME}/libs/parmetis-4.0.3

PARMETIS_INC = -I${PARMETIS_DIR}/include \
               -I${PARMETIS_DIR}/metis/include

PARMETIS_HEA = ${PARMETIS_DIR}/include/parmetis.h ${PARMETIS_DIR}/metis/include/metis.h

LDFLAG = ${HOME}/libs/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a \
         ${HOME}/libs/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a       \
         ${HOME}/libs/slepc-3.7.4/arch-linux-opt/lib/libslepc.so                  \
	 -lgsl -lgslcblas 

MIC_LDFLAG = -lgsl -lgslcblas -lm

INC = -I${DEP_DIR}

CFLAGS+= ${PARMETIS_INC} -I${SPU_INC_DIR}  -I${MAC_INC_DIR} -I${MIC_INC_DIR}

.PHONY: clean_ all

##############################
# LINK
all: ${MAC_DIR}/macro ${MIC_DIR}/micro

##############################
# MACRO
${MAC_DIR}/macro: ${MAC_OBJ} ${SPU_OBJ}
	gcc -o ${MAC_DIR}/macro $^ ${PETSC_KSP_LIB} -lm ${LDFLAG} ${SLEPC_EPS_LIB}
	@echo "MACRO great :) !" 

##############################
# MICRO
${MIC_DIR}/micro: ${MIC_OBJ}
	gcc -o ${MIC_DIR}/micro $^ ${PETSC_KSP_LIB} ${MIC_LDFLAG}
	@echo "MICRO great :) !" 

##############################
# SPUTNIK OBJECTS (do not work)
${SPU_OBJ_DIR}/%.o: ${SPU_SRC_DIR}/%.c ${DEPS_SPUTNIK} ${PARMETIS_HEA} ${SLEPC_EPS_LIB}
	${PETSC_COMPILE} -DPETSC -c ${CFLAGS} -o $@ $< 
	@echo ">>> "$@

##############################
# MACRO OBJECTS
${MAC_OBJ_DIR}/%.o: ${MAC_SRC_DIR}/%.c ${DEPS_MACMIC}
	${PETSC_COMPILE} -DSLEPC -DPETSC -c ${CFLAGS} -o $@ $< 	${SLEPC_EPS_LIB}
	@echo ">>> "$@

##############################
# MICRO OBJECTS
${MIC_OBJ_DIR}/%.o: ${MIC_SRC_DIR}/%.c ${DEPS_MIC}
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	${SLEPC_EPS_LIB}
	@echo ">>> "$@

vars:
	@echo "DEPS = " ${DEPS}
	@echo "INC = " ${INC}
	@echo "PWD = " ${PWD}
	@echo "HOME = " ${HOME}
	@echo "LDFLAG = " ${LDFLAG}
	@echo "PETSC_COMPILE = " ${PETSC_COMPILE}

clean_:	    
	rm -f obj/* macro/obj/* micro/obj/* macro/macro micro/micro

#include ${PETSC_DIR}/lib/petsc/conf/variables	
#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
