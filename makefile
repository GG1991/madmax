###################################################
# Makefile utilities
# make --quite : prints only "echos" instruction
#
export PETSC_DIR=/home/guido/libs/petsc-3.7.5
export PETSC_ARCH=arch-linux2-c-debug

PWD:= $(shell pwd)

SRC_DIR= ./src
OBJ_DIR= ./obj
DEP_DIR= ./inc

CFLAGS=-g -O0 
	
DEPS = ${DEP_DIR}/macro.h       \
       ${DEP_DIR}/list.h

OBJ  = ${OBJ_DIR}/spu_parser.o    
       #${OBJ_DIR}/spu_mesh.o      
       #${OBJ_DIR}/list.o

##############################
# PARMETIS VARIABLES

PARMETIS_DIR = ${HOME}/libs/parmetis-4.0.3

PARMETIS_INC = -I${PARMETIS_DIR}/include       \
               -I${PARMETIS_DIR}/metis/include

PARMETIS_HEA = ${PARMETIS_DIR}/include/parmetis.h ${PARMETIS_DIR}/metis/include/metis.h

INC = -I${DEP_DIR}

.PHONY: clean_

##############################
# LINK
all: ${OBJ} 
	${MAKE} -C macro 


##############################
# SPU_PARSE.O
obj/spu_parser.o: src/spu_parser.c inc/sputnik.h
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I./inc
	@echo "spu_parser.o" 


##############################
# SPU_MESH.O
obj/spu_mesh.o: src/spu_mesh.c ${DEPS} ${PARMETIS_HEA}
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $<  ${INC}  ${PARMETIS_INC} 
	@echo "spu_mesh.o" 


##############################
# LIST.O
obj/list.o: src/list.c inc/list.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc

vars:
	@echo "DEPS = " ${DEPS}
	@echo "INC = " ${INC}
	@echo "PWD = " ${PWD}
	@echo "HOME = " ${HOME}
	@echo "LDFLAG = " ${LDFLAG}
	@echo "PETSC_COMPILE = " ${PETSC_COMPILE}

clean_:	    
	rm -f obj/* 

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${PETSC_DIR}/lib/petsc/conf/rules
