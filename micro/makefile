
export PETSC_DIR=/home/guido/libs/petsc-3.7.5
export PETSC_ARCH=arch-linux2-c-debug

PWD:= $(shell pwd)

SRC_DIR= ./src
OBJ_DIR= ./obj
DEP_DIR= ./inc
SPU_OBJ_DIR= ../obj
SPU_DEP_DIR= ../inc

CFLAGS=-g -O0 
	
DEPS = ${DEP_DIR}/micro.h          \
       ${SPU_DEP_DIR}/sputnik.h

DEPS_DIR = ${DEP_DIR}              \
           ${SPU_DEP_DIR}

OBJ  = ${OBJ_DIR}/mic_main.o       \
       ${OBJ_DIR}/mic_comm.o       \
       ${SPU_OBJ_DIR}/spu_parser.o 
#       ${SPU_OBJ_DIR}/spu_mesh.o   

.PHONY: clean_

##############################
# PARMETIS VARIABLES

LDFLAG = ${HOME}/libs/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a \
         ${HOME}/libs/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a 

PARMETIS_DIR = ${HOME}/libs/parmetis-4.0.3

PARMETIS_INC = -I${PARMETIS_DIR}/include       \
               -I${PARMETIS_DIR}/metis/include

##############################
# LINK
all: ${OBJ} 
	gcc -o micro $^ ${PETSC_KSP_LIB} -lgsl -lgslcblas -lm ${LDFLAG}
	@echo "micro.o great :) !" 


##############################
# MAC_MAIN.O
obj/mic_main.o: src/mic_main.c ${DEPS}
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I${DEP_DIR} -I${SPU_DEP_DIR}
	@echo "mic_main.o" 


##############################
# MIC_COMM.O
obj/mic_comm.o: src/mic_comm.c ../inc/sputnik.h inc/micro.h
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I../inc  -I./inc
	@echo "mic_comm.o" 


clean_:	    
	rm -f obj/* micro

vars:
	echo ${DEPS}

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${PETSC_DIR}/lib/petsc/conf/rules


	
