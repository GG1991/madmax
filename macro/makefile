
export PETSC_DIR=/home/guido/libs/petsc-3.7.5
export PETSC_ARCH=arch-linux2-c-debug

PWD:= $(shell pwd)

SRC_DIR= ./src
OBJ_DIR= ./obj
DEP_DIR= ./inc
SPU_OBJ_DIR= ../obj
SPU_DEP_DIR= ../inc

CFLAGS=-g -O0 
	
DEPS = ${DEP_DIR}/macro.h          \
       ${SPU_DEP_DIR}/sputnik.h

DEPS_DIR = ${DEP_DIR}              \
           ${SPU_DEP_DIR}

OBJ  = ${OBJ_DIR}/mac_main.o       \
       ${OBJ_DIR}/mac_color.o      \
       ${SPU_OBJ_DIR}/spu_parser.o \
       ${SPU_OBJ_DIR}/spu_mesh.o   \
       ${SPU_OBJ_DIR}/list.o  

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
	gcc -o macro $^ ${PETSC_KSP_LIB} -lgsl -lgslcblas -lm ${LDFLAG}
	@echo "macro.o great :) !" 


##############################
# MAC_MAIN.O
obj/mac_main.o: src/mac_main.c ${DEPS}
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I${DEP_DIR} -I${SPU_DEP_DIR}
	@echo "mac_main.o" 


##############################
# MAC_UTIL.O
obj/mac_util.o: src/mac_util.c inc/macro.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc 
	@echo "mac_util.o" 



##############################
# MAC_COLOR.O
obj/mac_color.o: src/mac_color.c ../inc/sputnik.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I../inc 
	@echo "mac_color.o" 


clean_:	    
	rm -f obj/* macro

vars:
	echo ${DEPS}

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${PETSC_DIR}/lib/petsc/conf/rules


	
