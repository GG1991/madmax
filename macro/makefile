
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

LDFLAG = -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libparmetis \
         -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libmetis 

.PHONY: clean_

##############################
# LINK
all: ${OBJ} 
	gcc -o macro $^ ${PETSC_KSP_LIB} -lgsl -lgslcblas -lm



##############################
# MAC_MAIN.O
obj/mac_main.o: src/mac_main.c ${DEPS}
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I${DEP_DIR} -I${SPU_DEP_DIR}



##############################
# MAC_MESH.O
obj/mac_mesh.o: src/mac_mesh.c ${SPU_DEP_DIR}/list.h \
                ~/libs/parmetis-4.0.3/include/parmetis.h \
                ~/libs/parmetis-4.0.3/metis/include/metis.h
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc -I../../../libs/parmetis-4.0.3/include -I../../../libs/parmetis-4.0.3/metis/include 

##############################
# MAC_UTIL.O
obj/mac_util.o: src/mac_util.c inc/macro.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc 



##############################
# MAC_COLOR.O
obj/mac_color.o: src/mac_color.c ../inc/sputnik.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I../inc 


clean_:	    
	rm -f obj/* macro

vars:
	echo ${DEPS}

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${PETSC_DIR}/lib/petsc/conf/rules


	
