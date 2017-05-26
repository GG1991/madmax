
export PETSC_DIR=/home/guido/libs/petsc-3.7.5
export PETSC_ARCH=arch-linux2-c-debug

PWD:= $(shell pwd)

SRC_DIR= ./src
OBJ_DIR= ./obj
DEP_DIR= ./inc                     

CFLAGS=-g -O0 
	
DEPS = ${DEP_DIR}/monts.h                   

OBJ  = ${OBJ_DIR}/mon_parse.o    \
       ${OBJ_DIR}/list.o

LDFLAG = -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libparmetis \
         -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libmetis 

.PHONY: clean_

##############################
# LINK
all: ${OBJ} 
#	gcc -o peuge $^ ${PETSC_KSP_LIB} -lgsl -lgslcblas -lm



##############################
# MON_PARSE.O
obj/mon_parse.o: src/mon_parse.c inc/monts.h
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I./inc



##############################
# LIST.O
obj/list.o: src/list.c inc/list.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc




var:
	echo "PWD = " ${PWD}
	echo "PETSC_COMPILE = " ${PETSC_COMPILE}

clean_:	    
	rm -f obj/* peunt

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${PETSC_DIR}/lib/petsc/conf/rules


	
