
export PETSC_DIR=/home/guido/libs/petsc-3.7.5
export PETSC_ARCH=arch-linux2-c-debug

PWD:= $(shell pwd)

SRC_DIR= ./src
OBJ_DIR= ./obj
DEP_DIR= ./inc                     
MON_OBJ_DIR= ../obj
MON_DEP_DIR= ../inc                     

CFLAGS=-g -O0 
	
DEPS = ${DEP_DIR}/giant.h                   

OBJ  = ${OBJ_DIR}/gia_main.o      \
       ${OBJ_DIR}/gia_color.o     \
       ${MON_OBJ_DIR}/mon_parse.o \
       ${MON_OBJ_DIR}/list.o  

LDFLAG = -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libparmetis \
         -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libmetis 

.PHONY: clean_

##############################
# LINK
all: ${OBJ} 
	gcc -o giant $^ ${PETSC_KSP_LIB} -lgsl -lgslcblas -lm



##############################
# GIA_MAIN.O
obj/gia_main.o: src/gia_main.c inc/giant.h
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I./inc -I${MON_DEP_DIR}





##############################
# GIA_UTIL.O
obj/gia_util.o: src/gia_util.c inc/giant.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc 




##############################
# GIA_COLOR.O
obj/gia_color.o: src/gia_color.c ../inc/monts.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I../inc


var:
	echo "PWD = " ${PWD}
	echo "PETSC_COMPILE = " ${PETSC_COMPILE}

clean_:	    
	rm -f obj/* giant

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${PETSC_DIR}/lib/petsc/conf/rules


	
