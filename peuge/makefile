
export PETSC_DIR=/home/guido/libs/petsc-3.7.5
export PETSC_ARCH=arch-linux2-c-debug

PWD:= $(shell pwd)

SRC_DIR= ./src
OBJ_DIR= ./obj
DEP_DIR= ./inc                     

CFLAGS=-g -O0 
	
DEPS = ${DEP_DIR}/peuge.h                   



OBJ  = ${OBJ_DIR}/peu_main.o  \
       ${OBJ_DIR}/peu_color.o  \
       ${OBJ_DIR}/list.o  

LDFLAG = -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libparmetis \
         -L../../../libs/parmetis-4.0.3/build/Linux-x86_64/libmetis 

.PHONY: clean_

##############################
# LINK
all: ${OBJ} 
	gcc -o peunt $^ ${PETSC_KSP_LIB} -lgsl -lgslcblas -lm



##############################
# PEU_MAIN.O
obj/peu_main.o: src/peu_main.c inc/peunt.h
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I./inc



##############################
# PEU_MESH.O
obj/peu_mesh.o: src/peu_mesh.c inc/list.h \
                ~/libs/parmetis-4.0.3/include/parmetis.h \
                ~/libs/parmetis-4.0.3/metis/include/metis.h
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc -I../../../libs/parmetis-4.0.3/include -I../../../libs/parmetis-4.0.3/metis/include 

##############################
# PEU_UTIL.O
obj/peu_util.o: src/peu_util.c inc/peuge.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I./inc 



##############################
# PEU_COLOR.O
obj/peu_color.o: src/peu_color.c ../inc/monts.h 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I../inc 



##############################
# PEU_PARSER.O
obj/peu_parser.o: src/peu_parser.c 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< 	-I../inc 



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


	
