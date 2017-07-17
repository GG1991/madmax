
matpre: matpre.o
	gcc ${CFLAGS} -o $@ $< ${PETSC_KSP_LIB} 	

matpre.o: matpre.c
	${PETSC_COMPILE} ${CFLAGS} -o $@ $< ${PETSC_KSP_LIB} 	

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${PETSC_DIR}/lib/petsc/conf/rules
