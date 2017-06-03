EXECS=prog_1
MPICC?=mpicc

all: ${EXECS}

prog_1: prog_1.c
	${MPICC} -o prog_1 -g prog_1.c

clean:
	rm -f ${EXECS}
