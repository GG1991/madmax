EXECS=prog_1 prog_2
MPICC?=mpicc

all: ${EXECS}

prog_1: prog_1.c
	${MPICC} -o prog_1 -g prog_1.c

prog_2: prog_2.c
	${MPICC} -o prog_2 -g prog_2.c

clean:
	rm -f ${EXECS}
