#MPICC= /home/guido/libs/openmpi-install/bin/mpicc
MPICC=mpicc
MPIRUN=mpirun
SRC_DIR= ../src
OBJ_DIR= ../obj
HEA_DIR= ../inc

all: unity_1 unity_2

unity_1: unity_1.o $(OBJ_DIR)/myio.o  
	$(MPICC) -o unity_1 $^
	@echo "unity_1 generated..." 

unity_2: unity_2.o $(OBJ_DIR)/myio.o $(OBJ_DIR)/gmsh.o $(OBJ_DIR)/list.o
	$(MPICC) -o unity_2 $^
	@echo "unity_2 generated..." 

unity_1.o: unity_1.c $(HEA_DIR)/myio.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11
	@echo ">>> "$@

unity_2.o: unity_2.c $(HEA_DIR)/gmsh.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11
	@echo ">>> "$@

$(OBJ_DIR)/myio.o: ../src/myio.c $(HEA_DIR)/myio.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -I/usr/lib/openmpi/include -std=gnu11
	@echo ">>> "$@

$(OBJ_DIR)/list.o: ../src/list.c $(HEA_DIR)/list.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11
	@echo ">>> "$@

$(OBJ_DIR)/gmsh.o: ../src/gmsh.c $(HEA_DIR)/gmsh.h $(HEA_DIR)/myio.h $(HEA_DIR)/list.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11
	@echo ">>> "$@

clean:	    
	rm -f $(OBJ_DIR)/* unity_1 unity_1.o

run:
	$(MPIRUN) -np 1 ./unity_2
	$(MPIRUN) -np 2 ./unity_2
	$(MPIRUN) -np 3 ./unity_2
