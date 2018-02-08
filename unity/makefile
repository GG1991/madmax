MPICC=mpicc
MPIRUN=mpirun
SRC_DIR= ../src
OBJ_DIR= ../obj
HEA_DIR= ../inc

all: unity_1 unity_2 unity_3 unity_ell

unity_1: unity_1.o $(OBJ_DIR)/myio.o  
	$(MPICC) -o unity_1 $^
	@echo "unity_1 generated..." 

unity_2: unity_2.o $(OBJ_DIR)/myio.o $(OBJ_DIR)/gmsh.o $(OBJ_DIR)/list.o
	$(MPICC) -o unity_2 $^
	@echo "unity_2 generated..." 

unity_3: unity_3.o $(OBJ_DIR)/util.o
	$(MPICC) -o unity_3 $^ -lm
	@echo "unity_3 generated..." 

unity_4: unity_4.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/list.o $(OBJ_DIR)/myio.o $(OBJ_DIR)/util.o
	$(MPICC) -o unity_4 $^ -lm
	@echo "unity_4 generated..." 

unity_ell: unity_ell.o $(OBJ_DIR)/ell.o
	$(MPICC) -o unity_ell $^ -lm
	@echo "unity_ell generated..." 

unity_1.o: unity_1.c $(HEA_DIR)/myio.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

unity_2.o: unity_2.c $(HEA_DIR)/gmsh.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

unity_3.o: unity_3.c $(HEA_DIR)/util.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

unity_4.o: unity_4.c $(HEA_DIR)/mesh.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

unity_ell.o: unity_ell.c $(HEA_DIR)/ell.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

$(OBJ_DIR)/myio.o: ../src/myio.c $(HEA_DIR)/myio.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -I/usr/lib/openmpi/include -std=gnu11

$(OBJ_DIR)/list.o: ../src/list.c $(HEA_DIR)/list.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

$(OBJ_DIR)/gmsh.o: ../src/gmsh.c $(HEA_DIR)/gmsh.h $(HEA_DIR)/myio.h $(HEA_DIR)/list.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

$(OBJ_DIR)/util.o: ../src/util.c $(HEA_DIR)/util.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

$(OBJ_DIR)/mesh.o: ../src/mesh.c $(HEA_DIR)/mesh.h $(HEA_DIR)/util.h $(HEA_DIR)/myio.h $(HEA_DIR)/list.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

$(OBJ_DIR)/ell.o: ../src/ell.c $(HEA_DIR)/ell.h
	$(MPICC) -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11

clean:	    
	rm -f $(OBJ_DIR)/* *.o unity_1 unity_2 unity_3 unity_4 unity_ell

run:
	$(MPIRUN) -np 1 ./unity_1
	$(MPIRUN) -np 1 ./unity_2
	$(MPIRUN) -np 2 ./unity_2
	$(MPIRUN) -np 3 ./unity_2
	$(MPIRUN) -np 1 ./unity_3
	$(MPIRUN) -np 1 ./unity_ell
