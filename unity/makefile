
SRC_DIR= ../src
OBJ_DIR= ../obj
HEA_DIR= ../inc

all: unity_1

unity_1: unity_1.o $(OBJ_DIR)/myio.o  
	gcc -o unity_1 $^
	@echo "unity_1 generated..." 

unity_1.o: unity_1.c $(HEA_DIR)/myio.h
	gcc -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11
	@echo ">>> "$@

$(OBJ_DIR)/myio.o: ../src/myio.c $(HEA_DIR)/myio.h
	gcc -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11
	@echo ">>> "$@

clean:	    
	rm -f $(OBJ_DIR)/* unity_1 unity_1.o

run:
	./unity_1
