# sputnik

This is a high performance algorithm for solving multi-scale problems.

To use the PETSc profiling do:

```bash
 mpirun -np <np_mac> macro/macro input.mac -log_view ascii:log_summary_mac.dat :
        -np <np_mic> micro/micro input.mic -log_view ascii:log_summary_mic.dat 
```

To plot PETSc trace do:

```bash
 mpirun -np <np_mac> macro/macro input.mac -log_trace ascii:trace_mac.dat :
        -np <np_mic> micro/micro input.mic -log_trace ascii:trace_mic.dat 
```

It will gives you one file per processor

# ParMETIS

```bash
cmake -DGKLIB_PATH=$HOME/libs/parmetis-4.0.3/metis/GKlib 
      -DMETIS_PATH=$HOME/libs/parmetis-4.0.3/metis 
      -DCMAKE_C_COMPILER=$HOME/libs/openmpi-install/bin/mpicc .
```
