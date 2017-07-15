# sputnik

This is a high performance algorithm for solving multi-scale problems.

To use the PETSc profiling do:

```bash
 mpirun -np <np_mac> macro/macro input.mac -log_view ascii:log_summary_mac.dat :
        -np <np_mic> micro/micro input.mic -log_view ascii:log_summary_mic.dat 
```
