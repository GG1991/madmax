#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{

    int       rank;
    int       size;
    MPI_Comm  world;

    world = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(world, &size);
    MPI_Comm_rank(world, &rank);

    printf("rank %d from prog_2 size %d\n",rank, size);

    MPI_Finalize();

    return 0;
}
