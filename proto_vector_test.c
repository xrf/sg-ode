#include <stdio.h>
#include <mpi.h>
#include "proto_vector.h"

extern struct ProtoVector MPI_PROTO_VECTOR;

struct MpiProtoVector {
    MPI_Comm comm;
    size_t len;
};

#define N 42

int main(void)
{
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    struct MpiProtoVector pv = {MPI_COMM_WORLD, N};

    Vector v = (*MPI_PROTO_VECTOR.create)(&pv);

    // how to initialize vector without access to the indices ???
//    (*MPI_PROTO_VECTOR.create)(&pv);

    for (int i = 0; i < N; ++i) {
        ((double *)v)[i] = N * rank + i;
    }

    printf("%f\n", vector_sum(&pv, &MPI_PROTO_VECTOR, v));

    (*MPI_PROTO_VECTOR.destroy)(&pv, v);

    MPI_Finalize();
}
