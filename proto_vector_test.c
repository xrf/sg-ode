#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <mpi.h>
#include "mpi_proto_vector.h"
#include "proto_vector.h"
#include "proto_vector_macros.h"

/* Initialize the vectors with the numbers 0.0, 1.0, 2.0, … */
DEF_MAP_FN_1(static, init_f, v, { *v = index; })

#define N 42

int main(void)
{
    static const struct ProtoVectorVtable *const vt = &MPI_PROTO_VECTOR_VTABLE;
    int rank, np;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    {
        struct MpiProtoVector pv = {MPI_COMM_WORLD, N * (unsigned)rank, N};
        Vector v = (*vt->create)(&pv);

        (*vt->fold_map)(&pv, NULL, 0, init_f, NULL, 0, &v, 1);

        assert(vector_sum(&pv, vt, v) ==
               (N * np) * (N * np - 1) / 2.0);

        (*vt->destroy)(&pv, v);
    }

    MPI_Finalize();
}
