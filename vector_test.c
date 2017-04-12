#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <mpi.h>
#include "mpi_vector.h"
#include "vector.h"
#include "vector_macros.h"

/* Initialize the vectors with the numbers 0.0, 1.0, 2.0, … */
SG_DEFINE_VECTOR_MAP_1(static, init_f, v, { *v = index; })

#define N 42

int main(void)
{
    int np;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    {
        struct SgMpiVectorDriver mdrv =
            sg_mpi_vector_driver_new(MPI_COMM_WORLD, N);
        struct SgVectorDriver drv = sg_mpi_vector_driver_get(&mdrv);

        SgVector *v = sg_vector_new(drv);

        sg_vector_operate(drv, NULL, 0, init_f, NULL, 0, &v, 1);

        assert(sg_vector_sum(drv, v) == (N * np) * (N * np - 1) / 2.0);

        sg_vector_del(drv, v);
    }

    MPI_Finalize();
}
