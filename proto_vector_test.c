#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <mpi.h>
#include "mpi_proto_vector.h"
#include "proto_vector.h"

/* we define a separate inner function because Clang and GCC seem to ignore
   'restrict' on the array pointer parameters and non-parameter variables */

#define DEF_MAP_FN_1(name, var1, block)                 \
    static void name##_inner(size_t _offset,            \
                             double *restrict _v1,      \
                             size_t _num_elems)         \
    {                                                   \
        size_t _i;                                      \
        for (_i = 0; _i < _num_elems; ++_i) {           \
            size_t index = _offset + _i;                \
            double *const var1 = _v1 + _i;              \
            (void)index;                                \
            { block }                                   \
        }                                               \
    }                                                   \
    static void name(void *f_ctx,                       \
                     Accum accum,                       \
                     Accum val,                         \
                     size_t offset,                     \
                     double *restrict *data,            \
                     size_t num_elems)                  \
    {                                                   \
        (void)f_ctx;                                    \
        (void)accum;                                    \
        (void)val;                                      \
        if (num_elems) {                                \
            name##_inner(offset, data[0], num_elems);   \
        }                                               \
    }

/* Initialize the vectors with the numbers 0.0, 1.0, 2.0, … */
DEF_MAP_FN_1(init_f, v, { *v = index; })

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
