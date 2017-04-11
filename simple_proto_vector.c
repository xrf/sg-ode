#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "proto_vector.h"
#include "utils.h"
#include "simple_proto_vector.h"

static Vector simple_vector_create(void *self)
{
    const size_t len = *(const size_t *)self;
    return check_oom(malloc(len * sizeof(double)));
}

static void simple_vector_destroy(void *self, Vector vector)
{
    (void)self;
    free(vector);
}

static void simple_vector_fold_map(void *self,
                                   Accum accum,
                                   AccumType accum_type,
                                   FoldMapFn f,
                                   void *f_ctx,
                                   size_t offset,
                                   Vector *vectors,
                                   size_t num_vectors)
{
    const size_t len = *(const size_t *)self;
    size_t i, accum_size;
    int accum_len;
    Accum zero;
    double **data;

    if (accum_type >= 0) {
        accum_len = accum_type;
        accum_size = (size_t)accum_len;
    } else {
        if (INT_MIN < -INT_MAX &&
            accum_type == INT_MIN &&
            /* convoluted comparison to avoid
               -Wtautological-constant-out-of-range-compare */
            INT_MAX > (size_t)(-1) / sizeof(double) &&
            (unsigned)(-accum_len) > (size_t)(-1) / sizeof(double)) {
            fprintf(stderr, "Integer overflow in accum_len\n");
            fflush(stderr);
            abort();
        }
        accum_len = -accum_type;
        accum_size = (size_t)accum_len * sizeof(double);
    }

    zero = check_oom(malloc(accum_size));
    memcpy(zero, accum, accum_size);

    data = (double **)check_oom(malloc(num_vectors * sizeof(*data)));
    for (i = 0; i < num_vectors; ++i) {
        data[i] = (double *)vectors[i];
    }

    (*f)(f_ctx, accum, zero, offset, data, len);

    free(data);
    free(zero);
}

struct ProtoVectorVtable SIMPLE_PROTO_VECTOR_VTABLE = {
    simple_vector_create,
    simple_vector_destroy,
    simple_vector_fold_map
};
