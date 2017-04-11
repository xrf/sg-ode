#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "proto_vector.h"

static void vector_sum_f(void *f_ctx,
                         Accum accum,
                         Accum val,
                         double **data,
                         size_t num_elems)
{
    const double *restrict v;
    size_t i;
    double acc = *(const double *)accum;
    (void)f_ctx;
    acc += *(const double *)val;
    if (num_elems) {
        v = data[0];
        for (i = 0; i < num_elems; ++i) {
            acc += v[i];
        }
    }
    *(double *)accum = acc;
}

double vector_sum(void *self, const struct ProtoVector *self_vtable, Vector v)
{
    double accum;
    (*self_vtable->fold_map)(self, &accum, -1, vector_sum_f, NULL, &v, 1);
    return accum;
}



static void *check_oom(void *ptr) {
    if (!ptr) {
        fprintf(stderr, "Out of memory\n");
        fflush(stderr);
        abort();
    }
    return ptr;
}

// TODO: move the MPI stuff into its own header to avoid polluting with MPI stuff

struct MpiProtoVector {
    MPI_Comm comm;
    size_t len;
};

extern struct ProtoVector MPI_PROTO_VECTOR;

static Vector mpi_vector_create(void *self)
{
    struct MpiProtoVector *const inner = (struct MpiProtoVector *)self;
    return (Vector)check_oom(malloc(inner->len * sizeof(double)));
}

static void mpi_vector_destroy(void *self, Vector vector)
{
    (void)self;
    free(vector);
}

static _Thread_local void *fptr_ctx;
static _Thread_local FoldMapFn fptr;

static void op_wrapper(void *invec, void *inoutvec,
                       int *len, MPI_Datatype *datatype)
{
    (void)len;
    (void)datatype;
    (*fptr)(fptr_ctx, inoutvec, invec, NULL, 0);
}

static void mpi_vector_fold_map(void *self,
                                Accum accum,
                                AccumType accum_type,
                                FoldMapFn f,
                                void *f_ctx,
                                Vector *vectors,
                                size_t num_vectors)
{
    struct MpiProtoVector *const inner = (struct MpiProtoVector *)self;
    size_t i, accum_size;
    int e, accum_len;
    Accum zero;
    MPI_Datatype datatype;
    MPI_Op op;
    double **data;

    if (accum_type >= 0) {
        accum_len = accum_type;
        accum_size = (size_t)accum_len;
        datatype = MPI_UNSIGNED_CHAR;
    } else {
        if (INT_MIN < -INT_MAX && accum_type == INT_MIN) {
            fprintf(stderr, "Integer overflow in accum_len\n");
            fflush(stderr);
            abort();
        }
        accum_len = -accum_type;
        accum_size = (size_t)accum_len * sizeof(double);
        datatype = MPI_DOUBLE;
    }

    zero = (Accum)check_oom(malloc(accum_size));
    memcpy(zero, accum, accum_size);

    e = MPI_Op_create(&op_wrapper, 0, &op);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
        abort();
    }

    data = (double **)check_oom(malloc(num_vectors * sizeof(*data)));
    for (i = 0; i < num_vectors; ++i) {
        data[i] = (double *)vectors[i];
    }

    (*f)(f_ctx, accum, zero, data, inner->len);

    fptr = f;
    fptr_ctx = f_ctx;

    e = MPI_Allreduce(MPI_IN_PLACE, accum, accum_len,
                      datatype, op, inner->comm);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
        abort();
    }

    e = MPI_Op_free(&op);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
        abort();
    }

    free(data);
    free(zero);
}

struct ProtoVector MPI_PROTO_VECTOR = {
    mpi_vector_create,
    mpi_vector_destroy,
    mpi_vector_fold_map
};
