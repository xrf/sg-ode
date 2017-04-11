#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "proto_vector.h"
#include "simple_proto_vector.h"
#include "utils.h"
#include "mpi_proto_vector.h"

static Vector mpi_vector_create(void *self)
{
    struct MpiProtoVector *const inner = (struct MpiProtoVector *)self;
    return check_oom(malloc(inner->len * sizeof(double)));
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
    (*fptr)(fptr_ctx, inoutvec, invec, 0, NULL, 0);
}

static void mpi_vector_fold_map(void *self,
                                Accum accum,
                                AccumType accum_type,
                                FoldMapFn f,
                                void *f_ctx,
                                size_t offset,
                                Vector *vectors,
                                size_t num_vectors)
{
    struct MpiProtoVector *const inner = (struct MpiProtoVector *)self;
    int e, accum_len;
    MPI_Datatype datatype;
    MPI_Op op;

    if (accum_type >= 0) {
        accum_len = accum_type;
        datatype = MPI_UNSIGNED_CHAR;
    } else {
        accum_len = -accum_type;
        datatype = MPI_DOUBLE;
    }

    e = MPI_Op_create(&op_wrapper, 0, &op);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
        abort();
    }

    (*SIMPLE_PROTO_VECTOR_VTABLE.fold_map)(&inner->len,
                                           accum,
                                           accum_type,
                                           f,
                                           f_ctx,
                                           offset + inner->offset,
                                           vectors,
                                           num_vectors);

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
}

struct ProtoVectorVtable MPI_PROTO_VECTOR_VTABLE = {
    mpi_vector_create,
    mpi_vector_destroy,
    mpi_vector_fold_map
};
