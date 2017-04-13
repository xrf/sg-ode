#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "vector.h"
#include "utils.h"
#include "mpi_vector.h"
#include "thread_local_begin.h"

static void encode_size(size_t n, unsigned char *buf, size_t len)
{
    size_t i;
    for (i = 0; i < len; ++i) {
        buf[i] = (unsigned char)((n >> (i * CHAR_BIT)) & ((1 << CHAR_BIT) - 1));
    }
    for (i = len; i < sizeof(n); ++i) {
        if ((n >> (i * CHAR_BIT)) & ((1 << CHAR_BIT) - 1)) {
            fprintf(stderr,
                    "mpi_vector.c:encode_size: size_t can't fit in buffer\n");
            fflush(stderr);
            abort();
        }
    }
}

static size_t decode_size(const unsigned char *buf, size_t len)
{
    size_t i, n = 0;
    for (i = 0; i < sizeof(n); ++i) {
        n |= (size_t)buf[i] << (i * CHAR_BIT);
    }
    for (i = sizeof(n); i < len; ++i) {
        if (buf[i]) {
            fprintf(stderr,
                    "mpi_vector.c:decode_size: integer too large for size_t\n");
            fflush(stderr);
            abort();
        }
    }
    return n;
}

struct SgMpiVectorDriver sg_mpi_vector_driver_new(MPI_Comm comm, size_t len)
{
    struct SgMpiVectorDriver d;
    int e, rank, size;
    /* use a fixed-width buffer to ensure nodes with differing size_t can
       still communicate */
    unsigned char buf[8];

    e = MPI_Comm_rank(comm, &rank);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
    }

    e = MPI_Comm_size(comm, &size);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
    }

    if (rank == 0) {
        d.offset = 0;
    } else {
        e = MPI_Recv(buf, sizeof(buf), MPI_UNSIGNED_CHAR, rank - 1,
                     0, comm, MPI_STATUS_IGNORE);
        d.offset = decode_size(buf, sizeof(buf));
        if (e) {
            MPI_Abort(MPI_COMM_WORLD, e);
        }
    }

    if (d.offset > (size_t)(-1) - len) {
        fprintf(stderr, "sg_mpi_vector_driver_new: size_t overflow\n");
        fflush(stderr);
        abort();
    }
    encode_size(d.offset + len, buf, sizeof(buf));

    if (rank != size - 1) {
        e = MPI_Send(buf, sizeof(buf), MPI_UNSIGNED_CHAR, rank + 1, 0, comm);
        if (e) {
            MPI_Abort(MPI_COMM_WORLD, e);
        }
    }

    e = MPI_Bcast(buf, sizeof(buf), MPI_UNSIGNED_CHAR, size - 1, comm);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
    }

    d.base = decode_size(buf, sizeof(buf));
    d.comm = comm;
    d.local_len = len;
    return d;
}

struct SgVectorDriver sg_mpi_vector_driver_get(struct SgMpiVectorDriver *d)
{
    struct SgVectorDriver drv;
    drv.data = &d->base;
    drv.vtable = &SG_MPI_VECTOR_DRIVER_VT;
    return drv;
}

static SgVector *vector_try_new(SgVectorDriverBase *self)
{
    struct SgMpiVectorDriver *const inner = (struct SgMpiVectorDriver *)self;
    struct SgVectorDriver drv = sg_basic_vector_driver(&inner->local_len);
    return sg_vector_try_new(drv);
}

static void vector_del(SgVectorDriverBase *self, SgVector *vector)
{
    struct SgMpiVectorDriver *const inner = (struct SgMpiVectorDriver *)self;
    struct SgVectorDriver drv = sg_basic_vector_driver(&inner->local_len);
    sg_vector_del(drv, vector);
}

static thread_local void *fptr_ctx;
static thread_local SgVectorOperation *fptr;

static void op_wrapper(void *invec, void *inoutvec,
                       int *len, MPI_Datatype *datatype)
{
    (void)len;
    (void)datatype;
    (*fptr)(fptr_ctx, inoutvec, invec, 0, NULL, 0);
}

static void vector_operate(SgVectorDriverBase *self,
                           SgVectorAccum *accum,
                           SgVectorAccumType accum_type,
                           SgVectorOperation *f,
                           void *f_ctx,
                           size_t offset,
                           SgVector **vectors,
                           size_t num_vectors)
{
    struct SgMpiVectorDriver *const inner = (struct SgMpiVectorDriver *)self;
    struct SgVectorDriver drv = sg_basic_vector_driver(&inner->local_len);
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
    }

    sg_vector_operate(drv,
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
    }

    e = MPI_Op_free(&op);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
    }
}

struct SgVectorDriverVt SG_MPI_VECTOR_DRIVER_VT = {
    vector_try_new,
    vector_del,
    vector_operate
};
