#include <limits.h>
#include <stdio.h>
#include <stdint.h> /* for uint64_t */
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "vector.h"
#include "utils.h"
#include "mpi_vector.h"

struct SgMpiVectorDriver sg_mpi_vector_driver_new(MPI_Comm comm, size_t len)
{
    struct SgMpiVectorDriver d;
    uint64_t offset = 0, total_len;
    int e, rank, size;

    e = MPI_Comm_rank(comm, &rank);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
    }

    e = MPI_Comm_size(comm, &size);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
    }

    if (rank != 0) {
        e = MPI_Recv(&offset, 1, MPI_UINT64_T, rank - 1,
                     0, comm, MPI_STATUS_IGNORE);
        if (e) {
            MPI_Abort(MPI_COMM_WORLD, e);
        }
    }

    total_len = offset + len;

    if (rank != size - 1) {
        e = MPI_Send(&total_len, 1, MPI_UINT64_T, rank + 1, 0, comm);
        if (e) {
            MPI_Abort(MPI_COMM_WORLD, e);
        }
    }

    e = MPI_Bcast(&total_len, 1, MPI_UINT64_T, size - 1, comm);
    if (e) {
        MPI_Abort(MPI_COMM_WORLD, e);
    }

    d.base.len = (size_t)total_len;
    d.comm = comm;
    d.offset = (size_t)offset;
    d.driver = sg_basic_vector_driver_new(len);
    return d;
}

struct SgVectorDriver sg_mpi_vector_driver_get(struct SgMpiVectorDriver *d)
{
    struct SgVectorDriver drv = {&d->base, &SG_MPI_VECTOR_DRIVER_VT};
    return drv;
}

static SgVector *vector_try_new(struct SgVectorDriverBase *self)
{
    struct SgMpiVectorDriver *const inner = (struct SgMpiVectorDriver *)self;
    struct SgVectorDriver drv = sg_basic_vector_driver_get(&inner->driver);
    return sg_vector_try_new(drv);
}

static void vector_del(struct SgVectorDriverBase *self, SgVector *vector)
{
    struct SgMpiVectorDriver *const inner = (struct SgMpiVectorDriver *)self;
    struct SgVectorDriver drv = sg_basic_vector_driver_get(&inner->driver);
    sg_vector_del(drv, vector);
}

static _Thread_local void *fptr_ctx;
static _Thread_local SgVectorOperation *fptr;

static void op_wrapper(void *invec, void *inoutvec,
                       int *len, MPI_Datatype *datatype)
{
    (void)len;
    (void)datatype;
    (*fptr)(fptr_ctx, inoutvec, invec, 0, NULL, 0);
}

static void vector_operate(struct SgVectorDriverBase *self,
                           SgVectorAccum *accum,
                           SgVectorAccumType accum_type,
                           SgVectorOperation *f,
                           void *f_ctx,
                           size_t offset,
                           SgVector **vectors,
                           size_t num_vectors)
{
    struct SgMpiVectorDriver *const inner = (struct SgMpiVectorDriver *)self;
    struct SgVectorDriver drv = sg_basic_vector_driver_get(&inner->driver);
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
