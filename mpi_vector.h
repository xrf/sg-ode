#ifndef G_S34QMP94V6RNSEMXM74ZIPLM8ID4U
#define G_S34QMP94V6RNSEMXM74ZIPLM8ID4U
#include <stddef.h>
#include <mpi.h>
#include "vector.h"
#ifdef __cplusplus
extern "C" {
#endif

struct SgMpiVectorDriver {
    /** Contains the total length */
    struct SgVectorDriverBase base;

    MPI_Comm comm;

    /** Offset of this subvector relative to global vector */
    size_t offset;

    /** Contains the local length */
    struct SgBasicVectorDriver driver;
};

/** Create a driver that distributes a vector over multiple nodes.  This is a
    collective call, therefore every node within the communicator is expected
    to participate.  The `len` parameter should be the length of the current
    node’s subvector, not the total length. */
struct SgMpiVectorDriver sg_mpi_vector_driver_new(MPI_Comm comm, size_t len);

struct SgVectorDriver sg_mpi_vector_driver_get(struct SgMpiVectorDriver *);

extern struct SgVectorDriverVt SG_MPI_VECTOR_DRIVER_VT;

#ifdef __cplusplus
}
#endif
#endif
