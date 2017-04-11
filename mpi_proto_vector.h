#ifndef G_S34QMP94V6RNSEMXM74ZIPLM8ID4U
#define G_S34QMP94V6RNSEMXM74ZIPLM8ID4U
#include "proto_vector.h"
#ifdef __cplusplus
extern "C" {
#endif

struct MpiProtoVector {
    MPI_Comm comm;
    size_t offset;
    size_t len;
};

extern struct ProtoVectorVtable MPI_PROTO_VECTOR_VTABLE;

#ifdef __cplusplus
}
#endif
#endif
