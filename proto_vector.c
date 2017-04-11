#include <string.h>
#include "proto_vector.h"

static void vector_sum_f(void *f_ctx,
                         Accum accum,
                         Accum val,
                         size_t offset,
                         double *restrict *data,
                         size_t num_elems)
{
    double s = *(const double *)accum + *(const double *)val;
    (void)f_ctx;
    (void)offset;
    if (num_elems) {
        size_t i;
        const double *const restrict v = data[0];
        for (i = 0; i < num_elems; ++i) {
            s += v[i];
        }
    }
    *(double *)accum = s;
}

double vector_sum(void *self,
                  const struct ProtoVectorVtable *vtable,
                  Vector v)
{
    double accum;
    (*vtable->fold_map)(self, &accum, -1, vector_sum_f, NULL, 0, &v, 1);
    return accum;
}
