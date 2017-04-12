#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "vector.h"

SgVector *sg_vector_try_new(struct SgVectorDriver drv)
{
    return (*drv.vtable->try_new)(drv.data);
}

SgVector *sg_vector_new(struct SgVectorDriver drv)
{
    SgVector *v = sg_vector_try_new(drv);
    check_oom(v);
    return v;
}

void sg_vector_del(struct SgVectorDriver drv, SgVector *vector)
{
    (*drv.vtable->del)(drv.data, vector);
}

size_t sg_vector_len(struct SgVectorDriver drv)
{
    return drv.data->len;
}

void sg_vector_operate(struct SgVectorDriver drv,
                       SgVectorAccum *accum,
                       SgVectorAccumType accum_type,
                       SgVectorOperation *f,
                       void *f_ctx,
                       size_t offset,
                       SgVector **vectors,
                       size_t num_vectors)
{
    (*drv.vtable->operate)(drv.data, accum, accum_type, f, f_ctx,
                           offset, vectors, num_vectors);
}

void sg_vector_fill(struct SgVectorDriver drv,
                    const SgVector *vector, double value)
{
    SgVector *v = (SgVector *)vector;
    sg_vector_operate(drv, NULL, 0,
                      &sg_vector_fill_operation, &value, 0, &v, 1);
}

void sg_vector_fill_operation(void *f_ctx,
                              SgVectorAccum *accum,
                              const SgVectorAccum *val,
                              size_t offset,
                              double **data,
                              size_t num_elems)
{
    double c = *(const double *)f_ctx;
    (void)accum;
    (void)val;
    (void)offset;
    if (num_elems) {
        size_t i;
        double *v = data[0];
        if (c == 0.0) {
            for (i = 0; i < num_elems; ++i) {
                v[i] = 0.0;
            }
        } else {
            for (i = 0; i < num_elems; ++i) {
                v[i] = c;
            }
        }
    }
}

double sg_vector_sum(struct SgVectorDriver drv, const SgVector *vector)
{
    double accum = 0.0;
    SgVector *v = (SgVector *)vector;
    sg_vector_operate(drv, &accum, -1,
                      &sg_vector_sum_operation, NULL, 0, &v, 1);
    return accum;
}

void sg_vector_sum_operation(void *f_ctx,
                             SgVectorAccum *accum,
                             const SgVectorAccum *val,
                             size_t offset,
                             double **data,
                             size_t num_elems)
{
    double s = *(const double *)accum + *(const double *)val;
    (void)f_ctx;
    (void)offset;
    if (num_elems) {
        size_t i;
        const double *v = data[0];
        for (i = 0; i < num_elems; ++i) {
            s += v[i];
        }
    }
    *(double *)accum = s;
}

struct SgBasicVectorDriver sg_basic_vector_driver_new(size_t len)
{
    struct SgBasicVectorDriver d = {{len}};
    return d;
}

struct SgVectorDriver sg_basic_vector_driver_get(struct SgBasicVectorDriver *d)
{
    struct SgVectorDriver drv = {&d->base, &SG_BASIC_VECTOR_DRIVER_VT};
    return drv;
}

static SgVector *basic_vector_try_new(struct SgVectorDriverBase *self)
{
    const size_t len = *(const size_t *)self;
    return check_oom(malloc(len * sizeof(double)));
}

static void basic_vector_del(struct SgVectorDriverBase *self, SgVector *vector)
{
    (void)self;
    free(vector);
}

static void basic_vector_operate(struct SgVectorDriverBase *self,
                                 SgVectorAccum *accum,
                                 SgVectorAccumType accum_type,
                                 SgVectorOperation f,
                                 void *f_ctx,
                                 size_t offset,
                                 SgVector **vectors,
                                 size_t num_vectors)
{
    const size_t len = self->len;
    size_t i, accum_size;
    int accum_len;
    SgVectorAccum *zero;
    double **data;

    if (accum_type >= 0) {
        accum_len = accum_type;
        accum_size = (size_t)accum_len;
    } else {
        /* trying very hard to avoid warnings such as
           -Wtautological-constant-out-of-range-compare */
        static const size_t max = (size_t)(-1) / sizeof(double);
        if (INT_MIN < -INT_MAX &&
            accum_type == INT_MIN &&
            (unsigned)(-accum_type) > max) {
            fprintf(stderr, "Integer overflow in accum_len\n");
            fflush(stderr);
            abort();
        }
        accum_len = -accum_type;
        accum_size = (size_t)accum_len * sizeof(double);
    }

    zero = check_oom(malloc(accum_size));
    if (accum) {
        memcpy(zero, accum, accum_size);
    }

    data = (double **)check_oom(malloc(num_vectors * sizeof(*data)));
    for (i = 0; i < num_vectors; ++i) {
        data[i] = (double *)vectors[i];
    }

    (*f)(f_ctx, accum, zero, offset, data, len);

    free(data);
    free(zero);
}

struct SgVectorDriverVt SG_BASIC_VECTOR_DRIVER_VT = {
    &basic_vector_try_new,
    &basic_vector_del,
    &basic_vector_operate
};
