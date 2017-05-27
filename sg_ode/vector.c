#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "vector_macros.h"
#include "vector.h"
#include "restrict_begin.h"

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
    return (*drv.vtable->len)(drv.data);
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
                    double value,
                    SgVector *vector)
{
    sg_vector_operate(drv, NULL, 0,
                      &sg_vector_fill_operation, &value,
                      0, &vector, 1);
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
            /* compilers should be able to optimize this to a simple memset
               (can't use a real memset because C standard doesn't guarantee
               that would work) */
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

void sg_vector_copy(struct SgVectorDriver drv,
                    const SgVector *src,
                    SgVector *dest)
{
    SgVector *v[2];
    v[0] = (SgVector *)src;
    v[1] = dest;
    sg_vector_operate(drv, NULL, 0, &sg_vector_copy_operation, NULL,
                      0, v, sizeof(v) / sizeof(*v));
}

SG_DEFINE_VECTOR_MAP_2(extern, sg_vector_copy_operation, src, dest, {
        *dest = *src;
    })

void sg_vector_neg_assign(struct SgVectorDriver drv, SgVector *z)
{
    sg_vector_operate(drv, NULL, 0,
                      &sg_vector_neg_assign_operation, NULL,
                      0, &z, 1);
}

SG_DEFINE_VECTOR_MAP_1(extern, sg_vector_neg_assign_operation, z, {
        *z *= -1.0;
    })

void sg_vector_neg(struct SgVectorDriver drv, const SgVector *x, SgVector *z)
{
    SgVector *v[2];
    v[0] = (SgVector *)x;
    v[1] = z;
    sg_vector_operate(drv, NULL, 0,
                      &sg_vector_neg_assign_operation, NULL,
                      0, v, sizeof(v) / sizeof(*v));
}

SG_DEFINE_VECTOR_MAP_2(extern, sg_vector_neg_operation, z, x, {
        *z = -*x;
    })

void sg_vector_scale_assign(struct SgVectorDriver drv,
                            double alpha,
                            SgVector *z)
{
    if (alpha == 0.0) {
        sg_vector_fill(drv, 0.0, z);
    } else if (alpha == 1.0) {
        /* do nothing */
    } else if (alpha == -1.0) {
        sg_vector_neg_assign(drv, z);
    } else {
        sg_vector_operate(drv, NULL, 0,
                          &sg_vector_scale_assign_operation, &alpha,
                          0, &z, 1);
    }
}

SG_DEFINE_VECTOR_MAP_1(extern, sg_vector_scale_assign_operation, z, {
        const double c = *(const double *)ctx;
        *z *= c;
    })

void sg_vector_scale(struct SgVectorDriver drv,
                     double alpha,
                     const SgVector *x,
                     SgVector *z)
{
    if (x == z) {
        sg_vector_scale_assign(drv, alpha, z);
    } else if (alpha == 0.0) {
        sg_vector_fill(drv, 0.0, z);
    } else if (alpha == 1.0) {
        /* do nothing */
    } else if (alpha == -1.0) {
        sg_vector_neg(drv, x, z);
    } else {
        SgVector *v[2];
        v[0] = (SgVector *)x;
        v[1] = z;
        sg_vector_operate(drv, NULL, 0,
                          &sg_vector_scale_operation, &alpha,
                          0, v, sizeof(v) / sizeof(*v));
    }
}

SG_DEFINE_VECTOR_MAP_2(extern, sg_vector_scale_operation, x, z, {
        const double c = *(const double *)ctx;
        *z = c * *x;
    })

void sg_vector_linear_assign(struct SgVectorDriver drv,
                             double alpha,
                             double beta,
                             const SgVector *y,
                             SgVector *z)
{
    if (y == z || beta == 0.0) {
        sg_vector_scale_assign(drv, alpha + beta, z);
    } else if (alpha == 0.0) {
        sg_vector_scale(drv, beta, y, z);
    } else {
        double c[2];
        SgVector *v[2];
        c[0] = alpha;
        c[1] = beta;
        v[0] = (SgVector *)y;
        v[1] = z;
        sg_vector_operate(drv, NULL, 0,
                          &sg_vector_linear_assign_operation, &c,
                          0, v, sizeof(v) / sizeof(*v));
    }
}

SG_DEFINE_VECTOR_MAP_2(extern, sg_vector_linear_assign_operation, y, z, {
        const double *const c = (const double *)ctx;
        *z = c[0] * *z + c[1] * *y;
    })

void sg_vector_linear(struct SgVectorDriver drv,
                      double alpha,
                      const SgVector *x,
                      double beta,
                      const SgVector *y,
                      SgVector *z)
{
    if (x == z) {
        sg_vector_linear_assign(drv, alpha, beta, y, z);
    } else if (y == z) {
        sg_vector_linear_assign(drv, beta, alpha, x, z);
    } else if (alpha == 0.0) {
        sg_vector_scale(drv, beta, y, z);
    } else if (beta == 0.0) {
        sg_vector_scale(drv, alpha, x, z);
    } else {
        double c[2];
        SgVector *v[3];
        c[0] = alpha;
        c[1] = beta;
        v[0] = (SgVector *)x;
        v[1] = (SgVector *)y;
        v[2] = z;
        sg_vector_operate(drv, NULL, 0,
                          &sg_vector_linear_operation, &c,
                          0, v, sizeof(v) / sizeof(*v));
    }
}

SG_DEFINE_VECTOR_MAP_3(extern, sg_vector_linear_operation, x, y, z, {
        const double *const c = (const double *)ctx;
        *z = c[0] * *x + c[1] * *y;
    })

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

static size_t basic_vector_len(const SgVectorDriverBase *self)
{
    return *(const size_t *)self;
}

static SgVector *basic_vector_try_new(const SgVectorDriverBase *self)
{
    return malloc(basic_vector_len(self) * sizeof(double));
}

static void basic_vector_del(const SgVectorDriverBase *self, SgVector *vector)
{
    (void)self;
    free(vector);
}

static void basic_vector_operate(const SgVectorDriverBase *self,
                                 SgVectorAccum *accum,
                                 SgVectorAccumType accum_type,
                                 SgVectorOperation f,
                                 void *f_ctx,
                                 size_t offset,
                                 SgVector **vectors,
                                 size_t num_vectors)
{
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

    (*f)(f_ctx, accum, zero, offset, data, basic_vector_len(self));

    free(data);
    free(zero);
}

struct SgVectorDriverVt SG_BASIC_VECTOR_DRIVER_VT = {
    &basic_vector_len,
    &basic_vector_try_new,
    &basic_vector_del,
    &basic_vector_operate
};
