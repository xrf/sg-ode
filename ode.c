/* Original: http://www.netlib.org/ode/ode.f */
/* Converted through f2c and then manually refactored. */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "vector.h"
#include "vector_macros.h"
#include "ode.h"

typedef void (*fn_type)(void *f_ctx,
                        double t,
                        const SgVector *restrict y,
                        SgVector *restrict yp);

static void vector_div_normsq_operation_inner(double *restrict sum,
                                              const double *restrict numer,
                                              const double *restrict denom,
                                              size_t num_elems)
{
    size_t i;
    for (i = 0; i < num_elems; ++i) {
        *sum += pow(numer[i] / denom[i], 2.0);
    }
}

static void vector_div_normsq_operation(void *f_ctx,
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
        vector_div_normsq_operation_inner(&s, data[0], data[1], num_elems);
    }
    *(double *)accum = s;
}

double vector_div_normsq(struct SgVectorDriver drv,
                         const SgVector *restrict numer,
                         const SgVector *restrict denom)
{
    double accum = 0.0;
    SgVector *v[] = {(SgVector *)numer, (SgVector *)denom};
    sg_vector_operate(drv, &accum, -1, &vector_div_normsq_operation, NULL,
                      0, v, sizeof(v) / sizeof(*v));
    return accum;
}

SG_DEFINE_VECTOR_MAP_4(static, vector_taup_operation, y, phi14, phi15, p, {
        /* 𝛕 = h 𝐩 - 𝛗[14]
           𝐩 ← 𝐲 + 𝛕
           𝛗[15] ← (𝐩 - 𝐲) - 𝛕 */
        const double h = *(const double *)ctx;
        const double tau = h * *p - *phi14;
        *p = *y + tau;
        *phi15 = (*p - *y) - tau;
    })

void vector_taup(struct SgVectorDriver drv,
                 double h,
                 const SgVector *restrict y,
                 const SgVector *restrict phi14,
                 const SgVector *restrict phi15,
                 SgVector *restrict p)
{
    SgVector *v[] = {(SgVector *)y, (SgVector *)phi14, (SgVector *)phi15, p};
    sg_vector_operate(drv, NULL, 0, &vector_taup_operation, &h,
                      0, v, sizeof(v) / sizeof(*v));
}

typedef void VectorErkOperation(double *restrict erk,
                                double *restrict erkm1,
                                double *restrict erkm2,
                                const double *restrict wt,
                                const double *restrict yp,
                                const double *restrict phi0,
                                const double *restrict phikm1,
                                const double *restrict phikm2,
                                size_t num_elems);

static void vector_erk_operation_inner1(double *restrict erk,
                                        double *restrict erkm1,
                                        double *restrict erkm2,
                                        const double *restrict wt,
                                        const double *restrict yp,
                                        const double *restrict phi0,
                                        const double *restrict phikm1,
                                        const double *restrict phikm2,
                                        size_t num_elems)
{
    size_t i;
    (void)erkm1;
    (void)erkm2;
    (void)phikm1;
    (void)phikm2;
    for (i = 0; i < num_elems; ++i) {
        /* erk = ‖(𝐲′ - 𝛗[0]) / 𝛚‖² */
        const double iwt = 1.0 / wt[i];
        const double ypmphi = yp[i] - phi0[i];
        *erk += pow(ypmphi * iwt, 2.0);
    }
}

static void vector_erk_operation_inner2(double *restrict erk,
                                        double *restrict erkm1,
                                        double *restrict erkm2,
                                        const double *restrict wt,
                                        const double *restrict yp,
                                        const double *restrict phi0,
                                        const double *restrict phikm1,
                                        const double *restrict phikm2,
                                        size_t num_elems)
{
    size_t i;
    (void)erkm2;
    (void)phikm2;
    for (i = 0; i < num_elems; ++i) {
        /* erkm1 = ‖(𝛗[k - 1] + 𝐲′ - 𝛗[0]) / 𝛚‖²
           erk = ‖(𝐲′ - 𝛗[0]) / 𝛚‖² */
        const double iwt = 1.0 / wt[i];
        const double ypmphi = yp[i] - phi0[i];
        *erkm1 += pow((phikm1[i] + ypmphi) * iwt, 2.0);
        *erk += pow(ypmphi * iwt, 2.0);
    }
}

static void vector_erk_operation_inner3(double *restrict erk,
                                        double *restrict erkm1,
                                        double *restrict erkm2,
                                        const double *restrict wt,
                                        const double *restrict yp,
                                        const double *restrict phi0,
                                        const double *restrict phikm1,
                                        const double *restrict phikm2,
                                        size_t num_elems)
{
    size_t i;
    for (i = 0; i < num_elems; ++i) {
        /* erkm2 = ‖(𝛗[k - 2] + 𝐲′ - 𝛗[0]) / 𝛚‖²
           erkm1 = ‖(𝛗[k - 1] + 𝐲′ - 𝛗[0]) / 𝛚‖²
           erk = ‖(𝐲′ - 𝛗[0]) / 𝛚‖² */
        const double iwt = 1.0 / wt[i];
        const double ypmphi = yp[i] - phi0[i];
        *erkm2 += pow((phikm2[i] + ypmphi) * iwt, 2.0);
        *erkm1 += pow((phikm1[i] + ypmphi) * iwt, 2.0);
        *erk += pow(ypmphi * iwt, 2.0);
    }
}

static void vector_erk_operation(void *ctx,
                                 SgVectorAccum *vaccum,
                                 const SgVectorAccum *vval,
                                 size_t offset,
                                 double **data,
                                 size_t num_elems)
{
    VectorErkOperation *inner = *(VectorErkOperation *const *)ctx;
    double *accum = (double *)vaccum;
    const double *val = (double *)vval;
    double erk = accum[0] + val[0];
    double erkm1 = accum[1] + val[1];
    double erkm2 = accum[2] + val[2];
    (void)offset;
    if (num_elems) {
        (*inner)(&erk, &erkm1, &erkm2,
                 data[0], data[1], data[2], data[3], data[4], num_elems);
    }
    accum[0] = erk;
    accum[1] = erkm1;
    accum[2] = erkm2;
}

void vector_erk(struct SgVectorDriver drv,
                unsigned k,
                const SgVector *restrict wt,
                const SgVector *restrict yp,
                SgVector *const restrict *phi,
                double *restrict erk,
                double *restrict erkm1,
                double *restrict erkm2)
{
    double accum[3] = {0.0, 0.0, 0.0};
    SgVector *v[] = {
        (SgVector *)wt,
        (SgVector *)yp,
        (SgVector *)phi[0],
        NULL,
        NULL,
    };
    VectorErkOperation *ctx;
    if (k > 2) {
        ctx = &vector_erk_operation_inner3;
        v[3] = (SgVector *)phi[k - 1];
        v[4] = (SgVector *)phi[k - 2];
    } else if (k >= 2) {
        ctx = &vector_erk_operation_inner2;
        v[3] = (SgVector *)phi[k - 1];
    } else {
        ctx = &vector_erk_operation_inner1;
    }
    sg_vector_operate(drv, accum, -3,
                      &vector_erk_operation, &ctx,
                      0, v, sizeof(v) / sizeof(*v));
    *erk = accum[0];
    *erkm1 = accum[1];
    *erkm2 = accum[2];
}

SG_DEFINE_VECTOR_MAP_2(static, vector_ibeta_diff_operation, phiip1, phii, {
        const double ibetai = *(const double *)ctx;
        /* 𝛗[i] ← β[i]⁻¹ (𝛗[i] - 𝛗[i + 1]) */
        *phii = ibetai * (*phii - *phiip1);
    })

void vector_ibeta_diff(struct SgVectorDriver drv,
                       double ibetai,
                       const SgVector *restrict phiip1,
                       SgVector *restrict phii)
{
    SgVector *v[] = {(SgVector *)phiip1, phii};
    sg_vector_operate(drv, NULL, 0, &vector_ibeta_diff_operation, &ibetai,
                      0, v, sizeof(v) / sizeof(*v));
}

SG_DEFINE_VECTOR_MAP_6(
    static, vector_rhop_operation,
    yp, phi0, phi15, p, y, phi14, {
        /* 𝛒 = h g[k] (𝐲′ - 𝛗[0]) - 𝛗[15]
           𝐲 ← 𝐩 + 𝛒
           𝛗[14] ← (𝐲 - 𝐩) - 𝛒 */
        const double hgk = *(const double *)ctx;
        const double rho = hgk * (*yp - *phi0) - *phi15;
        *y = *p + rho;
        *phi14 = (*y - *p) - rho;
    })

void vector_rhop(struct SgVectorDriver drv,
                 double hgk,
                 const SgVector *restrict yp,
                 const SgVector *restrict phi0,
                 const SgVector *restrict phi15,
                 const SgVector *restrict p,
                 SgVector *restrict y,
                 SgVector *restrict phi14)
{
    SgVector *v[] = {
        (SgVector *)yp,
        (SgVector *)phi0,
        (SgVector *)phi15,
        (SgVector *)p,
        y,
        phi14
    };
    sg_vector_operate(drv, NULL, 0, &vector_rhop_operation, &hgk,
                      0, v, sizeof(v) / sizeof(*v));
}

SG_DEFINE_VECTOR_MAP_4(static, vector_eval_y_operation, p, yp, phi0, y, {
        /* 𝐲 ← 𝐩 + h g[k] (𝐲′ - 𝛗[0]) */
        const double hgk = *(const double *)ctx;
        *y = *p + hgk * (*yp - *phi0);
    })

void vector_eval_y(struct SgVectorDriver drv,
                   const SgVector *restrict p,
                   double hgk,
                   const SgVector *restrict yp,
                   const SgVector *restrict phi0,
                   SgVector *restrict y)
{
    SgVector *v[] = {(SgVector *)p, (SgVector *)yp, (SgVector *)phi0, y};
    sg_vector_operate(drv, NULL, 0, &vector_eval_y_operation, &hgk,
                      0, v, sizeof(v) / sizeof(*v));
}

SG_DEFINE_VECTOR_MAP_4(
    static, vector_upd_diffs_operation,
    yp, phi0, phik, phikp1, {
        /* 𝛗[k] ← 𝐲′ - 𝛗[0]
           𝛗[k + 1] ← 𝛗[k] - 𝛗[k + 1] */
        *phik = *yp - *phi0;
        *phikp1 = *phik - *phikp1;
    })

void vector_upd_diffs(struct SgVectorDriver drv,
                      const SgVector *restrict yp,
                      const SgVector *restrict phi0,
                      SgVector *restrict phik,
                      SgVector *restrict phikp1)
{
    SgVector *v[] = {(SgVector *)yp, (SgVector *)phi0, phik, phikp1};
    sg_vector_operate(drv, NULL, 0, &vector_upd_diffs_operation, NULL,
                      0, v, sizeof(v) / sizeof(*v));
}

SG_DEFINE_VECTOR_MAP_3(static, vector_intrp_yout_operation, phii, yout, ypout, {
        /* 𝐲° ← 𝐲° + g[i] 𝛗[i]
           𝐲°′ ← 𝐲°′ + ρ[i] 𝛗[i] */
        const double gi = ((const double *)ctx)[0];
        const double rhoi = ((const double *)ctx)[1];
        *yout += gi * *phii;
        *ypout += rhoi * *phii;
    })

void vector_intrp_yout(struct SgVectorDriver drv,
                       double gi,
                       double rhoi,
                       const SgVector *restrict phii,
                       SgVector *restrict yout,
                       SgVector *restrict ypout)
{
    double ctx[] = {gi, rhoi};
    SgVector *v[] = {(SgVector *)phii, yout, ypout};
    sg_vector_operate(drv, NULL, 0, &vector_intrp_yout_operation, &ctx,
                      0, v, sizeof(v) / sizeof(*v));
}

SG_DEFINE_VECTOR_MAP_2(static, vector_update_wt_operation, yy, wt, {
        /* 𝛚 ← releps |𝐘| + abseps */
        const double releps = ((const double *)ctx)[0];
        const double abseps = ((const double *)ctx)[1];
        *wt = releps * fabs(*yy) + abseps;
    })

void vector_update_wt(struct SgVectorDriver drv,
                      double releps,
                      const SgVector *restrict yy,
                      double abseps,
                      SgVector *restrict wt)
{
    double ctx[] = {releps, abseps};
    SgVector *v[] = {(SgVector *)yy, wt};
    sg_vector_operate(drv, NULL, 0, &vector_update_wt_operation, &ctx,
                      0, v, sizeof(v) / sizeof(*v));
}

/*
  Integrates a system of first order ordinary differential equations one step,
  normally from `x` to `x+h`, using a modified divided difference form of the
  Adams PECE formulas.  Local extrapolation is used to improve absolute
  stability and accuracy.  The code adjusts its order and step size to control
  the local error per unit step in a generalized sense.  Special devices are
  included to control roundoff error and to detect when the user is requesting
  too much accuracy.

  @param x
  Independent variable

  @param y
  Solution vector at `x` (length: `neqn`)

  @param yp
  Derivative of solution vector at `x` after successful step (length: `neqn`)

  @param h
  Appropriate step size for next step.  Normally determined by code

  @param eps
  Local error tolerance.  Must be variable (length: `1`)

  @param wt
  Vector of weights for error criterion (length: `neqn`)

  @param start
  Boolean variable set `true` for first step, `false` otherwise

  @param hold
  Step size used for last successful step

  @param k
  Appropriate order for next step (determined by code).
  Invariant: k >= 1 && k < 13

  @param kold
  Order used for last successful step

  @return
  Nonzero when no step can be taken, zero otherwise (`crash`).

  The arrays `phi`, `psi` are required for the interpolation function
  `intrp`.  The array `p` is internal to the code.

  # Input to `step`

  ## First call

  The user must provide storage in their driver program for all arrays in the
  call list, namely

      y[neqn], wt[neqn], phi[neqn * 16], p[neqn], yp[neqn], psi[12]

  The user must also declare the `start` Boolean variable and the function
  `f(f_ctx, x, y, yp)` to evaluate

      dy[i]/dx = yp[i] = f(x, y[0], y[1], ..., y[neqn - 1])

  and initialize only the following parameters:

    - `x`: Initial value of the independent variable
    - `y`: Vector of initial values of dependent variables
    - `h`: Nominal step size indicating direction of integration
           and maximum size of step.
    - `eps`: Local error tolerance per step.
    - `wt`: Vector of non-zero weights for error criterion
    - `start`: Set to `true`.

  `step` requires the L2 norm of the vector with components
  `local_error[l] / wt[l]` be less than `eps` for a successful step.  The
  array `wt` allows the user to specify an error test appropriate for their
  problem.  For example,

    - `wt[l] = 1.0` specifies absolute error,
    - `wt[l] = fabs(y[l])` error relative to the most recent value of the l-th
      component of the solution,
    - `wt[l] = fabs(yp[l])` error relative to the most recent value of the
      l-th component of the derivative,
    - `wt[l] = max(wt[l], fabs(y[l]))` error relative to the largest magnitude
      of l-th component obtained so far,
    - `wt[l] = fabs(y(l)) * relerr / eps + abserr / eps` specifies a mixed
      relative-absolute test where `relerr` is relative error, `abserr` is
      absolute error and `eps = max(relerr, abserr)`.

  ## Subsequent calls

  Function `step` is designed so that all information needed to continue the
  integration, including the step size `h` and the order `k`, is returned with
  each step.  With the exception of the step size, the error tolerance, and
  the weights, none of the parameters should be altered.  The array `wt` must
  be updated after each step to maintain relative error tests like those
  above.  Normally the integration is continued just beyond the desired
  endpoint and the solution interpolated there with `intrp`.  If it is
  impossible to integrate beyond the endpoint, the step size may be reduced to
  hit the endpoint since the code will not take a step larger than the `h`
  input.  Changing the direction of integration, i.e., the sign of h ,
  requires the user set `start = true` before calling `step` again.  This is
  the only situation in which `start` should be altered.

  # Output from `step`

  ## Successful step

  The function returns zero after each successful step with `start` set to
  `false`.  `x` represents the independent variable advanced one step of
  length hold from its value on input and `y` the solution vector at the new
  value of `x`.  All other parameters represent information corresponding to
  the new `x` needed to continue the integration.

  ## Unsuccessful step

  When the error tolerance is too small for the machine precision, the
  function returns nonzero without taking a step.  An appropriate step size
  and error tolerance for continuing are estimated and all other information
  is restored as upon input before returning.  To continue with the larger
  tolerance, the user just calls the code again.  A restart is neither
  required nor desirable.
*/
int step(double *const restrict y,
         const fn_type f,
         void *const restrict f_ctx,
         double *const restrict eps,
         double *const restrict wt,
         double *const restrict *const phi,
         double *const restrict p,
         double *const restrict yp,
         struct Ode *const self)
{
    struct SgVectorDriver drv = self->drv;
    static const double gstr[13] = {
        0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114,
        0.00936, 0.00789, 0.00679, 0.00592, 0.00524, 0.00468
    };

    const double p5eps = *eps * 0.5;
    bool *const start = &self->start;
    bool *const phase1 = &self->phase1;
    bool *const nornd = &self->nornd;
    unsigned *const k = &self->k;
    unsigned *const kold = &self->kold;
    unsigned *const ns = &self->ns;
    double *const x = &self->x;
    double *const h = &self->h;
    double *const hold = &self->hold;
    double *const psi = self->psi;
    double *const alpha = self->alpha;
    double *const beta = self->beta;
    double *const sig = self->sig;
    double *const v = self->v;
    double *const w = self->w;
    double *const g = self->g;

    unsigned knew;
    unsigned i, iq, j;
    int ifail;
    double absh, erk, erkm1, erkm2, erkp1, err, hnew, round, temp1;

    /**** begin block 0 ****/
    /* check if step size or error tolerance is too small for machine
       precision.  if first step, initialize phi array and estimate a starting
       step size. */

    /* if step size is too small, determine an acceptable one */
    if (fabs(*h) < 4.0 * DBL_EPSILON * fabs(*x)) {
        *h = copysign(4.0 * DBL_EPSILON * fabs(*x), *h);
        return 1;
    }

    /* if error tolerance is too small, increase it to an acceptable value */
    /* round = 2 ε ‖𝐲 / 𝛚‖ */
    round = 2.0 * DBL_EPSILON * sqrt(vector_div_normsq(drv, y, wt));
    if (p5eps < round) {
        *eps = round * 2.0 * (4.0 * DBL_EPSILON + 1.0);
        return 1;
    }
    g[0] = 1.0;
    g[1] = 0.5;
    sig[0] = 1.0;
    if (*start) {
        /* initialize.  compute appropriate step size for first step */
        (*f)(f_ctx, *x, y, yp);
        {
            double sum = 0.0;
            /* 𝛗[1] ← 𝟎 */
            sg_vector_fill(drv, 0.0, phi[1]);
            /* 𝛗[0] ← 𝐲′ */
            sg_vector_copy(drv, yp, phi[0]);
            /* sum = ‖𝐲′ / 𝛚‖ */
            sum = sqrt(vector_div_normsq(drv, yp, wt));
            absh = fabs(*h);
            if (*eps < sum * 16.0 * pow(*h, 2.0)) {
                absh = sqrt(*eps / sum) * 0.25;
            }
        }
        *h = copysign(max(absh, 4.0 * DBL_EPSILON * fabs(*x)), *h);
        *hold = 0.0;
        *k = 1;
        *kold = 0;
        *start = false;
        *phase1 = true;
        *nornd = true;
        if (p5eps <= round * 100.0) {
            *nornd = false;
            /* 𝛗[14] ← 𝟎 */
            sg_vector_fill(drv, 0.0, phi[14]);
        }
    }
    ifail = 0;
    /**** end block 0 ****/

    /**** begin block 1 ****/
    /* compute coefficients of formulas for this step.  avoid computing */
    /* those quantities not changed when step size is not changed. */
    while (1) {
        /* ns is the number of steps taken with size h, including the current
           one.  when k < ns, no coefficients change */
        if (*h != *hold) {
            *ns = 0;
        }
        if (*ns <= *kold) {
            ++*ns;
        }
        if (*k >= *ns) {
            /* PRE: ns ≥ 1 */

            temp1 = *h * *ns;
            /* compute those components of alpha, beta, psi, sig which are
               changed */
            beta[*ns - 1] = 1.0;
            alpha[*ns - 1] = 1.0 / *ns;
            sig[*ns] = 1.0;
            for (i = *ns; i < *k; ++i) {
                const double psiim1 = psi[i - 1];
                psi[i - 1] = temp1;
                beta[i] = beta[i - 1] * psi[i - 1] / psiim1;
                temp1 = psiim1 + *h;
                alpha[i] = *h / temp1;
                sig[i + 1] = (double)(i + 1) * alpha[i] * sig[i];
            }
            psi[*k - 1] = temp1;

            /* compute coefficients g */

            /* initialize v and set w.  g[1] is set in data statement */
            if (*ns <= 1) {
                for (iq = 0; iq < *k; ++iq) {
                    v[iq] = 1.0 / ((iq + 1) * (iq + 2));
                    w[iq] = v[iq];
                }
            } else {
                /* if order was raised, update diagonal part of v */
                if (*k > *kold) {
                    v[*k - 1] = 1.0 / (*k * (*k + 1));
                    for (j = 1; j < *ns - 1; ++j) {
                        const unsigned i = *k - j;
                        v[i - 1] -= alpha[j] * v[i];
                    }
                }
                /* update v and set w */
                const double temp5 = alpha[*ns - 1];
                for (iq = 0; iq <= *k - *ns; ++iq) {
                    v[iq] -= temp5 * v[iq + 1];
                    w[iq] = v[iq];
                }
                g[*ns] = w[0];
            }

            /* compute the g in the work vector w */
            for (i = *ns; i < *k; ++i) {
                for (iq = 0; iq < *k - i; ++iq) {
                    w[iq] -= alpha[i] * w[iq + 1];
                }
                g[i + 1] = w[0];
            }
        }
        /**** end block 1 ****/

        /**** begin block 2 ****/
        /* predict a solution p, evaluate derivatives using predicted
           solution, estimate local error at order k and errors at orders k,
           k-1, k-2 as if constant step size were used. */

        /* change phi to phi star */
        for (i = *ns; i < *k; ++i) {
            /* 𝛗[i] ← β[i] 𝛗[i] */
            sg_vector_scale_assign(drv, beta[i], phi[i]);
        }
        /* predict solution and differences */
        /* 𝛗[k + 1] ← 𝛗[k] */
        sg_vector_copy(drv, phi[*k], phi[*k + 1]);
        /* 𝛗[k] ← 𝟎 */
        sg_vector_fill(drv, 0.0, phi[*k]);
        /* 𝐩 ← 𝛗 𝐠 (matrix-vector) */
        sg_vector_fill(drv, 0.0, p);
        for (i = *k; i-- > 0;) {
            /* 𝐩 ← 𝐩 + g[i] 𝛗[i] */
            sg_vector_linear_assign(drv, 1.0, g[i], phi[i], p);
        }
        for (i = *k; i-- > 0;) {
            /* 𝛗[i] ← 𝛗[i] + 𝛗[i + 1] */
            sg_vector_linear_assign(drv, 1.0, 1.0, phi[i + 1], phi[i]);
        }
        if (!*nornd) {
            /* 𝐩 ← 𝛕′(h, 𝐲, 𝛗[14], 𝛗[15], 𝐩) */
            vector_taup(drv, *h, y, phi[14], phi[15], p);
        } else {
            /* 𝐩 ← 𝐲 + h 𝐩 */
            sg_vector_linear_assign(drv, *h, 1.0, y, p);
        }
        {
            const double xold = *x;
            *x += *h;
            absh = fabs(*h);
            (*f)(f_ctx, *x, p, yp);

            /* estimate errors at orders k, k-1, k-2 */
            /* erkm2 = ‖(𝛗[k - 2] + 𝐲′ - 𝛗[0]) / 𝛚‖²
               erkm1 = ‖(𝛗[k - 1] + 𝐲′ - 𝛗[0]) / 𝛚‖²
               erk = ‖(𝐲′ - 𝛗[0]) / 𝛚‖² */
            // TODO: remove cast
            vector_erk(drv, *k, wt, yp, (SgVector *const *)phi, &erk, &erkm1, &erkm2);
            if (*k > 2) {
                erkm2 = absh * sig[*k - 2] * gstr[*k - 3] * sqrt(erkm2);
            }
            if (*k >= 2) {
                erkm1 = absh * sig[*k - 1] * gstr[*k - 2] * sqrt(erkm1);
            }
            const double temp5 = absh * sqrt(erk);
            err = temp5 * (g[*k - 1] - g[*k]);
            erk = temp5 * sig[*k] * gstr[*k - 1];
            knew = *k;

            /* test if order should be lowered */
            if (*k == 2) {
                if (erkm1 <= erk * 0.5) {
                    knew = *k - 1;
                }
            } else if (*k > 2) {
                if (max(erkm1, erkm2) <= erk) {
                    knew = *k - 1;
                }
            }

            /* test if step successful */
            if (err <= *eps) {
                break;
            }
            /**** end block 2 ****/

            /**** begin block 3 ****/
            /* the step is unsuccessful.  restore x, phi, psi.  if third
               consecutive failure, set order to one.  if step fails more than
               three times, consider an optimal step size.  double error
               tolerance and return if estimated step size is too small for
               machine precision. */

            /* restore x, phi and psi */
            *phase1 = false;
            *x = xold;
        }
        for (i = 0; i < *k; ++i) {
            vector_ibeta_diff(drv, 1.0 / beta[i], phi[i + 1], phi[i]);
        }
        for (i = 1; i < *k; ++i) {
            psi[i - 1] = psi[i] - *h;
        }
        /* on third failure, set order to one.  thereafter, use optimal step
           size */
        ++ifail;
        {
            double temp2 = 0.5;
            if (ifail >= 3) {
                if (ifail != 3 && p5eps < erk * 0.25) {
                    temp2 = sqrt(p5eps / erk);
                }
                knew = 1;
            }
            *h = temp2 * *h;
        }
        *k = knew;
        if (fabs(*h) < 4.0 * DBL_EPSILON * fabs(*x)) {
            *h = copysign(4.0 * DBL_EPSILON * fabs(*x), *h);
            *eps += *eps;
            return 1;
        }
    }
    /**** end block 3 ****/

    /**** begin block 4 ****/
    /* the step is successful.  correct the predicted solution, evaluate the
       derivatives using the corrected solution and update the differences.
       determine best order and step size for next step. */
    *kold = *k;
    *hold = *h;

    /* correct and evaluate */
    {
        const double hgk = *h * g[*k];
        if (!*nornd) {
            vector_rhop(drv, hgk, yp, phi[0], phi[15], p, y, phi[14]);
        } else {
            vector_eval_y(drv, p, hgk, yp, phi[0], y);
        }
    }
    (*f)(f_ctx, *x, y, yp);

    /* update differences for next step */
    /* 𝛗[k] ← 𝐲′ - 𝛗[0]
       𝛗[k + 1] ← 𝛗[k] - 𝛗[k + 1] */
    vector_upd_diffs(drv, yp, phi[0], phi[*k], phi[*k + 1]);
    /* ∀ i ∈ [0, k).  𝛗[i] ← 𝛗[i] + 𝛗[k] */
    for (i = 0; i < *k; ++i) {
        sg_vector_linear_assign(drv, 1.0, 1.0, phi[*k], phi[i]);
    }

    /* estimate error at order k+1 unless: in first phase when always raise
       order, already decided to lower order, step size not constant so
       estimate unreliable */

    erkp1 = 0.0;
    if (knew == *k - 1 || *k == 12) {
        *phase1 = false;
    }

    if (*phase1) {
        ++*k;
        erk = erkp1;
    } else if (knew == *k - 1) {
        --*k;
        erk = erkm1;
    } else if (*k < *ns) {
        /* erkp1 = ‖𝛗[k + 1] / 𝛚‖² */
        erkp1 = absh * gstr[*k] * sqrt(vector_div_normsq(drv, phi[*k + 1], wt));

        /* using estimated error at order k+1, determine appropriate order
           for next step */

        if (*k > 1) {
            if (erkm1 <= min(erk, erkp1)) {
                --*k;
                erk = erkm1;
            } else if (erkp1 < erk && *k != 12) {
                ++*k;
                erk = erkp1;
            }
        } else if (erkp1 < erk * 0.5) {
            ++*k;
            erk = erkp1;
        }
    }

    /* with new order determine appropriate step size for next step */
    hnew = 2.0 * *h;
    if (!*phase1 && p5eps < erk * pow(2.0, *k + 1)) {
        hnew = *h;
        if (p5eps >= erk) {
            return 0;
        }
        hnew = absh * max(0.5, min(0.9, pow(p5eps / erk, 1.0 / (*k + 1))));
        hnew = copysign(max(hnew, 4.0 * DBL_EPSILON * fabs(*x)), *h);
    }
    *h = hnew;
    /**** end block 4 ****/
    return 0;
}

/*
  The methods in function `step` approximate the solution near `x` by a
  polynomial.  Function `intrp` approximates the solution at `xout` by
  evaluating the polynomial there.  Information defining this polynomial is
  passed from `step` so `intrp` cannot be used alone.

  # Input to `intrp`

  The user provides storage in the calling program for the arrays in the call
  list and defines `xout`, point at which solution is desired.

  The remaining parameters are defined in `step` and passed to `intrp` from
  that function.

  # Output from `intrp`

    - `yout`: Solution at `xout`
    - `ypout`: Derivative of solution at `xout`

  The remaining parameters are returned unaltered from their input values.
  Integration with `step` may be continued.
*/
void intrp(const struct SgVectorDriver drv,
           const double x,
           const double *const restrict y,
           const double xout,
           double *const restrict yout,
           double *const restrict ypout,
           const unsigned kold,
           double *const restrict *const phi,
           const double *const restrict psi)
{
    const double hi = xout - x;
    const unsigned ki = (unsigned)(kold + 1);

    unsigned i, j;
    double term = 0.0;
    double g[13] = {1.0};
    double rho[13] = {1.0};
    double w[13] = {0.0};

    if (kold >= 13) {
        fprintf(stderr, "invalid kold\n");
        fflush(stderr);
        abort();
    }

    /* initialize w for computing g */
    for (i = 0; i < ki; ++i) {
        w[i] = 1.0 / (i + 1);
    }

    /* compute g */
    for (j = 1; j < ki; ++j) {
        const double psijm1 = psi[j - 1];
        const double gamma = (hi + term) / psijm1;
        const double eta = hi / psijm1;
        for (i = 0; i < ki - j; ++i) {
            w[i] = gamma * w[i] - eta * w[i + 1];
        }
        g[j] = w[0];
        rho[j] = gamma * rho[j - 1];
        term = psijm1;
    }

    /* interpolate */
    /* 𝐲° ← 𝛗 𝐠 (matrix-vector)
       𝐲°′ ← 𝛗 𝛒 (matrix-vector) */
    sg_vector_fill(drv, 0.0, ypout);
    sg_vector_fill(drv, 0.0, yout);
    for (i = ki; i-- > 0;) {
        vector_intrp_yout(drv, g[i], rho[i], phi[i], yout, ypout);
    }
    /* 𝐲° ← 𝐲 + hi 𝐲° */
    sg_vector_linear_assign(drv, hi, 1.0, y, yout);
}

/*
  `ode` merely allocates storage for `de` to relieve the user of the
  inconvenience of a long call list.  Consequently `de` is used as described
  in the comments for `ode`.

  The constant `maxnum` is the maximum number of steps allowed in one call to
  `de`.
*/
void de(struct Ode *const self,
        const fn_type f,
        void *const restrict f_ctx,
        double *const restrict y,
        double *const restrict t,
        const double tout,
        double *const restrict relerr,
        double *const restrict abserr,
        const unsigned maxnum,
        int *const restrict iflag)
{
    struct SgVectorDriver drv = self->drv;
    const bool isn = *iflag >= 0;
    const double del = tout - *t;
    const double absdel = fabs(del);
    const double tend = isn ? *t + del * 10.0 : tout;

    bool stiff;
    double abseps, eps, releps;
    int kle4;
    unsigned nostep;

    /* trivial case */
    if (sg_vector_len(self->drv) == 0) {
        *iflag = 2;
        *t = tout;
        return;
    }

    /* test for improper parameters */
    eps = max(*relerr, *abserr);
    *iflag = abs(*iflag);
    if (*t == tout ||
        *relerr < 0.0 || *abserr < 0.0 ||
        eps <= 0.0 || *iflag == 0 ||
        (*iflag != 1 && (*t != self->told || *iflag < 2 || *iflag > 5))) {
        *iflag = 6;
        return;
    }

    /* on each call set interval of integration and counter for number of
       steps.  adjust input error tolerances to define weight vector for
       `step` */
    kle4 = 0;
    stiff = false;
    releps = *relerr / eps;
    abseps = *abserr / eps;
    if (*iflag == 1 || !self->isnold || self->delsgn * del <= 0.0) {
        /* on start and restart also set work variables x and yy, store the
           direction of integration and initialize the step size */
        self->start = true;
        self->x = *t;
        /* 𝐘 ← 𝐲 */
        sg_vector_copy(drv, y, self->yy);
        self->delsgn = copysign(1.0, del);
        self->h = copysign(max(fabs(tout - self->x),
                               4.0 * DBL_EPSILON * fabs(self->x)),
                           tout - self->x);
    }

    /* if already past output point, interpolate and return */
    for (nostep = 0;; ++nostep) {

        if (fabs(self->x - *t) >= absdel) {
            /* (𝐲, 𝐲°′) ← intrp(x, 𝐘, tout, kold, 𝛗, 𝛙) */
            intrp(self->drv, self->x, self->yy, tout, y, self->ypout,
                  self->kold, self->phi, self->psi);
            *iflag = 2;
            *t = tout;
            self->told = *t;
            self->isnold = isn;
            return;
        }

        /* if cannot go past output point and sufficiently close, extrapolate and
           return */
        if (!isn || fabs(tout - self->x) < 4.0 * DBL_EPSILON * fabs(self->x)) {
            self->h = tout - self->x;
            (*f)(f_ctx, self->x, self->yy, self->yp);
            /* 𝐲 ← 𝐘 + h 𝐲′ */
            sg_vector_linear(drv, 1.0, self->yy, self->h, self->yp, y);
            *iflag = 2;
            *t = tout;
            self->told = *t;
            self->isnold = isn;
            return;
        }

        /* test for too many steps */
        if (nostep >= maxnum) {
            *iflag = isn ? 4 : -4;
            if (stiff) {
                *iflag = isn ? 5 : -5;
            }
            /* 𝐲 ← 𝐘 */
            sg_vector_copy(drv, self->yy, y);
            *t = self->x;
            self->told = *t;
            self->isnold = true;
            return;
        }

        /* limit step size, set weight vector and take a step */
        self->h = copysign(min(fabs(self->h), fabs(tend - self->x)), self->h);
        /* 𝛚 ← releps |𝐘| + abseps */
        vector_update_wt(drv, releps, self->yy, abseps, self->wt);

        /* test for tolerances too small */
        if (step(self->yy, f, f_ctx, &eps, self->wt,
                 self->phi, self->p, self->yp, self)) {
            *iflag = isn ? 3 : -3;
            *relerr = eps * releps;
            *abserr = eps * abseps;
            /* 𝐲 ← 𝐘 */
            sg_vector_copy(drv, self->yy, y);
            *t = self->x;
            self->told = *t;
            self->isnold = true;
            return;
        }

        /* augment counter on number of steps and test for stiffness */
        ++kle4;
        if (self->kold > 4) {
            kle4 = 0;
        }
        if (kle4 >= 50) {
            stiff = true;
        }
    }
}

/*
  Integrates a system of `neqn` first order ordinary differential equations of
  the form:

      dy[i]/dt = f(t, y[0], y[1], ..., y[neqn - 1])
      y[i] given at `t`

  The function integrates from `t` to `tout`.  On return the parameters in the
  call list are set for continuing the integration.  The user has only to
  define a new value `tout` and call `ode` again.

  The differential equations are actually solved by a suite of codes `de`,
  `step`, and `intrp`.  `ode` allocates virtual storage in the `work` array
  and calls `de`.  `de` is a supervisor which directs the solution.  It calls
  on the routines `step` and `intrp` to advance the integration and to
  interpolate at output points.  `step` uses a modified divided difference
  form of the Adams PECE formulas and local extrapolation.  It adjusts the
  order and step size to control the local error per unit step in a
  generalized sense.  Normally each call to `step` advances the solution one
  step in the direction of `tout`.  For reasons of efficiency `de` integrates
  beyond `tout` internally, though never beyond `t + 10 * (tout - t)`, and
  calls `intrp` to interpolate the solution at `tout`.  An option is provided
  to stop the integration at `tout` but it should be used only if it is
  impossible to continue the integration beyond `tout`.

  This code is completely explained and documented in the text, Computer
  Solution of Ordinary Differential Equations: The Initial Value Problem
  L. F. Shampine and M. K. Gordon.

  @param f
  Function `f(f_ctx, t, y, yp)` to evaluate derivatives `yp[i] = dy[i]/dt`

  @param y
  Solution vector at `t` (length: `neqn`)

  @param t
  Independent variable (length: `1`)

  @param tout
  Point at which solution is desired

  @param relerr
  Relative error tolerance for local error test (length: `1`).
  At each step the code requires
  `fabs(local_error) <= fabs(y) * relerr + abserr`
  for each component of the local error and solution vectors

  @param abserr
  Absolute error tolerance for local error test (length: `1`).
  See `relerr`.

  @param iflag
  Indicates status of integration

  @param self
  Arrays to hold information internal to `de` which is necessary for
  subsequent calls (length: `5`)

  # First call to `ode`

  The user must supply a function `f(f_ctx, t, y, yp)` to evaluate

      dy[i]/dt = yp[i] = f(t, y[0], y[1], ..., y[neqn - 1])

  and initialize the parameters:

    - `y`: Vector of initial conditions
    - `t`: Starting point of integration
    - `tout`: Point at which solution is desired
    - `relerr, abserr`: Relative and absolute local error tolerances
    - `iflag`: `+1` or `-1`. Indicator to initialize the code.  Normal input
      is `+1`.  The user should set `iflag = -1` only if it is impossible to
      continue the integration beyond `tout`.

  All parameters except `f` and `tout` may be altered by the code on output so
  must be variables in the calling program.

  # Output from `ode`

    - `y`: Solution at `t`
    - `t`: Last point reached in integration.  normal return has `t == tout`.
    - `tout`: Unchanged
    - `relerr`, `abserr`: normal return has tolerances unchanged.  `iflag = 3`
      signals tolerances increased
    - `iflag`: (Note: The value of `iflag` is returned negative when the input
      value is negative and the integration does not reach `tout`, i.e., `-3`,
      `-4`, `-5`.)
        - `2`: Normal return.  Integration reached `tout`.
        - `3`: Integration did not reach `tout` because error tolerances too
          small.  `relerr`, `abserr` increased appropriately for continuing.
        - `4`: Integration did not reach `tout` because more than `maxnum`
           steps needed.
        - `5`: Integration did not reach `tout` because equations appear to be
          stiff.
        - `6`: Invalid input parameters (fatal error).
    - `work`, `self`: Information generally of no interest to the user but
      necessary for subsequent calls.

  # Subsequent calls to `ode`

  Function `ode` returns with all information needed to continue the
  integration.  If the integration reached `tout`, the user need only define a
  new `tout` and call again.  If the integration did not reach `tout` and the
  user wants to continue, simply call again.  The output value of `iflag` is
  the appropriate input value for subsequent calls.  The only situation in
  which it should be altered is to stop the integration internally at the new
  `tout`, i.e., change output `iflag = 2` to input `iflag = -2`.  Error
  tolerances may be changed by the user before continuing.  All other
  parameters must remain unchanged.
*/
void ode(struct Ode *const self,
         const fn_type f,
         void *const restrict f_ctx,
         SgVector *const restrict y,
         double *const restrict t,
         const double tout,
         double *const restrict relerr,
         double *const restrict abserr,
         const unsigned maxnum,
         int *const restrict iflag)
{
    de(self, f, f_ctx, y, t, tout, relerr, abserr, maxnum, iflag);
}

void ode_init(struct Ode *self, struct SgVectorDriver drv)
{
    size_t i;
    self->drv = drv;
    self->yy = sg_vector_new(self->drv);
    self->wt = sg_vector_new(self->drv);
    self->p = sg_vector_new(self->drv);
    self->yp = sg_vector_new(self->drv);
    self->ypout = sg_vector_new(self->drv);
    for (i = 0; i < sizeof(self->phi) / sizeof(*self->phi); ++i) {
        self->phi[i] = sg_vector_new(self->drv);
    }
}

void ode_del(struct Ode *self)
{
    size_t i;
    for (i = 0; i < sizeof(self->phi) / sizeof(*self->phi); ++i) {
        sg_vector_del(self->drv, self->phi[i]);
    }
    sg_vector_del(self->drv, self->ypout);
    sg_vector_del(self->drv, self->yp);
    sg_vector_del(self->drv, self->p);
    sg_vector_del(self->drv, self->wt);
    sg_vector_del(self->drv, self->yy);
}
