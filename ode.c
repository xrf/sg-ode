/* Original: http://www.netlib.org/ode/ode.f */
/* Converted through f2c and then manually refactored. */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ode.h"

typedef void (*fn_type)(void *restrict f_ctx,
                        double t,
                        const double *restrict y,
                        double *restrict yp);

static double min(double x, double y)
{
    return x < y ? x : y;
}

static double max(double x, double y)
{
    return x >= y ? x : y;
}

static void clear_double_array(double *a, size_t n) {
    // compilers should be able to optimize this to a simple memset (can't use
    // a real memset because C standard doesn't guarantee that would work)
    size_t i;
    for (i = 0; i < n; ++i) {
        a[i] = 0.0;
    }
}

static void copy_double_array(double *out, const double *in, size_t n) {
    memcpy(out, in, n * sizeof(*in));
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

  @param neqn
  Number of equations to be integrated

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

  The arrays `phi`, `psi` are required for the interpolation subroutine
  `intrp`.  The array `p` is internal to the code.

  # Input to `step`

  ## First call

  The user must provide storage in their driver program for all arrays in the
  call list, namely

      y[neqn], wt[neqn], phi[neqn * 16], p[neqn], yp[neqn], psi[12]

  The user must also declare the `start` Boolean variable and `f` an external
  subroutine, supply the subroutine `f(f_ctx, x, y, yp)` to evaluate

      dy[i]/dx = yp[i] = f(x, y[0], y[1], ..., y[neqn - 1])

  and initialize only the following parameters:

    - `x`: Initial value of the independent variable
    - `y`: Vector of initial values of dependent variables
    - `neqn`: Number of equations to be integrated
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

  Subroutine `step` is designed so that all information needed to continue the
  integration, including the step size `h` and the order `k`, is returned with
  each step.  With the exception of the step size, the error tolerance, and
  the weights, none of the parameters should be altered.  The array `wt` must
  be updated after each step to maintain relative error tests like those
  above.  Normally the integration is continued just beyond the desired
  endpoint and the solution interpolated there with subroutine `intrp`.  If it
  is impossible to integrate beyond the endpoint, the step size may be reduced
  to hit the endpoint since the code will not take a step larger than the `h`
  input.  Changing the direction of integration, i.e., the sign of h ,
  requires the user set `start = true` before calling `step` again.  This is
  the only situation in which `start` should be altered.

  # Output from `step`

  ## Successful step

  The subroutine returns zero after each successful step with `start`
  set to `false`.  `x` represents the independent variable advanced one step
  of length hold from its value on input and `y` the solution vector at the
  new value of `x`.  All other parameters represent information corresponding
  to the new `x` needed to continue the integration.

  ## Unsuccessful step

  When the error tolerance is too small for the machine precision, the
  subroutine returns nonzero without taking a step.  An appropriate step size
  and error tolerance for continuing are estimated and all other information
  is restored as upon input before returning.  To continue with the larger
  tolerance, the user just calls the code again.  A restart is neither
  required nor desirable.
*/
int step(double *const restrict y,
         const fn_type f,
         void *const restrict f_ctx,
         const size_t neqn,
         double *const restrict eps,
         double *const restrict wt,
         double *const restrict phi,
         double *const restrict p,
         double *const restrict yp,
         struct Ode *const self)
{
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

    size_t l;
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
    /* round = ‖𝐲 / 𝛚‖² */
    round = 0.0;
    for (l = 0; l < neqn; ++l) {
        round += pow(y[l] / wt[l], 2.0);
    }
    round = 2.0 * DBL_EPSILON * sqrt(round);
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
            clear_double_array(phi + neqn, neqn);
            /* 𝛗[0] ← 𝐲′ */
            copy_double_array(phi, yp, neqn);
            /* sum = ‖𝐲′ / 𝛚‖² */
            for (l = 0; l < neqn; ++l) {
                sum += pow(yp[l] / wt[l], 2.0);
            }
            sum = sqrt(sum);
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
            clear_double_array(&phi[neqn * 14], neqn);
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
            /* PRE: ns >= 1 */

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
            for (l = 0; l < neqn; ++l) {
                phi[l + i * neqn] *= beta[i];
            }
        }
        /* predict solution and differences */
        /* 𝛗[k + 1] ← 𝛗[k] */
        copy_double_array(phi + (*k + 1) * neqn, phi + *k * neqn, neqn);
        /* 𝛗[k] ← 𝟎 */
        clear_double_array(phi + *k * neqn, neqn);
        /* 𝐩 ← 𝛗 𝐠 (matrix-vector) */
        clear_double_array(p, neqn);
        for (i = *k; i-- > 0;) {
            for (l = 0; l < neqn; ++l) {
                p[l] += g[i] * phi[l + i * neqn];
            }
        }
        for (i = *k; i-- > 0;) {
            /* 𝛗[i] ← 𝛗[i] + 𝛗[i + 1] */
            for (l = 0; l < neqn; ++l) {
                phi[l + i * neqn] += phi[l + (i + 1) * neqn];
            }
        }
        if (!*nornd) {
            /* 𝛕 = h 𝐩 - 𝛗[14]
               𝐩 ← 𝐲 + 𝛕
               𝛗[15] ← (𝐩 - 𝐲) - 𝛕 */
            for (l = 0; l < neqn; ++l) {
                const double tau = *h * p[l] - phi[l + neqn * 14];
                p[l] = y[l] + tau;
                phi[l + neqn * 15] = (p[l] - y[l]) - tau;
            }
        } else {
            /* 𝐩 ← 𝐲 + h 𝐩 */
            for (l = 0; l < neqn; ++l) {
                p[l] = y[l] + *h * p[l];
            }
        }
        {
            const double xold = *x;
            *x += *h;
            absh = fabs(*h);
            (*f)(f_ctx, *x, p, yp);

            /* estimate errors at orders k, k-1, k-2 */
            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;
            /* erkm2 = ‖(𝛗[k - 2] + 𝐲′ - 𝛗[0]) / 𝛚‖²
               erkm1 = ‖(𝛗[k - 1] + 𝐲′ - 𝛗[0]) / 𝛚‖²
               erk = ‖(𝐲′ - 𝛗[0]) / 𝛚‖² */
            for (l = 0; l < neqn; ++l) {
                const double iwt = 1.0 / wt[l];
                const double ypmphi = yp[l] - phi[l];
                if (*k > 2) {
                    erkm2 += pow((phi[l + (*k - 2) * neqn] + ypmphi) * iwt, 2.0);
                }
                if (*k >= 2) {
                    erkm1 += pow((phi[l + (*k - 1) * neqn] + ypmphi) * iwt, 2.0);
                }
                erk += pow(ypmphi * iwt, 2.0);
            }
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
            /* 𝛗[i] ← β[i]⁻¹ (𝛗[i] - 𝛗[i + 1]) */
            for (l = 0; l < neqn; ++l) {
                phi[l + i * neqn] =
                    (1.0 / beta[i]) *
                    (phi[l + i * neqn] - phi[l + (i + 1) * neqn]);
            }
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
            /* 𝛒 = h g[k] (𝐲′ - 𝛗[0]) - 𝛗[15]
               𝐲 ← 𝐩 + 𝛒
               𝛗[14] ← (𝐲 - 𝐩) - 𝛒 */
            for (l = 0; l < neqn; ++l) {
                const double rho = hgk * (yp[l] - phi[l]) - phi[l + neqn * 15];
                y[l] = p[l] + rho;
                phi[l + neqn * 14] = (y[l] - p[l]) - rho;
            }
        } else {
            /* 𝐲 ← 𝐩 + h g[k] (𝐲′ - 𝛗[0]) */
            for (l = 0; l < neqn; ++l) {
                y[l] = p[l] + hgk * (yp[l] - phi[l]);
            }
        }
    }
    (*f)(f_ctx, *x, y, yp);

    /* update differences for next step */
    /* 𝛗[k] ← 𝐲′ - 𝛗[0]
       𝛗[k + 1] ← 𝛗[k] - 𝛗[k + 1] */
    for (l = 0; l < neqn; ++l) {
        phi[l + *k * neqn] = yp[l] - phi[l];
        phi[l + (*k + 1) * neqn] =
            phi[l + *k * neqn] -
            phi[l + (*k + 1) * neqn];
    }
    /* ∀ i ∈ [0, k).  𝛗[i] ← 𝛗[i] + 𝛗[k] */
    for (i = 0; i < *k; ++i) {
        for (l = 0; l < neqn; ++l) {
            phi[l + i * neqn] += phi[l + *k * neqn];
        }
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
        for (l = 0; l < neqn; ++l) {
            erkp1 += pow(phi[l + (*k + 1) * neqn] / wt[l], 2.0);
        }
        erkp1 = absh * gstr[*k] * sqrt(erkp1);

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
  The methods in subroutine `step` approximate the solution near `x` by a
  polynomial.  Subroutine `intrp` approximates the solution at `xout` by
  evaluating the polynomial there.  Information defining this polynomial is
  passed from `step` so `intrp` cannot be used alone.

  # Input to `intrp`

  The user provides storage in the calling program for the arrays in the call
  list and defines `xout`, point at which solution is desired.

  The remaining parameters are defined in `step` and passed to `intrp` from
  that subroutine.

  # Output from `intrp`

    - `yout`: Solution at `xout`
    - `ypout`: Derivative of solution at `xout`

  The remaining parameters are returned unaltered from their input values.
  Integration with `step` may be continued.
*/
void intrp(const double x,
           const double *const restrict y,
           const double xout,
           double *const restrict yout,
           double *const restrict ypout,
           const size_t neqn,
           const unsigned kold,
           const double *const restrict phi,
           const double *const restrict psi)
{
    const double hi = xout - x;
    const unsigned ki = (unsigned)(kold + 1);

    size_t l;
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
    clear_double_array(ypout, neqn);
    clear_double_array(yout, neqn);
    for (i = ki; i-- > 0;) {
        const double gi = g[i];
        const double rhoi = rho[i];
        for (l = 0; l < neqn; ++l) {
            yout[l] += gi * phi[l + i * neqn];
            ypout[l] += rhoi * phi[l + i * neqn];
        }
    }
    /* 𝐲° ← 𝐲 + hi 𝐲° */
    for (l = 0; l < neqn; ++l) {
        yout[l] = y[l] + hi * yout[l];
    }
}

/*
  `ode` merely allocates storage for `de` to relieve the user of the
  inconvenience of a long call list.  Consequently `de` is used as described
  in the comments for `ode`.

  The constant `maxnum` is the maximum number of steps allowed in one call to
  `de`.
*/
void de(const fn_type f,
        void *const restrict f_ctx,
        const size_t neqn,
        double *const restrict y,
        double *const restrict t,
        const double tout,
        double *const restrict relerr,
        double *const restrict abserr,
        int *const restrict iflag,
        double *const restrict yy,
        double *const restrict wt,
        double *const restrict p,
        double *const restrict yp,
        double *const restrict ypout,
        double *const restrict phi,
        struct Ode *const self,
        const int maxnum)
{
    const bool isn = *iflag >= 0;
    const double del = tout - *t;
    const double absdel = fabs(del);
    const double tend = isn ? *t + del * 10.0 : tout;

    bool stiff;
    double abseps, eps, releps;
    int kle4, nostep;
    size_t l;

    /* test for improper parameters */
    eps = max(*relerr, *abserr);
    *iflag = abs(*iflag);
    if (neqn < 1 || *t == tout ||
        *relerr < 0.0 || *abserr < 0.0 ||
        eps <= 0.0 || *iflag == 0 ||
        (*iflag != 1 && (*t != self->told || *iflag < 2 || *iflag > 5))) {
        *iflag = 6;
        return;
    }

    /* on each call set interval of integration and counter for number of
       steps.  adjust input error tolerances to define weight vector for
       subroutine step */
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
        copy_double_array(yy, y, neqn);
        self->delsgn = copysign(1.0, del);
        self->h = copysign(max(fabs(tout - self->x),
                               4.0 * DBL_EPSILON * fabs(self->x)),
                           tout - self->x);
    }

    /* if already past output point, interpolate and return */
    for (nostep = 0;; ++nostep) {

        if (fabs(self->x - *t) >= absdel) {
            /* (𝐲, 𝐲°′) ← intrp(x, 𝐘, tout, kold, 𝛗, 𝛙) */
            intrp(self->x, yy, tout, y, ypout, neqn,
                  self->kold, phi, self->psi);
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
            (*f)(f_ctx, self->x, yy, yp);
            /* 𝐲 ← 𝐘 + h 𝐲′ */
            for (l = 0; l < neqn; ++l) {
                y[l] = yy[l] + self->h * yp[l];
            }
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
            copy_double_array(y, yy, neqn);
            *t = self->x;
            self->told = *t;
            self->isnold = true;
            return;
        }

        /* limit step size, set weight vector and take a step */
        self->h = copysign(min(fabs(self->h), fabs(tend - self->x)), self->h);
        /* 𝛚 ← releps |𝐘| + abseps */
        for (l = 0; l < neqn; ++l) {
            wt[l] = releps * fabs(yy[l]) + abseps;
        }

        /* test for tolerances too small */
        if (step(yy, f, f_ctx, neqn, &eps, wt, phi, p, yp, self)) {
            *iflag = isn ? 3 : -3;
            *relerr = eps * releps;
            *abserr = eps * abseps;
            /* 𝐲 ← 𝐘 */
            copy_double_array(y, yy, neqn);
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

  The subroutine integrates from `t` to `tout`.  On return the parameters in
  the call list are set for continuing the integration.  The user has only to
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
  Subroutine `f(f_ctx, t, y, yp)` to evaluate derivatives `yp[i] = dy[i]/dt`

  @param neqn
  Number of equations to be integrated

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

  @param work
  Arrays to hold information internal to `de` which is necessary for
  subsequent calls (length: `21 * neqn`)

  @param self
  Arrays to hold information internal to `de` which is necessary for
  subsequent calls (length: `5`)

  # First call to `ode`

  The user must provide storage in their calling program for the arrays
  in the call list,

      y[neqn], work[21 * neqn], self

  Supply supply the subroutine `f(f_ctx, t, y, yp)` to evaluate

      dy[i]/dt = yp[i] = f(t, y[0], y[1], ..., y[neqn - 1])

  and initialize the parameters:

    - `neqn`: Number of equations to be integrated
    - `y`: Vector of initial conditions
    - `t`: Starting point of integration
    - `tout`: Point at which solution is desired
    - `relerr, abserr`: Relative and absolute local error tolerances
    - `iflag`: `+1` or `-1`. Indicator to initialize the code.  Normal input
      is `+1`.  The user should set `iflag = -1` only if it is impossible to
      continue the integration beyond `tout`.

  All parameters except `f`, `neqn` and `tout` may be altered by the code on
  output so must be variables in the calling program.

  # Output from `ode`

    - `neqn`: Unchanged
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

  Subroutine `ode` returns with all information needed to continue the
  integration.  If the integration reached `tout`, the user need only define a
  new `tout` and call again.  If the integration did not reach `tout` and the
  user wants to continue, simply call again.  The output value of `iflag` is
  the appropriate input value for subsequent calls.  The only situation in
  which it should be altered is to stop the integration internally at the new
  `tout`, i.e., change output `iflag = 2` to input `iflag = -2`.  Error
  tolerances may be changed by the user before continuing.  All other
  parameters must remain unchanged.
*/
void ode(const fn_type f,
         void *const restrict f_ctx,
         const size_t neqn,
         double *const restrict y,
         double *const restrict t,
         const double tout,
         double *const restrict relerr,
         double *const restrict abserr,
         int *const restrict iflag,
         double *const restrict work,
         struct Ode *const self,
         const int maxnum)
{
    const size_t iwt = neqn;
    const size_t ip = iwt + neqn;
    const size_t iyp = ip + neqn;
    const size_t iypout = iyp + neqn;
    const size_t iphi = iypout + neqn;
    de(f, f_ctx, neqn, y, t, tout, relerr, abserr, iflag,
       &work[0], &work[iwt], &work[ip], &work[iyp], &work[iypout],
       &work[iphi], self, maxnum);
}
