// origin: http://www.netlib.org/ode/ode.f
// f2c'ed

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
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
  Appropriate order for next step (determined by code)

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
int step(double *const restrict x,
         double *restrict y,
         const fn_type f,
         void *const restrict f_ctx,
         const int neqn,
         double *const restrict h,
         double *const restrict eps,
         double *wt,
         bool *const restrict start,
         double *const restrict hold,
         int *const restrict k,
         int *const restrict kold,
         double *restrict phi,
         double *restrict p,
         double *restrict yp,
         double *const restrict psi,
         double *const restrict alpha,
         double *const restrict beta,
         double *restrict sig,
         double *restrict v,
         double *restrict w,
         double *restrict g,
         bool *const restrict phase1,
         int *const restrict ns,
         bool *const restrict nornd)
{
    static const double gstr[13] = {
        .5, .0833, .0417, .0264, .0188, .0143, .0114,
        .00936, .00789, .00679, .00592, .00524, .00468
    };

    const double p5eps = *eps * .5;

    int i, ifail, iq, j, km1, km2, knew, kp1, kp2, l;
    double absh, erk, erkm1, erkm2, erkp1, err, hnew, round, temp1;

    /* Parameter adjustments */
    --yp;
    --p;
    phi -= 1 + neqn;
    --wt;
    --y;
    --sig;
    --v;
    --w;
    --g;

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
    round = 0.;
    for (l = 1; l <= neqn; ++l) {
        round += pow(y[l] / wt[l], 2.0);
    }
    round = 2.0 * DBL_EPSILON * sqrt(round);
    if (p5eps < round) {
        *eps = round * 2.f * (4.0 * DBL_EPSILON + 1.);
        return 1;
    }
    g[1] = 1.;
    g[2] = .5;
    sig[1] = 1.;
    if (*start) {
        /* initialize.  compute appropriate step size for first step */
        (*f)(f_ctx, *x, &y[1], &yp[1]);
        {
            double sum = 0.;
            for (l = 1; l <= neqn; ++l) {
                phi[l + neqn] = yp[l];
                phi[l + neqn * 2] = 0.;
                sum += pow(yp[l] / wt[l], 2.0);
            }
            sum = sqrt(sum);
            absh = fabs(*h);
            if (*eps < sum * 16. * pow(*h, 2.0)) {
                absh = sqrt(*eps / sum) * .25;
            }
        }
        *h = copysign(max(absh, 4.0 * DBL_EPSILON * fabs(*x)), *h);
        *hold = 0.;
        *k = 1;
        *kold = 0;
        *start = false;
        *phase1 = true;
        *nornd = true;
        if (p5eps <= round * 100.) {
            *nornd = false;
            clear_double_array(&phi[1 + neqn * 15], (size_t)neqn);
        }
    }
    ifail = 0;
    /**** end block 0 ****/

    /**** begin block 1 ****/
    /* compute coefficients of formulas for this step.  avoid computing */
    /* those quantities not changed when step size is not changed. */
    while (1) {
        int nsp1; /* single assignment */

        kp1 = *k + 1;
        kp2 = *k + 2;
        km1 = *k - 1;
        km2 = *k - 2;

        /* ns is the number of steps taken with size h, including the current
           one.  when k < ns, no coefficients change */
        if (*h != *hold) {
            *ns = 0;
        }
        if (*ns <= *kold) {
            ++*ns;
        }
        nsp1 = *ns + 1;
        if (*k >= *ns) {
            temp1 = *h * *ns;
            /* compute those components of alpha, beta, psi, sig which are
               changed */
            beta[*ns - 1] = 1.0;
            alpha[*ns - 1] = 1.0 / *ns;
            sig[nsp1] = 1.0;
            if (*k >= nsp1) {
                for (i = nsp1; i <= *k; ++i) {
                    const int im1 = i - 1;
                    const double temp2 = psi[im1 - 1];
                    psi[im1 - 1] = temp1;
                    beta[i - 1] = beta[im1 - 1] * psi[im1 - 1] / temp2;
                    temp1 = temp2 + *h;
                    alpha[i - 1] = *h / temp1;
                    sig[i + 1] = (double)i * alpha[i - 1] * sig[i];
                }
            }
            psi[*k - 1] = temp1;

            /* compute coefficients g */

            /* initialize v and set w.  g[1] is set in data statement */
            if (*ns <= 1) {
                for (iq = 1; iq <= *k; ++iq) {
                    v[iq] = 1.0 / (iq * (iq + 1));
                    w[iq] = v[iq];
                }
            } else {
                /* if order was raised, update diagonal part of v */
                if (*k > *kold) {
                    const int nsm2 = *ns - 2;
                    v[*k] = 1.0 / (*k * kp1);
                    if (nsm2 >= 1) {
                        for (j = 1; j <= nsm2; ++j) {
                            int i = *k - j;
                            v[i] -= alpha[j] * v[i + 1];
                        }
                    }
                }
                /* update v and set w */
                const double temp5 = alpha[*ns - 1];
                for (iq = 1; iq <= kp1 - *ns; ++iq) {
                    v[iq] -= temp5 * v[iq + 1];
                    w[iq] = v[iq];
                }
                g[nsp1] = w[1];
            }

            /* compute the g[] in the work vector w[] */
            if (kp1 >= *ns + 2) {
                for (i = *ns + 2; i <= kp1; ++i) {
                    const int limit2 = kp2 - i;
                    const double temp6 = alpha[i - 2];
                    for (iq = 1; iq <= limit2; ++iq) {
                        w[iq] -= temp6 * w[iq + 1];
                    }
                    g[i] = w[1];
                }
            }
        }
        /**** end block 1 ****/

        /**** begin block 2 ****/
        /* predict a solution p[], evaluate derivatives using predicted
           solution, estimate local error at order k and errors at orders k,
           k-1, k-2 as if constant step size were used. */

        /* change phi to phi star */
        if (*k >= nsp1) {
            for (i = nsp1; i <= *k; ++i) {
                for (l = 1; l <= neqn; ++l) {
                    phi[l + i * neqn] = beta[i - 1] * phi[l + i * neqn];
                }
            }
        }
        /* predict solution and differences */
        for (l = 1; l <= neqn; ++l) {
            phi[l + kp2 * neqn] = phi[l + kp1 * neqn];
            phi[l + kp1 * neqn] = 0.;
            p[l] = 0.;
        }
        for (j = 1; j <= *k; ++j) {
            const int i = kp1 - j;
            const int ip1 = i + 1;
            const double temp2 = g[i];
            for (l = 1; l <= neqn; ++l) {
                p[l] += temp2 * phi[l + i * neqn];
                phi[l + i * neqn] += phi[l + ip1 * neqn];
            }
        }
        if (!*nornd) {
            for (l = 1; l <= neqn; ++l) {
                const double tau = *h * p[l] - phi[l + neqn * 15];
                p[l] = y[l] + tau;
                phi[l + (neqn * 16)] = p[l] - y[l] - tau;
            }
        } else {
            for (l = 1; l <= neqn; ++l) {
                p[l] = y[l] + *h * p[l];
            }
        }
        {
            const double xold = *x;
            *x += *h;
            absh = fabs(*h);
            (*f)(f_ctx, *x, &p[1], &yp[1]);

            /* estimate errors at orders k, k-1, k-2 */
            erkm2 = 0.;
            erkm1 = 0.;
            erk = 0.;
            for (l = 1; l <= neqn; ++l) {
                const double temp3 = 1. / wt[l];
                const double temp4 = yp[l] - phi[l + neqn];
                if (km2 > 0) {
                    erkm2 += pow((phi[l + km1 * neqn] + temp4) * temp3, 2.0);
                }
                if (km2 >= 0) {
                    erkm1 += pow((phi[l + *k * neqn] + temp4) * temp3, 2.0);
                }
                erk += pow(temp4 * temp3, 2.0);
            }
            if (km2 > 0) {
                erkm2 = absh * sig[km1] * gstr[km2 - 1] * sqrt(erkm2);
            }
            if (km2 >= 0) {
                erkm1 = absh * sig[*k] * gstr[km1 - 1] * sqrt(erkm1);
            }
            const double temp5 = absh * sqrt(erk);
            err = temp5 * (g[*k] - g[kp1]);
            erk = temp5 * sig[kp1] * gstr[*k - 1];
            knew = *k;

            /* test if order should be lowered */
            if (km2 == 0) {
                if (erkm1 <= erk * .5) {
                    knew = km1;
                }
            } else if (km2 > 0) {
                if (max(erkm1, erkm2) <= erk) {
                    knew = km1;
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
               three times, consider an optimal step size.  double error tolerance
               and return if estimated step size is too small for machine
               precision. */

            /* restore x, phi and psi */
            *phase1 = false;
            *x = xold;
        }
        for (i = 1; i <= *k; ++i) {
            const int ip1 = i + 1;
            for (l = 1; l <= neqn; ++l) {
                phi[l + i * neqn] =
                    (1.0 / beta[i - 1]) *
                    (phi[l + i * neqn] - phi[l + ip1 * neqn]);
            }
        }
        if (*k >= 2) {
            int i;
            for (i = 2; i <= *k; ++i) {
                psi[i - 2] = psi[i - 1] - *h;
            }
        }
        /* on third failure, set order to one.  thereafter, use optimal step
           size */
        ++ifail;
        {
            double temp2 = .5;
            if (ifail >= 3) {
                if (ifail != 3 && p5eps < erk * .25) {
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
    temp1 = *h * g[kp1];
    if (!*nornd) {
        for (l = 1; l <= neqn; ++l) {
            const double rho =
                temp1 * (yp[l] - phi[l + neqn]) -
                phi[l + (neqn * 16)];
            y[l] = p[l] + rho;
            phi[l + neqn * 15] = y[l] - p[l] - rho;
        }
    } else {
        for (l = 1; l <= neqn; ++l) {
            /* L415: */
            y[l] = p[l] + temp1 * (yp[l] - phi[l + neqn]);
        }
    }
    (*f)(f_ctx, *x, &y[1], &yp[1]);

    /* update differences for next step */

    for (l = 1; l <= neqn; ++l) {
        phi[l + kp1 * neqn] = yp[l] - phi[l + neqn];
        phi[l + kp2 * neqn] = phi[l + kp1 * neqn] - phi[l + kp2 * neqn];
    }
    for (i = 1; i <= *k; ++i) {
        for (l = 1; l <= neqn; ++l) {
            phi[l + i * neqn] += phi[l + kp1 * neqn];
        }
    }

    /* estimate error at order k+1 unless: in first phase when always raise
       order, already decided to lower order, step size not constant so
       estimate unreliable */

    erkp1 = 0.;
    if (knew == km1 || *k == 12) {
        *phase1 = false;
    }

    if (*phase1) {
        *k = kp1;
        erk = erkp1;
    } else if (knew == km1) {
        *k = km1;
        erk = erkm1;
    } else if (kp1 <= *ns) {
        for (l = 1; l <= neqn; ++l) {
            erkp1 += pow(phi[l + kp2 * neqn] / wt[l], 2.0);
        }
        erkp1 = absh * gstr[kp1 - 1] * sqrt(erkp1);

        /* using estimated error at order k+1, determine appropriate order
           for next step */

        if (*k > 1) {
            if (erkm1 <= min(erk, erkp1)) {
                *k = km1;
                erk = erkm1;
            } else if (erkp1 < erk && *k != 12) {
                *k = kp1;
                erk = erkp1;
            }
        } else if (erkp1 < erk * .5) {
            *k = kp1;
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

  yout[] -- Solution at `xout`
  ypout[] -- Derivative of solution at `xout`

  The remaining parameters are returned unaltered from their input values.
  Integration with `step` may be continued.
*/
void intrp(double *const restrict x,
           double *const restrict y,
           const double xout,
           double *const restrict yout,
           double *const restrict ypout,
           const int neqn,
           int *const restrict kold,
           double *const restrict phi,
           double *const restrict psi)
{
    const double hi = xout - *x;
    const int ki = *kold + 1;

    int i, j, l;
    double term = 0.0;
    double g[13] = {1.0};
    double rho[13] = {1.0};
    double w[13] = {0.0};

    /* initialize w for computing g */
    assert(ki <= 13);
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
    clear_double_array(ypout, (size_t)neqn);
    clear_double_array(yout, (size_t)neqn);
    for (i = ki; i-- > 0;) {
        const double temp2 = g[i];
        const double temp3 = rho[i];
        for (l = 0; l < neqn; ++l) {
            yout[l] += temp2 * phi[l + i * neqn];
            ypout[l] += temp3 * phi[l + i * neqn];
        }
    }
    for (l = 0; l < neqn; ++l) {
        yout[l] = y[l] + hi * yout[l];
    }
}

/*
  `ode` merely allocates storage for `de` to relieve the user of the
  inconvenience of a long call list.  Consequently `de` is used as described
  in the comments for `ode` .

  The constant `maxnum` is the maximum number of steps allowed in one call to
  `de`.
*/
void de(const fn_type f,
        void *const restrict f_ctx,
        const int neqn,
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
        double *const restrict alpha,
        double *const restrict beta,
        double *const restrict sig,
        double *const restrict v,
        double *const restrict w,
        double *const restrict g,
        bool *const restrict phase1,
        double *const restrict psi,
        double *const restrict x,
        double *const restrict h,
        double *const restrict hold,
        bool *const restrict start,
        double *const restrict told,
        double *const restrict delsgn,
        int *const restrict ns,
        bool *const restrict nornd,
        int *const restrict k,
        int *const restrict kold,
        int *const restrict isnold,
        const int maxnum)
{
    const int isn = *iflag >= 0 ? 1 : -1;
    const double del = tout - *t;
    const double absdel = fabs(del);
    const double tend = isn < 0 ? tout : *t + del * 10.0;

    bool stiff;
    double abseps, eps, releps;
    int kle4, l, nostep;

    /* test for improper parameters */
    eps = max(*relerr, *abserr);
    *iflag = abs(*iflag);
    if (neqn < 1 || *t == tout ||
        *relerr < 0. || *abserr < 0. ||
        eps <= 0. || *iflag == 0 ||
        (*iflag != 1 && (*t != *told || *iflag < 2 || *iflag > 5))) {
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
    if (*iflag == 1 || *isnold < 0 || *delsgn * del <= 0.) {
        /* on start and restart also set work variables x and yy, store the
           direction of integration and initialize the step size */
        *start = true;
        *x = *t;
        copy_double_array(yy, y, (size_t)neqn);
        *delsgn = copysign(1.0, del);
        *h = copysign(max(fabs(tout - *x), 4.0 * DBL_EPSILON * fabs(*x)),
                      tout - *x);
    }

    /* if already past output point, interpolate and return */
    for (nostep = 0;; ++nostep) {

        if (fabs(*x - *t) >= absdel) {
            intrp(x, yy, tout, y, ypout, neqn, kold, phi, psi);
            *iflag = 2;
            *t = tout;
            *told = *t;
            *isnold = isn;
            return;
        }

        /* if cannot go past output point and sufficiently close, extrapolate and
           return */
        if (isn <= 0 || fabs(tout - *x) < 4.0 * DBL_EPSILON * fabs(*x)) {
            *h = tout - *x;
            (*f)(f_ctx, *x, yy, yp);
            for (l = 0; l < neqn; ++l) {
                y[l] = yy[l] + *h * yp[l];
            }
            *iflag = 2;
            *t = tout;
            *told = *t;
            *isnold = isn;
            return;
        }

        /* test for too many steps */
        if (nostep >= maxnum) {
            *iflag = isn * 4;
            if (stiff) {
                *iflag = isn * 5;
            }
            copy_double_array(y, yy, (size_t)neqn);
            *t = *x;
            *told = *t;
            *isnold = 1;
            return;
        }

        /* limit step size, set weight vector and take a step */
        *h = copysign(min(fabs(*h), fabs(tend - *x)), *h);
        for (l = 0; l < neqn; ++l) {
            wt[l] = releps * fabs(yy[l]) + abseps;
        }

        /* test for tolerances too small */
        if (step(x, yy, f, f_ctx, neqn, h, &eps, wt, start, hold,
                 k, kold, phi, p, yp, psi, alpha, beta,
                 sig, v, w, g, phase1, ns, nornd)) {
            *iflag = isn * 3;
            *relerr = eps * releps;
            *abserr = eps * abseps;
            copy_double_array(y, yy, (size_t)neqn);
            *t = *x;
            *told = *t;
            *isnold = 1;
            return;
        }

        /* augment counter on number of steps and test for stiffness */
        ++kle4;
        if (*kold > 4) {
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
  `step`, and `intrp`.  `ode` allocates virtual storage in the arrays `work`
  and `iwork` and calls `de`.  `de` is a supervisor which directs the
  solution.  It calls on the routines `step` and `intrp` to advance the
  integration and to interpolate at output points.  `step` uses a modified
  divided difference form of the Adams PECE formulas and local extrapolation.
  It adjusts the order and step size to control the local error per unit step
  in a generalized sense.  Normally each call to `step` advances the solution
  one step in the direction of `tout`.  For reasons of efficiency `de`
  integrates beyond `tout` internally, though never beyond
  `t + 10 * (tout - t)`, and calls `intrp` to interpolate the solution at
  `tout`.  An option is provided to stop the integration at `tout` but it
  should be used only if it is impossible to continue the integration beyond
  `tout`.

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
  subsequent calls (length: `100 + 21 * neqn`)

  @param iwork
  Arrays to hold information internal to `de` which is necessary for
  subsequent calls (length: `5`)

  # First call to `ode`

  The user must provide storage in their calling program for the arrays
  in the call list,

      y[neqn], work[100 + 21 * neqn], iwork[5],

  Supply supply the subroutine `f(f_ctx, t, y, yp)` to evaluate

      dy[i]/dt = yp[i] = f(t, y[0], y[1], ..., y[neqn - 1])

  and initialize the parameters:

  - `neqn`: Number of equations to be integrated
  - `y`: Vector of initial conditions
  - `t`: Starting point of integration
  - `tout`: Point at which solution is desired
  - `relerr, abserr`: Relative and absolute local error tolerances
  - `iflag`: `+1` or `-1`. Indicator to initialize the code.  Normal input is
    `+1`.  The user should set `iflag = -1` only if it is impossible to
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
      - `4`: Integration did not reach `tout` because more than `maxnum` steps
         needed.
      - `5`: Integration did not reach `tout` because equations appear to be
        stiff.
      - `6`: Invalid input parameters (fatal error).
  - `work`, `iwork`: Information generally of no interest to the user but
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
         const int neqn,
         double *const restrict y,
         double *const restrict t,
         const double tout,
         double *const restrict relerr,
         double *const restrict abserr,
         int *const restrict iflag,
         double *const restrict work,
         int *const restrict iwork,
         const int maxnum)
{
    static const int ialpha = 0;
    static const int ih = 88;
    static const int ihold = 89;
    static const int istart = 90;
    static const int itold = 91;
    static const int idelsn = 92;
    static const int ibeta = 12;
    static const int isig = 24;
    static const int iv = 37;
    static const int iw = 49;
    static const int ig = 61;
    static const int iphase = 74;
    static const int ipsi = 75;
    static const int ix = 87;
    static const int iyy = 99;

    const int iwt = iyy + neqn;
    const int ip = iwt + neqn;
    const int iyp = ip + neqn;
    const int iypout = iyp + neqn;
    const int iphi = iypout + neqn;

    bool nornd, phase1, start;

    if (abs(*iflag) != 1) {
        start = work[istart] > 0.;
        phase1 = work[iphase] > 0.;
        nornd = iwork[1] != -1;
    }
    de(f, f_ctx, neqn, y, t, tout, relerr, abserr, iflag,
       &work[iyy], &work[iwt], &work[ip], &work[iyp], &work[iypout],
       &work[iphi], &work[ialpha], &work[ibeta], &work[isig], &work[iv],
       &work[iw], &work[ig], &phase1, &work[ipsi], &work[ix], &work[ih],
       &work[ihold], &start, &work[itold], &work[idelsn], &iwork[0],
       &nornd, &iwork[2], &iwork[3], &iwork[4], maxnum);
    if (start) {
        work[istart] = 1.0;
    } else {
        work[istart] = -1.0;
    }
    if (phase1) {
        work[iphase] = 1.0;
    } else {
        work[iphase] = -1.0;
    }
    if (nornd) {
        iwork[1] = 1;
    } else {
        iwork[1] = -1;
    }
}
