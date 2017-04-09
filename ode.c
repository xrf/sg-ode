// origin: http://www.netlib.org/ode/ode.f
// f2c'ed

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "ode.h"

typedef void (*func_type)(void *, double, const double *, double *);

static double min(double x, double y)
{
    return x < y ? x : y;
}

static double max(double x, double y)
{
    return x >= y ? x : y;
}

/*
  Integrates a system of first order ordinary differential equations one step,
  normally from `x` to `x+h`, using a modified divided difference form of the
  Adams PECE formulas.  Local extrapolation is used to improve absolute
  stability and accuracy.  The code adjusts its order and step size to control
  the local error per unit step in a generalized sense.  Special devices are
  included to control roundoff error and to detect when the user is requesting
  too much accuracy.

  @param x Independent variable
  @param y Solution vector at `x` (length: `neqn`)
  @param yp Derivative of solution vector at `x` after successful
            step (length: `neqn`)
  @param neqn Number of equations to be integrated
  @param h Appropriate step size for next step.  Normally determined by
           code
  @param eps Local error tolerance.  Must be variable (length: `1`)
  @param wt Vector of weights for error criterion (length: `neqn`)
  @param start Boolean variable set `true` for first step, `false`
               otherwise
  @param hold Step size used for last successful step
  @param k Appropriate order for next step (determined by code)
  @param kold Order used for last successful step
  @param crash Boolean variable set `true` when no step can be taken,
               `false` otherwise.

  The arrays `phi`, `psi` are required for the interpolation subroutine
  `intrp`.  The array `p` is internal to the code.

  # Input to `step`

  ## First call

  The user must provide storage in his driver program for all arrays in the
  call list, namely

      y[neqn], wt[neqn], phi[neqn * 16], p[neqn], yp[neqn], psi[12]

  The user must also declare `start` and `crash` Boolean variables and `f` an
  external subroutine, supply the subroutine `f(f_ctx, x, y, yp)` to evaluate

      dy[i]/dx = yp[i] = f(x, y[0], y[1], ..., y[neqn - 1])

  and initialize only the following parameters:

      x -- initial value of the independent variable
      y[] -- vector of initial values of dependent variables
      neqn -- number of equations to be integrated
      h -- nominal step size indicating direction of integration
           and maximum size of step.  must be variable
      eps -- local error tolerance per step.  must be variable
      wt[] -- vector of non-zero weights for error criterion
      start -- .true.

  `step` requires the L2 norm of the vector with components
  `local_error[l] / wt[l]` be less than `eps` for a successful step.  The
  array `wt` allows the user to specify an error test appropriate for his
  problem.  For example,

      wt[l] = 1.0  specifies absolute error,
            = fabs(y[l])  error relative to the most recent value of
                          the l-th component of the solution,
            = fabs(yp[l])  error relative to the most recent value of
                           the l-th component of the derivative,
            = max(wt[l], fabs(y[l]))  error relative to the largest magnitude
                                      of l-th component obtained so far,
            = fabs(y(l)) * relerr / eps + abserr / eps
                 specifies a mixed relative-absolute test where relerr is
                 relative error, abserr is absolute error and
                 eps = max(relerr, abserr) .

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

  The subroutine returns after each successful step with `start` and `crash`
  set to `false`.  `x` represents the independent variable advanced one step
  of length hold from its value on input and `y` the solution vector at the
  new value of `x`.  All other parameters represent information corresponding
  to the new `x` needed to continue the integration.

  ## unsuccessful step

  When the error tolerance is too small for the machine precision, the
  subroutine returns without taking a step and `crash = true`.  An appropriate
  step size and error tolerance for continuing are estimated and all other
  information is restored as upon input before returning.  To continue with
  the larger tolerance, the user just calls the code again.  A restart is
  neither required nor desirable.
*/
int step(double *x, double *y, func_type f, void *f_ctx, int neqn,
         double *h__, double *eps, double *wt, bool *start,
         double *hold, int *k, int *kold, bool *crash,
         double *phi, double *p, double *yp, double *psi,
         double *alpha, double *beta, double *sig, double *v,
         double *w, double *g, bool *phase1, int *ns, bool *nornd)
{
    static const double gstr[13] = {.5, .0833, .0417, .0264, .0188, .0143, .0114,
                                    .00936, .00789, .00679, .00592, .00524, .00468};

    /* System generated locals */
    int phi_dim1, phi_offset, i__1, i__2;
    double d__1, d__2, d__3;

    /* Local variables */
    static int i__, j, l;
    static double r__;
    static int iq, im1, km1, km2, ip1, kp1, kp2;
    static double erk, err, tau, rho, sum;
    static int nsm2, nsp1, nsp2;
    static double absh, hnew;
    static int knew;
    static double xold, twou, erkm1, erkm2, erkp1, temp1, temp2, temp3,
        temp4, temp5, temp6, p5eps;
    static int ifail;
    static double reali, round;
    static double fouru;
    static int limit1, limit2;
    static double realns;

    /* Parameter adjustments */
    --yp;
    --p;
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --wt;
    --y;
    --psi;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;

    /* Function Body */
    /*     data g(1),g(2)/1.0,0.5/,sig(1)/1.0/ */

    twou = 2.f * DBL_EPSILON;
    fouru = twou * 2.f;
    /*       ***     begin block 0     *** */
    /*   check if step size or error tolerance is too small for machine */
    /*   precision.  if first step, initialize phi array and estimate a */
    /*   starting step size. */
    /*                   *** */

    /*   if step size is too small, determine an acceptable one */

    *crash = true;
    if (fabs(*h__) < fouru * fabs(*x)) {
        d__1 = fouru * fabs(*x);
        *h__ = copysign(d__1, *h__);
        return 0;
    }
    p5eps = *eps * .5;

    /*   if error tolerance is too small, increase it to an acceptable value */

    round = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L10: */
        /* Computing 2nd power */
        d__1 = y[l] / wt[l];
        round += d__1 * d__1;
    }
    round = twou * sqrt(round);
    if (p5eps >= round) {
        goto L15;
    }
    *eps = round * 2.f * (fouru + 1.);
    return 0;
L15:
    *crash = false;
    g[1] = 1.;
    g[2] = .5;
    sig[1] = 1.;
    if (!(*start)) {
        goto L99;
    }

    /*   initialize.  compute appropriate step size for first step */

    (*f)(f_ctx, *x, &y[1], &yp[1]);
    sum = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        phi[l + phi_dim1] = yp[l];
        phi[l + (phi_dim1 << 1)] = 0.;
        /* L20: */
        /* Computing 2nd power */
        d__1 = yp[l] / wt[l];
        sum += d__1 * d__1;
    }
    sum = sqrt(sum);
    absh = fabs(*h__);
    if (*eps < sum * 16. * *h__ * *h__) {
        absh = sqrt(*eps / sum) * .25;
    }
    /* Computing MAX */
    d__2 = absh, d__3 = fouru * fabs(*x);
    d__1 = max(d__2, d__3);
    *h__ = copysign(d__1, *h__);
    *hold = 0.;
    *k = 1;
    *kold = 0;
    *start = false;
    *phase1 = true;
    *nornd = true;
    if (p5eps > round * 100.) {
        goto L99;
    }
    *nornd = false;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L25: */
        phi[l + phi_dim1 * 15] = 0.;
    }
L99:
    ifail = 0;
/*       ***     end block 0     *** */

/*       ***     begin block 1     *** */
/*   compute coefficients of formulas for this step.  avoid computing */
/*   those quantities not changed when step size is not changed. */
/*                   *** */

L100:
    kp1 = *k + 1;
    kp2 = *k + 2;
    km1 = *k - 1;
    km2 = *k - 2;

    /*   ns is the number of steps taken with size h, including the current */
    /*   one.  when k.lt.ns, no coefficients change */

    if (*h__ != *hold) {
        *ns = 0;
    }
    if (*ns <= *kold) {
        ++(*ns);
    }
    nsp1 = *ns + 1;
    if (*k < *ns) {
        goto L199;
    }

    /*   compute those components of alpha(*),beta(*),psi(*),sig(*) which */
    /*   are changed */

    beta[*ns] = 1.;
    realns = (double)(*ns);
    alpha[*ns] = 1. / realns;
    temp1 = *h__ * realns;
    sig[nsp1] = 1.;
    if (*k < nsp1) {
        goto L110;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
        im1 = i__ - 1;
        temp2 = psi[im1];
        psi[im1] = temp1;
        beta[i__] = beta[im1] * psi[im1] / temp2;
        temp1 = temp2 + *h__;
        alpha[i__] = *h__ / temp1;
        reali = (double)i__;
        /* L105: */
        sig[i__ + 1] = reali * alpha[i__] * sig[i__];
    }
L110:
    psi[*k] = temp1;

    /*   compute coefficients g(*) */

    /*   initialize v(*) and set w(*).  g(2) is set in data statement */

    if (*ns > 1) {
        goto L120;
    }
    i__1 = *k;
    for (iq = 1; iq <= i__1; ++iq) {
        temp3 = (double)(iq * (iq + 1));
        v[iq] = 1. / temp3;
        /* L115: */
        w[iq] = v[iq];
    }
    goto L140;

/*   if order was raised, update diagonal part of v(*) */

L120:
    if (*k <= *kold) {
        goto L130;
    }
    temp4 = (double)(*k * kp1);
    v[*k] = 1. / temp4;
    nsm2 = *ns - 2;
    if (nsm2 < 1) {
        goto L130;
    }
    i__1 = nsm2;
    for (j = 1; j <= i__1; ++j) {
        i__ = *k - j;
        /* L125: */
        v[i__] -= alpha[j + 1] * v[i__ + 1];
    }

/*   update v(*) and set w(*) */

L130:
    limit1 = kp1 - *ns;
    temp5 = alpha[*ns];
    i__1 = limit1;
    for (iq = 1; iq <= i__1; ++iq) {
        v[iq] -= temp5 * v[iq + 1];
        /* L135: */
        w[iq] = v[iq];
    }
    g[nsp1] = w[1];

/*   compute the g(*) in the work vector w(*) */

L140:
    nsp2 = *ns + 2;
    if (kp1 < nsp2) {
        goto L199;
    }
    i__1 = kp1;
    for (i__ = nsp2; i__ <= i__1; ++i__) {
        limit2 = kp2 - i__;
        temp6 = alpha[i__ - 1];
        i__2 = limit2;
        for (iq = 1; iq <= i__2; ++iq) {
            /* L145: */
            w[iq] -= temp6 * w[iq + 1];
        }
        /* L150: */
        g[i__] = w[1];
    }
L199:
    /*       ***     end block 1     *** */

    /*       ***     begin block 2     *** */
    /*   predict a solution p(*), evaluate derivatives using predicted */
    /*   solution, estimate local error at order k and errors at orders k, */
    /*   k-1, k-2 as if constant step size were used. */
    /*                   *** */

    /*   change phi to phi star */

    if (*k < nsp1) {
        goto L215;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
        temp1 = beta[i__];
        i__2 = neqn;
        for (l = 1; l <= i__2; ++l) {
            /* L205: */
            phi[l + i__ * phi_dim1] = temp1 * phi[l + i__ * phi_dim1];
        }
        /* L210: */
    }

/*   predict solution and differences */

L215:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1];
        phi[l + kp1 * phi_dim1] = 0.;
        /* L220: */
        p[l] = 0.;
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        i__ = kp1 - j;
        ip1 = i__ + 1;
        temp2 = g[i__];
        i__2 = neqn;
        for (l = 1; l <= i__2; ++l) {
            p[l] += temp2 * phi[l + i__ * phi_dim1];
            /* L225: */
            phi[l + i__ * phi_dim1] += phi[l + ip1 * phi_dim1];
        }
        /* L230: */
    }
    if (*nornd) {
        goto L240;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        tau = *h__ * p[l] - phi[l + phi_dim1 * 15];
        p[l] = y[l] + tau;
        /* L235: */
        phi[l + (phi_dim1 << 4)] = p[l] - y[l] - tau;
    }
    goto L250;
L240:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L245: */
        p[l] = y[l] + *h__ * p[l];
    }
L250:
    xold = *x;
    *x += *h__;
    absh = fabs(*h__);
    (*f)(f_ctx, *x, &p[1], &yp[1]);

    /*   estimate errors at orders k,k-1,k-2 */

    erkm2 = 0.;
    erkm1 = 0.;
    erk = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        temp3 = 1. / wt[l];
        temp4 = yp[l] - phi[l + phi_dim1];
        if (km2 < 0) {
            goto L265;
        } else if (km2 == 0) {
            goto L260;
        } else {
            goto L255;
        }
    L255:
        /* Computing 2nd power */
        d__1 = (phi[l + km1 * phi_dim1] + temp4) * temp3;
        erkm2 += d__1 * d__1;
    L260:
        /* Computing 2nd power */
        d__1 = (phi[l + *k * phi_dim1] + temp4) * temp3;
        erkm1 += d__1 * d__1;
    L265:
        /* Computing 2nd power */
        d__1 = temp4 * temp3;
        erk += d__1 * d__1;
    }
    if (km2 < 0) {
        goto L280;
    } else if (km2 == 0) {
        goto L275;
    } else {
        goto L270;
    }
L270:
    erkm2 = absh * sig[km1] * gstr[km2 - 1] * sqrt(erkm2);
L275:
    erkm1 = absh * sig[*k] * gstr[km1 - 1] * sqrt(erkm1);
L280:
    temp5 = absh * sqrt(erk);
    err = temp5 * (g[*k] - g[kp1]);
    erk = temp5 * sig[kp1] * gstr[*k - 1];
    knew = *k;

    /*   test if order should be lowered */

    if (km2 < 0) {
        goto L299;
    } else if (km2 == 0) {
        goto L290;
    } else {
        goto L285;
    }
L285:
    if (max(erkm1, erkm2) <= erk) {
        knew = km1;
    }
    goto L299;
L290:
    if (erkm1 <= erk * .5) {
        knew = km1;
    }

/*   test if step successful */

L299:
    if (err <= *eps) {
        goto L400;
    }
    /*       ***     end block 2     *** */

    /*       ***     begin block 3     *** */
    /*   the step is unsuccessful.  restore  x, phi(*,*), psi(*) . */
    /*   if third consecutive failure, set order to one.  if step fails more */
    /*   than three times, consider an optimal step size.  double error */
    /*   tolerance and return if estimated step size is too small for machine */
    /*   precision. */
    /*                   *** */

    /*   restore x, phi(*,*) and psi(*) */

    *phase1 = false;
    *x = xold;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        temp1 = 1. / beta[i__];
        ip1 = i__ + 1;
        i__2 = neqn;
        for (l = 1; l <= i__2; ++l) {
            /* L305: */
            phi[l + i__ * phi_dim1] = temp1 * (phi[l + i__ * phi_dim1] - phi[l + ip1 * phi_dim1]);
        }
        /* L310: */
    }
    if (*k < 2) {
        goto L320;
    }
    i__1 = *k;
    for (i__ = 2; i__ <= i__1; ++i__) {
        /* L315: */
        psi[i__ - 1] = psi[i__] - *h__;
    }

/*   on third failure, set order to one.  thereafter, use optimal step */
/*   size */

L320:
    ++ifail;
    temp2 = .5;
    if ((i__1 = ifail - 3) < 0) {
        goto L335;
    } else if (i__1 == 0) {
        goto L330;
    } else {
        goto L325;
    }
L325:
    if (p5eps < erk * .25) {
        temp2 = sqrt(p5eps / erk);
    }
L330:
    knew = 1;
L335:
    *h__ = temp2 * *h__;
    *k = knew;
    if (fabs(*h__) >= fouru * fabs(*x)) {
        goto L340;
    }
    *crash = true;
    d__1 = fouru * fabs(*x);
    *h__ = copysign(d__1, *h__);
    *eps += *eps;
    return 0;
L340:
    goto L100;
/*       ***     end block 3     *** */

/*       ***     begin block 4     *** */
/*   the step is successful.  correct the predicted solution, evaluate */
/*   the derivatives using the corrected solution and update the */
/*   differences.  determine best order and step size for next step. */
/*                   *** */
L400:
    *kold = *k;
    *hold = *h__;

    /*   correct and evaluate */

    temp1 = *h__ * g[kp1];
    if (*nornd) {
        goto L410;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        rho = temp1 * (yp[l] - phi[l + phi_dim1]) - phi[l + (phi_dim1 << 4)];
        y[l] = p[l] + rho;
        /* L405: */
        phi[l + phi_dim1 * 15] = y[l] - p[l] - rho;
    }
    goto L420;
L410:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L415: */
        y[l] = p[l] + temp1 * (yp[l] - phi[l + phi_dim1]);
    }
L420:
    (*f)(f_ctx, *x, &y[1], &yp[1]);

    /*   update differences for next step */

    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        phi[l + kp1 * phi_dim1] = yp[l] - phi[l + phi_dim1];
        /* L425: */
        phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1] - phi[l + kp2 * phi_dim1];
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        i__2 = neqn;
        for (l = 1; l <= i__2; ++l) {
            /* L430: */
            phi[l + i__ * phi_dim1] += phi[l + kp1 * phi_dim1];
        }
        /* L435: */
    }

    /*   estimate error at order k+1 unless: */
    /*     in first phase when always raise order, */
    /*     already decided to lower order, */
    /*     step size not constant so estimate unreliable */

    erkp1 = 0.;
    if (knew == km1 || *k == 12) {
        *phase1 = false;
    }
    if (*phase1) {
        goto L450;
    }
    if (knew == km1) {
        goto L455;
    }
    if (kp1 > *ns) {
        goto L460;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L440: */
        /* Computing 2nd power */
        d__1 = phi[l + kp2 * phi_dim1] / wt[l];
        erkp1 += d__1 * d__1;
    }
    erkp1 = absh * gstr[kp1 - 1] * sqrt(erkp1);

    /*   using estimated error at order k+1, determine appropriate order */
    /*   for next step */

    if (*k > 1) {
        goto L445;
    }
    if (erkp1 >= erk * .5) {
        goto L460;
    }
    goto L450;
L445:
    if (erkm1 <= min(erk, erkp1)) {
        goto L455;
    }
    if (erkp1 >= erk || *k == 12) {
        goto L460;
    }

/*   here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have */
/*   been lowered in block 2.  thus order is to be raised */

/*   raise order */

L450:
    *k = kp1;
    erk = erkp1;
    goto L460;

/*   lower order */

L455:
    *k = km1;
    erk = erkm1;

/*   with new order determine appropriate step size for next step */

L460:
    hnew = *h__ + *h__;
    if (*phase1) {
        goto L465;
    }
    if (p5eps >= erk * pow(2.0, *k + 1)) {
        goto L465;
    }
    hnew = *h__;
    if (p5eps >= erk) {
        goto L465;
    }
    temp2 = (double)(*k + 1);
    d__1 = p5eps / erk;
    d__2 = 1. / temp2;
    r__ = pow(d__1, d__2);
    /* Computing MAX */
    d__1 = .5, d__2 = min(.9, r__);
    hnew = absh * max(d__1, d__2);
    /* Computing MAX */
    d__2 = hnew, d__3 = fouru * fabs(*x);
    d__1 = max(d__2, d__3);
    hnew = copysign(d__1, *h__);
L465:
    *h__ = hnew;
    return 0;
    /*       ***     end block 4     *** */
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
int intrp(double *x, double *y, double xout,
          double *yout, double *ypout, int neqn, int *kold,
          double *phi, double *psi)
{
    /* Initialized data */

    static double g[13] = {1.};
    static double rho[13] = {1.};

    /* System generated locals */
    int phi_dim1, phi_offset, i__1, i__2;

    /* Local variables */
    static int i__, j, l;
    static double w[13], hi;
    static int ki, jm1;
    static double eta;
    static int kip1;
    static double term, temp1, temp2, temp3, gamma;
    static int limit1;
    static double psijm1;

    /* Parameter adjustments */
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --ypout;
    --yout;
    --y;
    --psi;

    /* Function Body */

    hi = xout - *x;
    ki = *kold + 1;
    kip1 = ki + 1;

    /*   initialize w(*) for computing g(*) */

    i__1 = ki;
    for (i__ = 1; i__ <= i__1; ++i__) {
        temp1 = (double)i__;
        /* L5: */
        w[i__ - 1] = 1. / temp1;
    }
    term = 0.;

    /*   compute g(*) */

    i__1 = ki;
    for (j = 2; j <= i__1; ++j) {
        jm1 = j - 1;
        psijm1 = psi[jm1];
        gamma = (hi + term) / psijm1;
        eta = hi / psijm1;
        limit1 = kip1 - j;
        i__2 = limit1;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L10: */
            w[i__ - 1] = gamma * w[i__ - 1] - eta * w[i__];
        }
        g[j - 1] = w[0];
        rho[j - 1] = gamma * rho[jm1 - 1];
        /* L15: */
        term = psijm1;
    }

    /*   interpolate */

    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        ypout[l] = 0.;
        /* L20: */
        yout[l] = 0.;
    }
    i__1 = ki;
    for (j = 1; j <= i__1; ++j) {
        i__ = kip1 - j;
        temp2 = g[i__ - 1];
        temp3 = rho[i__ - 1];
        i__2 = neqn;
        for (l = 1; l <= i__2; ++l) {
            yout[l] += temp2 * phi[l + i__ * phi_dim1];
            /* L25: */
            ypout[l] += temp3 * phi[l + i__ * phi_dim1];
        }
        /* L30: */
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L35: */
        yout[l] = y[l] + hi * yout[l];
    }
    return 0;
}

/*
  `ode` merely allocates storage for `de` to relieve the user of the
  inconvenience of a long call list.  Consequently `de` is used as described
  in the comments for `ode` .

  The constant `maxnum` is the maximum number of steps allowed in one call to
  `de`.
*/
int de(func_type f, void *f_ctx, int neqn, double *y, double *t,
       double tout, double *relerr, double *abserr, int *iflag,
       double *yy, double *wt, double *p, double *yp,
       double *ypout, double *phi, double *alpha, double *beta,
       double *sig, double *v, double *w, double *g,
       bool *phase1, double *psi, double *x, double *h__,
       double *hold, bool *start, double *told, double *delsgn,
       int *ns, bool *nornd, int *k, int *kold, int *isnold, int maxnum)
{
    /* System generated locals */
    int phi_dim1, phi_offset, i__1;
    double d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int l;
    static double del, eps;
    static int isn, kle4;
    static double tend;
    static bool crash, stiff;
    static double fouru, absdel, abseps, releps;
    static int nostep;

    /* Parameter adjustments */
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --ypout;
    --yp;
    --p;
    --wt;
    --yy;
    --y;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;
    --psi;

    /* Function Body */

    /*            ***            ***            *** */
    /*   test for improper parameters */

    fouru = 4.f * DBL_EPSILON;
    if (neqn < 1) {
        goto L10;
    }
    if (*t == tout) {
        goto L10;
    }
    if (*relerr < 0. || *abserr < 0.) {
        goto L10;
    }
    eps = max(*relerr, *abserr);
    if (eps <= 0.) {
        goto L10;
    }
    if (*iflag == 0) {
        goto L10;
    }
    if (*iflag >= 0) {
        isn = 1;
    } else {
        isn = -1;
    }
    *iflag = abs(*iflag);
    if (*iflag == 1) {
        goto L20;
    }
    if (*t != *told) {
        goto L10;
    }
    if (*iflag >= 2 && *iflag <= 5) {
        goto L20;
    }
L10:
    *iflag = 6;
    return 0;

/*   on each call set interval of integration and counter for number of */
/*   steps.  adjust input error tolerances to define weight vector for */
/*   subroutine  step */

L20:
    del = tout - *t;
    absdel = fabs(del);
    tend = *t + del * 10.;
    if (isn < 0) {
        tend = tout;
    }
    nostep = 0;
    kle4 = 0;
    stiff = false;
    releps = *relerr / eps;
    abseps = *abserr / eps;
    if (*iflag == 1) {
        goto L30;
    }
    if (*isnold < 0) {
        goto L30;
    }
    if (*delsgn * del > 0.) {
        goto L50;
    }

/*   on start and restart also set work variables x and yy(*), store the */
/*   direction of integration and initialize the step size */

L30:
    *start = true;
    *x = *t;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L40: */
        yy[l] = y[l];
    }
    *delsgn = copysign(1.0, del);
    /* Computing MAX */
    d__3 = (d__1 = tout - *x, fabs(d__1)), d__4 = fouru * fabs(*x);
    d__2 = max(d__3, d__4);
    d__5 = tout - *x;
    *h__ = copysign(d__2, d__5);

/*   if already past output point, interpolate and return */

L50:
    if ((d__1 = *x - *t, fabs(d__1)) < absdel) {
        goto L60;
    }
    intrp(x, &yy[1], tout, &y[1], &ypout[1], neqn, kold, &phi[phi_offset], &psi[1]);
    *iflag = 2;
    *t = tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   if cannot go past output point and sufficiently close, */
/*   extrapolate and return */

L60:
    if (isn > 0 || (d__1 = tout - *x, fabs(d__1)) >= fouru * fabs(*x)) {
        goto L80;
    }
    *h__ = tout - *x;
    (*f)(f_ctx, *x, &yy[1], &yp[1]);
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L70: */
        y[l] = yy[l] + *h__ * yp[l];
    }
    *iflag = 2;
    *t = tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   test for too many steps */

L80:
    if (nostep < maxnum) {
        goto L100;
    }
    *iflag = isn << 2;
    if (stiff) {
        *iflag = isn * 5;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L90: */
        y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   limit step size, set weight vector and take a step */

L100:
    /* Computing MIN */
    d__3 = fabs(*h__), d__4 = (d__1 = tend - *x, fabs(d__1));
    d__2 = min(d__3, d__4);
    *h__ = copysign(d__2, *h__);
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L110: */
        wt[l] = releps * (d__1 = yy[l], fabs(d__1)) + abseps;
    }
    step(x, &yy[1], f, f_ctx, neqn, h__, &eps, &wt[1], start, hold, k, kold, &crash, &phi[phi_offset], &p[1], &yp[1], &psi[1], &alpha[1], &beta[1], &sig[1], &v[1], &w[1], &g[1], phase1, ns, nornd);

    /*   test for tolerances too small */

    if (!crash) {
        goto L130;
    }
    *iflag = isn * 3;
    *relerr = eps * releps;
    *abserr = eps * abseps;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
        /* L120: */
        y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   augment counter on number of steps and test for stiffness */

L130:
    ++nostep;
    ++kle4;
    if (*kold > 4) {
        kle4 = 0;
    }
    if (kle4 >= 50) {
        stiff = true;
    }
    goto L50;
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

  The user must provide storage in his calling program for the arrays
  in the call list,

      y[neqn], work[100 + 21 * neqn], iwork[5],

  Supply supply the subroutine `f(f_ctx, t, y, yp)` to evaluate

      dy[i]/dt = yp[i] = f(t, y[0], y[1], ..., y[neqn - 1])

  and initialize the parameters:

      `neqn` -- number of equations to be integrated
      `y[]` -- vector of initial conditions
      `t` -- starting point of integration
      `tout` -- point at which solution is desired
      `relerr, abserr` -- relative and absolute local error tolerances
      `iflag` -- +1,-1.  indicator to initialize the code.  normal input
           is +1.  the user should set iflag=-1 only if it is
           impossible to continue the integration beyond `tout`.

  All parameters except `f`, `neqn` and `tout`  may be altered by the
  code on output so must be variables in the calling program.

  # Output from `ode`

      `neqn` -- unchanged
      `y[]` -- solution at `t`
      `t` -- last point reached in integration.  normal return has
           `t == tout`.
      `tout` -- unchanged
      `relerr`, `abserr` -- normal return has tolerances unchanged.  `iflag=3`
           signals tolerances increased
      iflag = 2 -- normal return.  integration reached  tout
            = 3 -- integration did not reach  tout  because error
                   tolerances too small.  relerr ,  abserr  increased
                   appropriately for continuing
            = 4 -- integration did not reach  tout  because more than
                   500 steps needed
            = 5 -- integration did not reach  tout  because equations
                   appear to be stiff
            = 6 -- invalid input parameters (fatal error)
           the value of  iflag  is returned negative when the input
           value is negative and the integration does not reach  tout ,
           i.e., -3, -4, -5.
      work[], iwork[] -- information generally of no interest to the
           user but necessary for subsequent calls.

  # Subsequent calls to `ode`

  Subroutine `ode` returns with all information needed to continue the
  integration.  If the integration reached `tout`, the user need only define a
  new `tout` and call again.  If the integration did not reach `tout` and the
  user wants to continue, he just calls again.  The output value of `iflag` is
  the appropriate input value for subsequent calls.  The only situation in
  which it should be altered is to stop the integration internally at the new
  `tout`, i.e., change output `iflag=2` to input `iflag=-2`.  Error tolerances
  may be changed by the user before continuing.  All other parameters must
  remain unchanged.
*/
int ode(const func_type f,
        void *const f_ctx,
        const int neqn,
        double *y,
        double *const t,
        const double tout,
        double *const relerr,
        double *const abserr,
        int *const iflag,
        double *work,
        int *iwork,
        const int maxnum)
{
    static const int ialpha = 1;
    static const int ih = 89;
    static const int ihold = 90;
    static const int istart = 91;
    static const int itold = 92;
    static const int idelsn = 93;
    static const int ibeta = 13;
    static const int isig = 25;
    static const int iv = 38;
    static const int iw = 50;
    static const int ig = 62;
    static const int iphase = 75;
    static const int ipsi = 76;
    static const int ix = 88;

    static int ip, iyp, iwt, iyy, iphi;
    static bool nornd, start, phase1;
    static int iypout;

    /* Parameter adjustments */
    --y;
    --work;
    --iwork;

    /* Function Body */
    iyy = 100;
    iwt = iyy + neqn;
    ip = iwt + neqn;
    iyp = ip + neqn;
    iypout = iyp + neqn;
    iphi = iypout + neqn;
    if (abs(*iflag) == 1) {
        goto L1;
    }
    start = work[istart] > 0.;
    phase1 = work[iphase] > 0.;
    nornd = iwork[2] != -1;
L1:
    de(f, f_ctx, neqn, &y[1], t, tout, relerr, abserr, iflag, &work[iyy], &work[iwt], &work[ip], &work[iyp], &work[iypout], &work[iphi], &work[ialpha], &work[ibeta], &work[isig], &work[iv], &work[iw], &work[ig], &phase1, &work[ipsi], &work[ix], &work[ih], &work[ihold], &start, &work[itold], &work[idelsn], &iwork[1], &nornd, &iwork[3], &iwork[4], &iwork[5], maxnum);
    work[istart] = -1.;
    if (start) {
        work[istart] = 1.;
    }
    work[iphase] = -1.;
    if (phase1) {
        work[iphase] = 1.;
    }
    iwork[2] = -1;
    if (nornd) {
        iwork[2] = 1;
    }
    return 0;
}
