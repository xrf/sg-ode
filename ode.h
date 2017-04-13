#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
/** @file
    ODE solver by Shampine and Gordon.
*/
#include <stdbool.h>
#include <math.h>
#include "vector.h"
#include "extern.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef void SgDerivFn(void *f_ctx,
                       double t,
                       const SgVector *restrict y,
                       SgVector *restrict yp);

struct SgOde {
    struct SgVectorDriver drv;

    /** Solution vector at `x` */
    SgVector *yy;
    /** Derivative of solution vector at `x` after successful step */
    SgVector *yp;
    /** Vector of weights for error criterion */
    SgVector *wt;
    SgVector *p, *ypout, *phi[16];

    double alpha[12], beta[12], sig[13], v[12], w[12], g[13], psi[12];

    /** Independent variable */
    double x;
    /** Appropriate step size for next step.  Normally determined by code */
    double h;
    /** Step size used for last successful step */
    double hold;
    double told, delsgn;
    unsigned ns;
    /** Appropriate order for next step (determined by code).
        Invariant: k >= 1 && k < 13 */
    unsigned k;
    /** Order used for last successful step */
    unsigned kold;
    /** "iflag_sign_old": whether the user-provided iflag was positive
        (controls whether solver is allowed to overshoot and interpolate) */
    bool isnold;
    /** probably means "no_round" and has something to do with rounding */
    bool nornd;
    bool phase1;
    /** `true` on first step, `false` otherwise */
    bool start;
};

void sg_ode_init(struct SgOde *self, struct SgVectorDriver drv);

void sg_ode_del(struct SgOde *self);

/**
  Integrates a system of first order ordinary differential equations one step,
  normally from `x` to `x+h`, using a modified divided difference form of the
  Adams PECE formulas.  Local extrapolation is used to improve absolute
  stability and accuracy.  The code adjusts its order and step size to control
  the local error per unit step in a generalized sense.  Special devices are
  included to control roundoff error and to detect when the user is requesting
  too much accuracy.

  @param[in,out] self   The solver state.
  @param[in]     f      The derivative function.
  @param[in,out] f_ctx  Arbitrary pointer for `f`.
  @param[in,out] eps    Local error tolerance. (length: `1`)

  @return
  Nonzero when no step can be taken, zero otherwise (`crash`).

  The arrays `phi`, `psi` are required for the interpolation function
  `interpolate`.  The array `p` is internal to the code.

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
  endpoint and the solution interpolated there with `interpolate`.  If it is
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
SG_EXTERN int sg_ode_step(struct SgOde *self,
                          SgDerivFn *f,
                          void *restrict f_ctx,
                          double *restrict eps);

/**
  The methods in function `step` approximate the solution near `x` by a
  polynomial.  Function `interpolate` approximates the solution at `xout` by
  evaluating the polynomial there.  Information defining this polynomial is
  passed from `step` so `interpolate` cannot be used alone.

  # Input to `interpolate`

  The user provides storage in the calling program for the arrays in the call
  list and defines `xout`, point at which solution is desired.

  The remaining parameters are defined in `step` and passed to `interpolate`
  from that function.

  # Output from `interpolate`

    - `yout`: Solution at `xout`
    - `ypout`: Derivative of solution at `xout`

  The remaining parameters are returned unaltered from their input values.
  Integration with `step` may be continued.
*/
SG_EXTERN void sg_ode_interpolate(struct SgVectorDriver drv,
                                  double x,
                                  const SgVector *restrict y,
                                  double xout,
                                  SgVector *restrict yout,
                                  SgVector *restrict ypout,
                                  unsigned kold,
                                  SgVector *const restrict *phi,
                                  const double *restrict psi);

/**
  Integrates a system of `neqn` first order ordinary differential equations.

  The equations are of the form:

      dy[i]/dt = f(t, y[0], y[1], ..., y[neqn - 1])    y[i] given at t

  The function integrates from `t` to `tout`.  On return the parameters in the
  call list are set for continuing the integration.  The user has only to
  define a new value `tout` and call `ode` again.

  The differential equations are actually solved by a suite of codes `de`,
  `step`, and `interpolate`.  `ode` allocates virtual storage in the `work`
  array and calls `de`.  `de` is a supervisor which directs the solution.  It
  calls on the routines `step` and `interpolate` to advance the integration
  and to interpolate at output points.  `step` uses a modified divided
  difference form of the Adams PECE formulas and local extrapolation.  It
  adjusts the order and step size to control the local error per unit step in
  a generalized sense.  Normally each call to `step` advances the solution one
  step in the direction of `tout`.  For reasons of efficiency `de` integrates
  beyond `tout` internally, though never beyond `t + 10 * (tout - t)`, and
  calls `interpolate` to interpolate the solution at `tout`.  An option is
  provided to stop the integration at `tout` but it should be used only if it
  is impossible to continue the integration beyond `tout`.

  This code is completely explained and documented in the text, Computer
  Solution of Ordinary Differential Equations: The Initial Value Problem
  L. F. Shampine and M. K. Gordon.

  @param f
  Function `f(f_ctx, t, y, yp)` to evaluate derivatives `yp[i] = dy[i]/dt`

  @param f_ctx
  Arbitrary pointer passed to `f`.

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

  @param maxnum
  Maximum number of steps allowed in one call to `de`.

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
SG_EXTERN void sg_ode_de(struct SgOde *self,
                         SgDerivFn *f,
                         void *restrict f_ctx,
                         SgVector *restrict y,
                         double *restrict t,
                         double tout,
                         double *restrict relerr,
                         double *restrict abserr,
                         unsigned maxnum,
                         int *restrict iflag);

#ifdef __cplusplus
}
#endif
#endif
