/* the original test from Netlib
   http://www.netlib.org/ode/ex/ode.f */
#include <math.h>
#include "jacobian_elliptic.inl"

const double relerr = 1.0e-09;

const double abserr = 1.0e-16;

const double ts[] = {
    0.0,
    5.0,
    10.0,
    15.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
    45.0,
    50.0,
    55.0,
    60.0,
    NAN
};
