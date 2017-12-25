#include <math.h>

const double relerr = 1.0e-8;

const double abserr = 1.0e-8;

const double ts[] = {
    0.0,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    NAN
};

const double y_initial[] = {1.0, 0.0, NAN};

void f(double t, const double *y, double *yp)
{
    (void)t;
    yp[0] = y[1];
    yp[1] = -y[0];
}
