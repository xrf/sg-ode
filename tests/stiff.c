#include <math.h>

const double relerr = 1.0e-8;

const double abserr = 1.0e-8;

const double ts[] = {
    0.0,
    1.0,
    NAN
};

const double y_initial[] = {1.0, 0.0, NAN};

void f(double t, const double *y, double *yp)
{
    (void)t;
    yp[0] = y[1];
    yp[1] = -1000.0 * y[0] - 1001.0 * y[1];
}
