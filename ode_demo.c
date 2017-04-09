// Demonstration program for ode/de/step package.
//
// reference:shampine and gordon's computer solution of
// ordinary differential equations:the initial value problem.
//
// the package is used to solve the defining equations for the
// jacobian elliptic functions:
//
//     y1' = y2*y3,      y1(0) = 0
//     y2' = -y1*y3,     y2(0) = 1
//     y3' = -ksq*y1*y2, y3(0) = 1
//
// which is solved for ti=i*5 for i=0, 1, ..., 12.
// the analytic solution is
//
//     y1 = sn(t/ksq)
//     y2 = cn(t/ksq)
//     y3 = dn(t/ksq)
//
// in this case we use ksq=k*k=0.51.

#include <stdio.h>
#include "ode.h"

void f(void *ctx, double t, const double *y, double *yp)
{
    static const double ksq = 0.51;
    (void)ctx;
    (void)t;
    yp[0] = y[1] * y[2];
    yp[1] = -y[0] * y[2];
    yp[2] = -ksq * y[0] * y[1];
}

int main(void)
{
    static const int maxnum = 500;
    static const size_t neqn = 3;
    double y[3], t, tout, work[63] = {0.0};
    int i;
    struct Ode iwork = ODE_INITIALIZER;
    double relerr = 1.0e-9;
    double abserr = 1.0e-16;
    int iflag = 1;
    printf("neqn=%zu relerr=%g abserr=%g iflag=%i\n",
           neqn, relerr, abserr, iflag);

#define dump() \
    printf("t=%2.0f y=[%17.10e %17.10e %17.10e] iflag=%i\n", \
           t, y[0], y[1], y[2], iflag);

    t = 0.0;
    y[0] = 0.0;
    y[1] = 1.0;
    y[2] = 1.0;
    dump();

    for (i = 1; i <= 12; ++i) {
        tout = 5.0 * i;
    retry:
        ode(f, NULL, neqn, y, &t, tout, &relerr, &abserr,
            &iflag, work, &iwork, maxnum);
        dump();
        switch (iflag) {
        case 1:
            return 1;
        case 2:
            continue;
        case 3:
            printf("tolerances too much and have been changed\n");
            goto retry;
        case 4:
            printf("too many steps...eqn. is hard\n");
            return 1;
        case 5:
            printf("eqn. appears to be stiff\n");
            return 1;
        case 6:
            printf("invalid input arguments\n");
            return 1;
        default:
            printf("unexpected iflag\n");
            return 1;
        }
    }
    return 0;
}
