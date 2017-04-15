/* Solves the defining equations of Jacobian elliptic functions

       y₀′ =     y₁ y₂    y₀(0) = 0
       y₁′ =    -y₀ y₂    y₁(0) = 1
       y₂′ = -k² y₀ y₁    y₂(0) = 1

   for t[i] ∈ {0, 5, 10, 15, …, 60}.  The analytic solution is

       y₀ = sn(t / k²)
       y₁ = cn(t / k²)
       y₂ = dn(t / k²)

   in this case we use k² = 0.51.
*/
#include <math.h>

const double y_initial[] = {0.0, 1.0, 1.0, NAN};

void f(double t, const double *y, double *yp)
{
    static const double ksq = 0.51;
    (void)t;
    yp[0] = y[1] * y[2];
    yp[1] = -y[0] * y[2];
    yp[2] = -ksq * y[0] * y[1];
}
