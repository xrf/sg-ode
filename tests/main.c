#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ode.h"

extern const double ts[], y_initial[], relerr, abserr;

void f(double, const double *, double *);

static void wrapper(void *ctx, double t, const double *y, double *yp)
{
    (void)ctx;
    f(t, y, yp);
}

/* define our own because isfinite causes spurious warnings on Clang */
static bool is_nan(double x) {
    return x != x;
}

static void print_double_array(const double *a, size_t n)
{
    size_t i;
    printf("[");
    for (i = 0; i < n; ++i) {
        if (i != 0) {
            printf(", ");
        }
        printf("%25.17g", a[i]);
    }
    printf("]");
}

static void print_state(const int *iwork, const double *work, size_t work_len)
{
    size_t i;
    printf("{\"iwork\": [");
    for (i = 0; i < 5; ++i) {
        if (i == 0) {
            printf(", ");
        }
        printf("%5i", iwork[i]);
    }
    printf("], \"work\": ");
    print_double_array(work, work_len);
    printf("}\n");
}

int main(int argc, char **argv)
{
    bool dump_state = false;
    size_t i, neqn, work_len;
    int iwork[5] = {0};
    double t = ts[0];
    double *y, *work;

    for (i = 1; i < (size_t)argc; ++i) {
        if (strcmp(argv[i], "--dump-state") == 0) {
            dump_state = true;
        } else {
            fprintf(stderr, "unknown option: %s\n", argv[i]);
            fflush(stderr);
            return 1;
        }
    }

    for (neqn = 0; !is_nan(y_initial[neqn]); ++neqn) {}
    work_len = 100 + 21 * neqn;

    y = malloc(neqn * sizeof(*y));
    work = malloc(work_len * sizeof(*work));

    if (dump_state) {
        /* mark uninitialized values as NAN; unfortunately this breaks the
           memory sanitizer so we only do this if necessary */
        for (i = 0; i < work_len; ++i) {
            work[i] = NAN;
        }
    }

    for (i = 0; i < neqn; ++i) {
        y[i] = y_initial[i];
    }

    for (i = 0; !is_nan(ts[i]); ++i) {
        int status = 0;
        if (i > 0) {
            status = sg_ode(NULL, &wrapper, neqn, y, &t, ts[i],
                            relerr, abserr, 0, work, iwork);
        }

        printf("{\"t\": %25.17g, \"y\": ", t);
        print_double_array(y, neqn);
        printf(", \"status\": %1i", status);
        if (dump_state) {
            printf(", \"state\": ");
            print_state(iwork, work, work_len);
        }
        printf("}\n");

        if (status != 0) {
            break;
        }
    }

    free(work);
    free(y);
    return 0;
}
