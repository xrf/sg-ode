# Parallelizable Shampine-Gordon ODE solver [![Build status][ci]][ca]

**Quick links:** [documentation][doc], [releases][rel].

Solves a general system of [ordinary differential equations (ODEs)][ode] using
an [algorithm by L. F. Shampine and M. K. Gordon][sg]:

    y′(t) = f(t, y(t))

The code was originally translated from the [Netlib version][nl] with the aid of f2c.

[ode]: https://en.wikipedia.org/wiki/Ordinary_differential_equation
[nl]:  http://www.netlib.org/ode/ode.f
[sg]:  http://books.google.com/books?id=3Yl2nQEACAAJ

## Installation

[Download][rel] the tarball and run:

    make PREFIX=/usr/local install

Replace `/usr/local` with wherever you want it to be installed.  With the default settings, the library is installed to `/usr/local/lib/libsgode.so` and the main header file is installed to `/usr/local/include/sg_ode.h`

If you are using Arch Linux, you can use the `PKGBUILD` script instead.

## Parallelization support

To support parallelization, the algorithm utilizes an *abstract vector interface* that is sufficiently generic to support single-threaded, multi-threaded, and distributed vector representations without significantly compromising performance.  The actual implementations are provided by a *vector driver*.

To keep the library lightweight and dependency-free, the library comes with only a driver for the trivial vector representation.  It should be straightforward to implement a driver for a distributed vector representation for say MPI.

The core of the abstract vector interface lies in a single higher-order function that we have uncreatively named `operate`.  Here is a sketch of its type signature with the boring parts omitted for clarity:

~~~c
void operate(Accum total,
             void f(Accum subtotal,
                    const Accum input,
                    double data[m][n],
                    size_t n),
             Vector vectors[m],
             size_t m);
~~~

The function `operate` performs a map-like (i.e. element-wise) operation over `m` vectors followed by a fold (accumulate/reduce) operation over all the elements.  The actual operation is specified by `f`, which performs the desired operation over a set of `m` subvectors of length `n` and also accumulates the value of `input` as well as the value of the subvector elements into the `subtotal`.

## References

  * L. K. Shampine and M. K. Gordon,
    [*Computer Solution of Ordinary Differential Equations:
      The Initial Value Problem*][sg]
    (Freeman, 1975), ISBN: 0716704617.

[ca]:  https://travis-ci.org/xrf/sg-ode
[ci]:  https://travis-ci.org/xrf/sg-ode.svg?branch=master
[doc]: https://xrf.github.io/sg-ode/sg__ode_8h.html
[rel]: https://github.com/xrf/sg-ode/releases
