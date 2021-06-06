# Parallelizable Shampine-Gordon ODE solver [![Build status](https://github.com/xrf/sg-ode/actions/workflows/build.yml/badge.svg)](https://github.com/xrf/sg-ode/actions/workflows/build.yml)

**Quick links:** [documentation](https://xrf.github.io/sg-ode), [releases][rel].

Solves a general system of [ordinary differential equations (ODEs)](https://en.wikipedia.org/wiki/Ordinary_differential_equation) using an [algorithm by L. F. Shampine and M. K. Gordon][sg]:

    y′(t) = f(t, y(t))

The code was originally translated from the [Netlib version](http://www.netlib.org/ode/ode.f) with the aid of [f2c](http://netlib.org/f2c/).

## Installation

[Download][rel] the latest tarball and run:

    make PREFIX=/usr/local install

Replace `/usr/local` (the default `PREFIX`) with wherever you want it to be installed.  Afterward, the library will be installed to `$PREFIX/lib/libsgode.so` and the header files will be installed to `$PREFIX/include/sg_ode/`.

Arch Linux users can use the `PKGBUILD` script instead.

## Usage

To include the headers in your source code, write:

~~~c
#include <sg_ode/ode.h>
~~~

Be sure that `$PREFIX/include` is part of your header search path (`-I $PREFIX/include`).

To link with the library, pass `-l sgode` and make sure `$PREFIX/lib` is part of your library search path (`-L $PREFIX/include`).

## Examples

The `tests/main.c` directory provides a generic test driver for integrating ODE equations using this library.  It must be linked with either `tests/harmonic_oscillator.c`, `tests/jacobian_elliptic_a.c`, or `tests/jacobian_elliptic_b.c` to create a full program.

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

[rel]: https://github.com/xrf/sg-ode/releases
[sg]:  http://books.google.com/books?id=3Yl2nQEACAAJ
