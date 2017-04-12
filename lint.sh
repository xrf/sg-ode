#!/bin/sh
set -eu
jobs=4
# Notes
#
#   - If the sanitizer fails to mmap, it's probably because of ulimit.
#   - Address and memory sanitizers require MPI (or any external library) to
#     be recompiled, so we're avoid sanitizing any tests that require MPI.
#   - If the linting fails, you should probably make clean before doing
#     anything else.
#
cflags="-g -O2 -std=c99 -Wall -Wconversion -Wextra -pedantic"
make -j "${jobs}" clean
make -j "${jobs}" CC=gcc CFLAGS="$cflags -O3" HARNESS="valgrind --error-exitcode=1 -q" check
make -j "${jobs}" clean
make -j "${jobs}" CC=clang CFLAGS="$cflags -fsanitize=address" bin/ode_test.ok
make -j "${jobs}" clean
make -j "${jobs}" CC=clang CFLAGS="$cflags -fsanitize=undefined" check
make -j "${jobs}" clean
make -j "${jobs}" CC=clang CFLAGS="$cflags -fsanitize=memory" bin/ode_test.ok
make -j "${jobs}" clean
