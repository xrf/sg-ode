#!/bin/sh
set -eux
CFLAGS="-std=c99 -Wall -Wextra -Wconversion -pedantic"

objs='ode.c ode_demo.c vector.c'

gcc -O3 -Wall -Wextra -Wconversion -pedantic $objs -lm

clang -fsanitize=undefined $CFLAGS $objs -lm
if ! ./a.out >ode_demo.out || ! diff ode_demo.out ode_demo.txt >/dev/null; then
    if [ "`wc -c ode_demo.out | cut -f 1 -d ' '`" -ne 0 ]; then
        git --no-pager diff --no-index ode_demo.out ode_demo.txt
    fi
    exit 1
fi

valgrind --error-exitcode=1 -q ./a.out >/dev/null
clang -fsanitize=memory $CFLAGS $objs -lm && ./a.out >/dev/null
clang -fsanitize=address $CFLAGS $objs -lm && ./a.out >/dev/null

OMPI_CC=clang mpicc $CFLAGS vector.c mpi_vector.c vector_test.c && mpiexec -np 4 ./a.out
