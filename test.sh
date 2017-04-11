#!/bin/sh
set -eux
CFLAGS="-std=c99 -Wall -Wextra -Wconversion -pedantic"

objs='ode.c ode_demo.c proto_vector.c'

OMPI_CC=gcc mpicc -O3 -Wall -Wextra -Wconversion -pedantic $objs -lm

export OMPI_CC=clang

mpicc -fsanitize=undefined $CFLAGS $objs -lm
if ! ./a.out >ode_demo.out || ! diff ode_demo.out ode_demo.txt >/dev/null; then
    if [ "`wc -c ode_demo.out | cut -f 1 -d ' '`" -ne 0 ]; then
        git --no-pager diff --no-index ode_demo.out ode_demo.txt
    fi
    exit 1
fi

valgrind --error-exitcode=1 -q ./a.out >/dev/null
mpicc -fsanitize=memory $CFLAGS $objs -lm && ./a.out >/dev/null
mpicc -fsanitize=address $CFLAGS $objs -lm && ./a.out >/dev/null

mpicc $CFLAGS proto_vector.c proto_vector_test.c && mpiexec -np 4 ./a.out
