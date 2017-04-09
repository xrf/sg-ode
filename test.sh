#!/bin/sh
set -eux
gcc -O3 -g -Wall -Wextra -Wconversion -pedantic ode.c ode_demo.c -lm
clang -g -std=c99 -Wall -Wextra -Wconversion -pedantic ode.c ode_demo.c -lm
if ! timeout 10 ./a.out >ode_demo.out || ! diff ode_demo.out ode_demo.txt >/dev/null; then
    git --no-pager diff --no-index ode_demo.out ode_demo.txt
    exit 1
fi
valgrind -q ./a.out >/dev/null
