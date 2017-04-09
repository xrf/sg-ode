#!/bin/sh
set -eu
gcc -O3 -g -Wall -Wextra -Wconversion -pedantic ode.c ode_demo.c -lm
clang -g -Wall -Wextra -Wconversion -pedantic ode.c ode_demo.c -lm
if ! ./a.out >ode_demo.out || ! diff ode_demo.out ode_demo.txt; then
    git diff --no-index ode_demo.out ode_demo.txt
    exit 1
fi
