#!/bin/sh
set -eux
CFLAGS="-std=c99 -Wall -Wextra -Wconversion -pedantic"

gcc -O3 -Wall -Wextra -Wconversion -pedantic ode.c ode_demo.c -lm

clang -fsanitize=undefined $CFLAGS ode.c ode_demo.c -lm
if ! ./a.out >ode_demo.out || ! diff ode_demo.out ode_demo.txt >/dev/null; then
    if [ "`wc -c ode_demo.out | cut -f 1 -d ' '`" -ne 0 ]; then
        git --no-pager diff --no-index ode_demo.out ode_demo.txt
    fi
    exit 1
fi

valgrind --error-exitcode=1 -q ./a.out >/dev/null
clang -fsanitize=memory $CFLAGS ode.c ode_demo.c -lm && ./a.out >/dev/null
clang -fsanitize=address $CFLAGS ode.c ode_demo.c -lm && ./a.out >/dev/null
