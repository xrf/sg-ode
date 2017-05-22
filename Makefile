CFLAGS?=-g -O2
LDLIBS?=-lm
MPICC?=mpicc
MPIEXEC?=mpiexec
DIFF=git --no-pager diff --exit-code --no-index
harness=$(HARNESS) timeout 15

CDEPFLAGS=-MD -MP -MT $@ -MF $*.dep
MPICDEPFLAGS=-MD -MP -MT $@ -MF $*.dep_mpi

-include config.mk

all:

clean:
	rm -fr bin target *.o *.o_mpi */*.o */*.o_mpi */*.out

check: check_ode bin/vector_test.ok

check_ode: bin/jacobian_elliptic_a_test.ok bin/jacobian_elliptic_b_test.ok

doc:
	doxygen

test: check

bin/%_test.ok: bin/%_test tests/%$(TESTSUFFIX).txt
	mkdir -p $(@D)
	$(harness) $(@:.ok=) $(TESTFLAGS) >$(@:.ok=.out)
	$(DIFF) $(@:.ok=.out) $(word 2,$^)
	@touch $@

bin/vector_test.ok: bin/vector_test
	mkdir -p $(@D)
	$(harness) $(MPIEXEC) -np 4 bin/vector_test
	@touch $@

bin/jacobian_elliptic_a_test: tests/main.o ode.o vector.o tests/jacobian_elliptic_a.o
	mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/jacobian_elliptic_b_test: tests/main.o ode.o vector.o tests/jacobian_elliptic_b.o
	mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/vector_test: vector_test.o_mpi mpi_vector.o_mpi vector.o
	mkdir -p $(@D)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDLIBS)

vector_macros.h: vector_macros.py
	mkdir -p $(@D)
	./$<

# ------------------------------------------------------------------------
# Lint rules
#
# Notes
#
#   - All sanitizers are run in subdirectories of 'target' to avoid polluting
#     the main tree.
#   - If the sanitizer fails to mmap, it's probably due to ulimit.
#   - Address and memory sanitizers require all external libraries to be
#     recompiled, so we're avoid sanitizing any tests that require MPI.

target/clang_asan: makeflags=CC=clang CFLAGS='$(CFLAGS) -fsanitize=address' check_ode
target/clang_msan: makeflags=CC=clang CFLAGS='$(CFLAGS) -fsanitize=memory' check_ode
target/clang_ubsan: makeflags=CC=clang CFLAGS='$(CFLAGS) -fsanitize=undefined' MPICC='OMPI_CC=$$(CC) MPICH_CC=$$(CC) $(MPICC)' check
target/dump_state: makeflags=TESTFLAGS=--dump-state TESTSUFFIX=_state check_ode
target/gcc_valgrind: makeflags=CC=gcc CFLAGS='$(CFLAGS) -O3' HARNESS='valgrind --error-exitcode=1 -q' MPICC='OMPI_CC=$$(CC) MPICH_CC=$$(CC) $(MPICC)' check

lints!=sed -n 's|^\(target/[^:%][^:%]*\).*|\1|p' Makefile

lint: $(lints)
	echo $^

# We must add a list of all patterns for source files for vpath builds
# If make complains about missing files when linting this is probably why
$(lints):
	@mkdir -p $@
	@printf '' >$@/Makefile
	@echo vpath %.c ../../ >>$@/Makefile
	@echo vpath %.h ../../ >>$@/Makefile
	@echo vpath %.inl ../../ >>$@/Makefile
	@echo vpath %.py ../../ >>$@/Makefile
	@echo vpath %.txt ../../ >>$@/Makefile
	@echo include ../../Makefile >>$@/Makefile
	+$(MAKE) -C $@ $(makeflags)

# ------------------------------------------------------------------------
# Generic rules

.c.o:
	mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CDEPFLAGS) -c -o $@ $<

.c.o_mpi:
	mkdir -p $(@D)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(MPICDEPFLAGS) -c -o $@ $<

.SUFFIXES: .c .o .o_mpi

.PHONY: Makefile all check check_ode doc lint $(lints) test

-include $(wildcard *.dep) $(wildcard *.dep_mpi) $(wildcard */*.dep) $(wildcard */*.dep_mpi)
