CFLAGS?=-g -O2 -fvisibility=hidden
LDLIBS?=-lm
PREFIX?=/usr/local

CDEPFLAGS=-MD -MP -MT $@ -MF $*.dep
RPATH_ORIGIN=-Wl,-rpath,'$$ORIGIN'

# These control how the shared library is built:
#
#   - SHAREDNAME is the filename of the shared library
#   - SHAREDFLAGS are the flags needed to link the shared library
#   - SHAREDLN are commands used to create the symbolic links
#
SHAREDCFLAGS=-fPIC
SHAREDNAME=lib$(1).so.$(2).$(3).$(4)
SHAREDFLAGS=-shared -Wl,-soname,lib$(1).so.$(2)
SHAREDLN=ln -fs lib$(1).so.$(2).$(3).$(4) lib$(1).so.$(2) && ln -fs lib$(1).so.$(2).$(3).$(4) lib$(1).so

DIFF=git --no-pager diff --exit-code --no-index
harness=$(HARNESS) timeout 15

-include config.mk

CFLAGS+=$(SHAREDCFLAGS)
top?=.
sg_ode/%.o: CPPFLAGS+=-DSG_BUILD
tests/%.o: CPPFLAGS+=-I$(top)

major=1
minor=2
patch=0
libsgode=$(call SHAREDNAME,sgode,$(major),$(minor),$(patch))

all: target/build/$(libsgode)

clean:
	rm -fr bin lib target *.dep *.o */*.dep */*.o */*.out

install: all
	install -d $(DESTDIR)$(PREFIX)/include/sg_ode $(DESTDIR)$(PREFIX)/lib
	install -m644 -t $(DESTDIR)$(PREFIX)/include sg_ode.h
	install -m644 -t $(DESTDIR)$(PREFIX)/include/sg_ode sg_ode/extern.h sg_ode/inline_begin.h sg_ode/inline_end.h sg_ode/ode.h sg_ode/restrict_begin.h sg_ode/restrict_end.h sg_ode/vector.h
	install -m755 -t $(DESTDIR)$(PREFIX)/lib target/build/$(libsgode)
	cd $(DESTDIR)$(PREFIX)/lib $(LN_SHARED) && $(call SHAREDLN,sgode,$(major),$(minor),$(patch))

target/build/$(libsgode): sg_ode/ode.o sg_ode/vector.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) $(call SHAREDFLAGS,sgode,$(major),$(minor),$(patch)) -o $@ $^ $(LDLIBS)
	cd $(@D) && $(call SHAREDLN,sgode,$(major),$(minor),$(patch))

check: target/build/arguments_test.ok target/build/jacobian_elliptic_a_test.ok target/build/jacobian_elliptic_b_test.ok target/build/stiff_test.ok

target/build/%_test.ok: target/build/%_test tests/%$(TESTSUFFIX).txt
	@mkdir -p $(@D)
	$(harness) $< $(TESTFLAGS) >$(@:.ok=.out)
	$(DIFF) $(word 2,$^) $(@:.ok=.out)
	touch $@

target/build/%_test: tests/main.o tests/%.o target/build/$(libsgode)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) $(RPATH_ORIGIN) -Ltarget/build -o $@ $(wordlist 1,2,$^) $(LDLIBS) -lsgode

target/build/arguments_test.ok: TESTFLAGS+=2>$@.log

target/build/arguments_test: tests/arguments.o target/build/$(libsgode)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) $(RPATH_ORIGIN) -Ltarget/build -o $@ $< $(LDLIBS) -lsgode

tests/%$(TESTSUFFIX).txt:
	touch $@

sg_ode/vector_macros.h: sg_ode/vector_macros.py
	@mkdir -p $(@D)
	./$<

doc:
	doxygen

deploy-doc: doc/html/.git/config doc
	cd $(<D)/.. && \
	git add -A && \
	git commit --amend -q -m Autogenerated && \
	git push -f origin master:gh-pages

doc/html/.git/config:
	mkdir -p $(@D)
	url=`git remote -v | grep origin | awk '{ printf "%s", $$2; exit }'` && \
	cd $(@D)/.. && \
	git init && \
	git config user.name Bot && \
	git config user.email "<>" && \
	git commit -m _ --allow-empty && \
	git remote add origin "$$url"

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

target/clang_asan: makeflags=CC=clang CFLAGS='$(CFLAGS) -fsanitize=address' check
target/clang_msan: makeflags=CC=clang CFLAGS='$(CFLAGS) -fsanitize=memory' check
target/clang_ubsan: makeflags=CC=clang CFLAGS='$(CFLAGS) -fsanitize=undefined' check
target/dump_state: makeflags=TESTFLAGS=--dump-state TESTSUFFIX=_state check
target/gcc_valgrind: makeflags=CC=gcc CFLAGS='$(CFLAGS) -O3' HARNESS='valgrind --error-exitcode=1 -q' check

lints!=sed -n 's|^\(target/[^:%][^:%]*\): makeflags.*|\1|p' Makefile

lint: $(lints)
	echo $^

# We must add a list of all patterns for source files for vpath builds
# If make complains about missing files when linting this is probably why
$(lints):
	@mkdir -p $@
	@echo top=../../ >$@/Makefile
	@echo 'vpath %.c $$(top)' >>$@/Makefile
	@echo 'vpath %.h $$(top)' >>$@/Makefile
	@echo 'vpath %.inl $$(top)' >>$@/Makefile
	@echo 'vpath %.py $$(top)' >>$@/Makefile
	@echo 'vpath %.txt $$(top)' >>$@/Makefile
	@echo 'include $$(top)Makefile' >>$@/Makefile
	+$(MAKE) -C $@ $(makeflags)

# ------------------------------------------------------------------------
# Generic rules

.c.o:
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CDEPFLAGS) -c -o $@ $<

.SUFFIXES: .c .o

.PHONY: Makefile all check deploy-doc doc install lint $(lints)

.SECONDARY:

.DELETE_ON_ERROR:

-include $(wildcard *.dep) $(wildcard */*.dep)
