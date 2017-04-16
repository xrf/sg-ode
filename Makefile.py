#!/usr/bin/env python
import glob, itertools, os, re, sys
import makegen2 as mk

def make_src(path, suffix="_out"):
    '''Source files for non-inference rules should be wrapped using this rule
    to ensure VPATH works properly.'''
    _, ext = os.path.splitext(path)
    if not re.match("\.\w+\Z", ext):
        raise ValueError("weird path extension: {!r}".format(ext))
    if not re.match("\w+\Z", suffix):
        raise ValueError("invalid suffix: {!r}".format(ext))
    return mk.Make(path + "_out",
                   [mk.InferenceRule(ext, ext + suffix, ["ln $< $@"]).make()])

# ------------------------------------------------------------------------
# Main rules

makefile = mk.make()

default = makefile.append(mk.make("all", phony=True))
makefile |= mk.make("doc", [], ["doxygen"], phony=True)
makefile |= mk.make("Makefile", [],
                    ["@[ ! -x ./config.status ] || ./config.status"],
                    phony=True)
makefile |= mk.make("vector_macros.h", [make_src("vector_macros.py")],
                    lambda script: ["./{}".format(script)],
                    clean=False)

# ------------------------------------------------------------------------
# Test rules

check_ode = makefile.append(mk.make("check_ode", phony=True))
check = makefile.append(mk.make("check", [check_ode], phony=True))
makefile |= mk.make("test", [check], phony=True) # alias

t = mk.make_bin("bin/vector_test",
                [
                    mk.make_c_obj("vector_test.c", via=mk.C_MPI_INFERENCE_RULE),
                    mk.make_c_obj("mpi_vector.c", via=mk.C_MPI_INFERENCE_RULE),
                    "vector.c",
                ],
                ld="$(MPICC) $(CFLAGS)")
check |= mk.make(t.name + ".ok",
                 [t],
                 ["$(harness) $(MPIEXEC) -np 4 {}".format(t.name),
                  "@touch $@"])

for path in sorted(glob.glob("tests/*.c")):
    if path in ["tests/main.c"]:
        continue
    name, _ = os.path.splitext(os.path.basename(path))
    exe_name = "bin/{}_test".format(name)
    check_ode |= mk.make(
        exe_name + ".ok",
        [
            mk.make_bin(exe_name, ["tests/main.c", "ode.c", "vector.c", path]),
            make_src("tests/{}$(TESTSUFFIX).txt".format(name)),
        ],
        lambda _, txt: [
            "$(harness) $(@:.ok=) $(TESTFLAGS) >$(@:.ok=.out)",
            "$(DIFF) $(@:.ok=.out) {}".format(txt),
            "@touch $@",
        ],
        macros={
            "DIFF": "git --no-pager diff --exit-code --no-index",
            "harness": "$(HARNESS) timeout 15",
        })

# ------------------------------------------------------------------------
# Lint rules

# list of all patterns for source files
# (used for vpath builds)
# if make complains about missing files when linting this is probably why
vpath_patterns = ["%.c", "%.h", "%.inl", "%.py", "%.txt"]

# Notes
#
#   - All sanitizers are run in subdirectories of 'target' to avoid polluting
#     the main tree.
#   - If the sanitizer fails to mmap, it's probably due to ulimit.
#   - Address and memory sanitizers require all external libraries to be
#     recompiled, so we're avoid sanitizing any tests that require MPI.
#

lint = makefile.append(mk.make("lint", phony=True))

def make_lint(name, targets, macros):
    return mk.make_vpath(os.path.join("target", name),
                         targets=targets,
                         macros=macros,
                         vpath_patterns=vpath_patterns)

lint |= make_lint("gcc_valgrind",
                  targets=["check"],
                  macros={
                      "CC": "gcc",
                      "CFLAGS": "$(CFLAGS) -O3",
                      "HARNESS": "valgrind --error-exitcode=1 -q",
                      "MPICC": "OMPI_CC=$$(CC) MPICH_CC=$$(CC) $(MPICC)",
                  })
lint |= make_lint("clang_asan",
                  targets=["check_ode"],
                  macros={
                      "CC": "clang",
                      "CFLAGS": "$(CFLAGS) -fsanitize=address",
                  })
lint |= make_lint("clang_msan",
                  targets=["check_ode"],
                  macros={
                      "CC": "clang",
                      "CFLAGS": "$(CFLAGS) -fsanitize=memory",
                  })
lint |= make_lint("clang_ubsan",
                  targets=["check"],
                  macros={
                      "CC": "clang",
                      "CFLAGS": "$(CFLAGS) -fsanitize=undefined",
                      "MPICC": "OMPI_CC=$$(CC) MPICH_CC=$$(CC) $(MPICC)",
                  })
lint |= make_lint("dump_state",
                  targets=["check_ode"],
                  macros={
                      "TESTSUFFIX": "_state",
                      "TESTFLAGS": "--dump-state",
                  })

# ------------------------------------------------------------------------
# Finish

makefile |= mk.make(macros=(s.split("=", 1) for s in sys.argv[1:]))
makefile.save("Makefile.in", out_name="Makefile", cache=".Makefile.cache~")
