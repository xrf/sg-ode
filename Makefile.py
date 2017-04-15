#!/usr/bin/env python
import glob, itertools, os, subprocess, sys, tempfile

subst = {}

bins = {
    "bin/vector_test": {
        "deps": ["mpi_vector.o_mpi", "vector.o", "vector_test.o_mpi"],
        "cc": "MPICC",
    },
}

tests = []
test_targets = []

for fn in sorted(glob.glob("tests/*.c")):
    if fn in ["tests/main.c"]:
        continue
    name, _ = os.path.splitext(os.path.basename(fn))
    bins["bin/{name}_test".format(**locals())] = {
        "deps": ["tests/main.o", "ode.o", "vector.o",
                 os.path.splitext(fn)[0] + ".o"],
    }
    target = "bin/{name}_test.ok".format(**locals())
    test_targets.append(target)
    tests.append("""
{target}: bin/{name}_test tests/{name}$(TESTSUFFIX).txt_out
	@mkdir -p $(@D)
	$(harness) bin/{name}_test $(TESTFLAGS) >tests/{name}.out
	@git --no-pager diff --exit-code --no-index $(GITDIFFFLAGS) tests/{name}.out tests/{name}$(TESTSUFFIX).txt_out
	@touch $@
"""[1:].format(**locals()))

subst["tests"] = "\n".join(tests)
subst["test_targets"] = " ".join(test_targets)

def bin_rule(name, b):
    deps = " ".join(b["deps"])
    cc = b.get("cc", "CC")
    return """
{name}: {deps}
	mkdir -p $(@D)
	$({cc}) $(CFLAGS) -o $@ {deps} $(LIBS)
"""[1:].format(**locals())

subst["bins"] = "\n".join(bin_rule(n, b)
                          for n, b in sorted(bins.items())).rstrip("\n")

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
lints = {
    "gcc_valgrind": {
        "targets": ["check"],
        "macros": {
            "CC": "gcc",
            "CFLAGS": "$(CFLAGS) -O3",
            "HARNESS": "valgrind --error-exitcode=1 -q",
            "MPICC": "OMPI_CC=$$(CC) MPICH_CC=$$(CC) $(MPICC)",
        },
    },
    "clang_asan": {
        "targets": ["check_ode"],
        "macros": {
            "CC": "clang",
            "CFLAGS": "$(CFLAGS) -fsanitize=address",
        },
    },
    "clang_msan": {
        "targets": ["check_ode"],
        "macros": {
            "CC": "clang",
            "CFLAGS": "$(CFLAGS) -fsanitize=memory",
        },
    },
    "clang_ubsan": {
        "targets": ["check"],
        "macros": {
            "CC": "clang",
            "CFLAGS": "$(CFLAGS) -fsanitize=undefined",
            "MPICC": "OMPI_CC=$$(CC) MPICH_CC=$$(CC) $(MPICC)",
        },
    },
    "dump_state": {
        "targets": ["check_ode"],
        "macros": {
            "TESTSUFFIX": "_state",
            "TESTFLAGS": "--dump-state",
        },
    },
}

def lint_target(target):
    return "target/{target}/.ok".format(target=target)

def lint_rule(target, lint, vpath_patterns):
    # requires GNU Make;
    # use vpath instead of VPATH to avoid interference from non-VPATH builds
    vpaths = "\n".join(
        "\t@echo 'vpath {pattern} $$(src)' >>$(@D)/Makefile"
        .format(**locals())
        for pattern in sorted(vpath_patterns))
    tar = lint_target(target)
    targets = " ".join(sorted(lint["targets"]))
    macros = " ".join("{}='{}'".format(k, v)
                      for k, v in sorted(lint["macros"].items()))
    return """
{tar}:
	@mkdir -p $(@D)
	@echo 'src=../../' >$(@D)/Makefile
{vpaths}
	@echo '_all: {targets}' >>$(@D)/Makefile
	@echo 'include $$(src)Makefile' >>$(@D)/Makefile
	$(MAKE) -C $(@D) {macros}
"""[1:].format(**locals())

subst["lints"] = "\n".join(lint_rule(k, v, vpath_patterns)
                           for k, v in sorted(lints.items()))
subst["lint_targets"] = " ".join(lint_target(k) for k in sorted(lints))

def get_deps(objs):
    deps = set()
    for obj in objs:
        src = os.path.splitext(obj)[0] + ".c"
        with tempfile.NamedTemporaryFile(mode="r") as f:
            e = subprocess.call(["cc", "-o", "/dev/null", "-M", "-MG", "-MM",
                                 "-MT", obj, "-MF", f.name, src])
            if e:
                sys.stderr.write("warning: dependency generation for "
                                 "{obj} failed\n")
                sys.stderr.flush()
            deps.update(s.strip() for s in f.read().split("\n\n"))
    return "\n\n".join(sorted(deps))

subst["deps"] = get_deps(list(itertools.chain(*(b["deps"]
                                                for b in bins.values()))))

subst["macros"] = "\n".join(sys.argv[1:])

with open("Makefile.in") as f:
    template = f.read()
with open("Makefile", "w") as f:
    f.write(template.format(**subst))
