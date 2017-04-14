#!/usr/bin/env python
import itertools, os, subprocess, sys, tempfile

bins = [
    {
        "name": "bin/ode_test",
        "deps": ["ode.o", "ode_demo.o", "vector.o"],
    },
    {
        "name": "bin/vector_test",
        "deps": ["mpi_vector.o_mpi", "vector.o", "vector_test.o_mpi"],
        "cc": "MPICC",
    },
]

# list of all patterns for source files
# (used for vpath builds)
vpath_patterns = ["%.c", "%.h", "%.py", "%.txt"]

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
        "targets": ["bin/ode_test.ok"],
        "macros": {
            "CC": "clang",
            "CFLAGS": "$(CFLAGS) -fsanitize=address",
        },
    },
    "clang_msan": {
        "targets": ["bin/ode_test.ok"],
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
        "targets": ["bin/ode_test_state.ok"],
        "macros": {
            "CPPFLAGS": "-DDUMP_STATE",
        },
    },
}

def bin_rule(b):
    return """
{name}: {deps}
	mkdir -p $(@D)
	$({cc}) $(CFLAGS) -o $@ {deps} $(LIBS)
"""[1:].format(name=b["name"],
               deps=" ".join(b["deps"]),
               cc=b.get("cc", "CC"))

def lint_target(target):
    return "target/{target}/.ok".format(target=target)

def lint_rule(target, lint, vpath_patterns):
    # requires GNU Make;
    # use vpath instead of VPATH to avoid interference from non-VPATH builds
    vpaths = "\n".join(
        "\t@echo 'vpath {pattern} $$(src)' >>$(@D)/Makefile"
        .format(target=target, pattern=pattern)
        for pattern in sorted(vpath_patterns))
    return """
{lint_target}:
	@mkdir -p $(@D)
	@echo 'src=../../' >$(@D)/Makefile
{vpaths}
	@echo '_all: {targets}' >>$(@D)/Makefile
	@echo 'include $$(src)Makefile' >>$(@D)/Makefile
	$(MAKE) -C $(@D) {macros}
"""[1:].format(lint_target=lint_target(target),
               vpaths=vpaths,
               targets=" ".join(sorted(lint["targets"])),
               macros=" ".join("{}='{}'".format(k, v)
                               for k, v in sorted(lint["macros"].items())))

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

subst = {
    "macros": "\n".join(sys.argv[1:]),
    "bins": "\n".join(bin_rule(b) for b in bins).rstrip("\n"),
    "deps": get_deps(list(itertools.chain(*(b["deps"] for b in bins)))),
    "lints": "\n".join(lint_rule(k, v, vpath_patterns)
                       for k, v in sorted(lints.items())),
    "lint_targets": " ".join(lint_target(k) for k in sorted(lints)),
}

with open("Makefile.in") as f:
    template = f.read()
with open("Makefile", "w") as f:
    f.write(template.format(**subst))
