#!/usr/bin/env python3
import os

# we define a separate inner function because Clang and GCC seem to ignore
# 'restrict' on the array pointer parameters and non-parameter variables

def write_macro(f, n):
    f.write("""
/** Helper macro for defining map functions (no folding is done). */
"""[1:])
    f.write(f"#define DEF_MAP_FN_{n}(prefix, name, ")
    for i in range(n):
        f.write(f"var{i}, ")
    f.write("block) \\\n")
    f.write("    static void name##_inner(const size_t _offset, ")
    for i in range(n):
        f.write(f"double *const restrict _v{i}, ")
    f.write("const size_t _num_elems) \\\n")
    f.write("""
    { \\
        size_t _i; \\
        for (_i = 0; _i < _num_elems; ++_i) { \\
            size_t index = _offset + _i; \\
"""[1:])
    for i in range(n):
        f.write(f"            double *const var{i} = _v{i} + _i; \\\n")
    f.write("""
            (void)index; \\
            { block } \\
        } \\
    } \\
    prefix void name(void *f_ctx, \\
                     Accum accum, \\
                     Accum val, \\
                     size_t offset, \\
                     double *restrict *data, \\
                     size_t num_elems) \\
    { \\
        (void)f_ctx; \\
        (void)accum; \\
        (void)val; \\
        if (num_elems) { \\
            name##_inner(offset, """[1:])
    for i in range(n):
        f.write(f"data[{i}], ")
    f.write("""num_elems); \\
        } \\
    }

""")

with open(os.path.splitext(__file__)[0] + ".h", "w") as f:
    f.write("""
#ifndef G_WTDRPPR5HCLZDRBPLUWMIW34N5794
#define G_WTDRPPR5HCLZDRBPLUWMIW34N5794
/** @file
    Helper macros */
#ifdef __cplusplus
extern "C" {
#endif

"""[1:])
    for n in range(8):
        write_macro(f, n)
    f.write("""
#ifdef __cplusplus
}
#endif
#endif
"""[1:])
