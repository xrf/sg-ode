#ifndef G_WTDRPPR5HCLZDRBPLUWMIW34N5794
#define G_WTDRPPR5HCLZDRBPLUWMIW34N5794
/** @file
    Helper macros for defining vector operations. */
#ifdef __cplusplus
extern "C" {
#endif

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_1(prefix, name, var0, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], num_elems); \
        } \
    }

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_2(prefix, name, var0, var1, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, double *const restrict _v1, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            double *const var1 = _v1 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], data[1], num_elems); \
        } \
    }

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_3(prefix, name, var0, var1, var2, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, double *const restrict _v1, double *const restrict _v2, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            double *const var1 = _v1 + _i; \
            double *const var2 = _v2 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], data[1], data[2], num_elems); \
        } \
    }

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_4(prefix, name, var0, var1, var2, var3, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, double *const restrict _v1, double *const restrict _v2, double *const restrict _v3, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            double *const var1 = _v1 + _i; \
            double *const var2 = _v2 + _i; \
            double *const var3 = _v3 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], data[1], data[2], data[3], num_elems); \
        } \
    }

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_5(prefix, name, var0, var1, var2, var3, var4, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, double *const restrict _v1, double *const restrict _v2, double *const restrict _v3, double *const restrict _v4, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            double *const var1 = _v1 + _i; \
            double *const var2 = _v2 + _i; \
            double *const var3 = _v3 + _i; \
            double *const var4 = _v4 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], data[1], data[2], data[3], data[4], num_elems); \
        } \
    }

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_6(prefix, name, var0, var1, var2, var3, var4, var5, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, double *const restrict _v1, double *const restrict _v2, double *const restrict _v3, double *const restrict _v4, double *const restrict _v5, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            double *const var1 = _v1 + _i; \
            double *const var2 = _v2 + _i; \
            double *const var3 = _v3 + _i; \
            double *const var4 = _v4 + _i; \
            double *const var5 = _v5 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], data[1], data[2], data[3], data[4], data[5], num_elems); \
        } \
    }

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_7(prefix, name, var0, var1, var2, var3, var4, var5, var6, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, double *const restrict _v1, double *const restrict _v2, double *const restrict _v3, double *const restrict _v4, double *const restrict _v5, double *const restrict _v6, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            double *const var1 = _v1 + _i; \
            double *const var2 = _v2 + _i; \
            double *const var3 = _v3 + _i; \
            double *const var4 = _v4 + _i; \
            double *const var5 = _v5 + _i; \
            double *const var6 = _v6 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], data[1], data[2], data[3], data[4], data[5], data[6], num_elems); \
        } \
    }

/** Helper macro for defining map functions (no folding is done). */
#define SG_DEFINE_VECTOR_MAP_8(prefix, name, var0, var1, var2, var3, var4, var5, var6, var7, block) \
    static void name##_inner(const void *restrict ctx, const size_t _offset, double *const restrict _v0, double *const restrict _v1, double *const restrict _v2, double *const restrict _v3, double *const restrict _v4, double *const restrict _v5, double *const restrict _v6, double *const restrict _v7, const size_t _num_elems) \
    { \
        size_t _i; \
        (void)ctx; \
        for (_i = 0; _i < _num_elems; ++_i) { \
            size_t index = _offset + _i; \
            double *const var0 = _v0 + _i; \
            double *const var1 = _v1 + _i; \
            double *const var2 = _v2 + _i; \
            double *const var3 = _v3 + _i; \
            double *const var4 = _v4 + _i; \
            double *const var5 = _v5 + _i; \
            double *const var6 = _v6 + _i; \
            double *const var7 = _v7 + _i; \
            (void)index; \
            { block } \
        } \
    } \
    prefix void name(void *f_ctx, \
                     SgVectorAccum *accum, \
                     const SgVectorAccum *val, \
                     size_t offset, \
                     double **data, \
                     size_t num_elems) \
    { \
        (void)accum; \
        (void)val; \
        if (num_elems) { \
            name##_inner(f_ctx, offset, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], num_elems); \
        } \
    }

#ifdef __cplusplus
}
#endif
#endif
