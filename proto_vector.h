#ifndef G_LD1RHWBBXD615K1PEKWNYBSHDF5Q8
#define G_LD1RHWBBXD615K1PEKWNYBSHDF5Q8
/** @file
    Proto vector
*/
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/** A opaque handle to a vector. */
typedef void *Vector;

/** A pointer to an accumulator array. */
typedef void *Accum;

/** The type of an accumulator.

      - A positive integer indicates that the accumulator is an array of
        `unsigned char`, with length given by the absolute value of the type.

      - A negative integer indicates that the accumulator is an array of
        `double`, with length given by the absolute value of the type.

      - Zero indicates that the accumulator is empty and should be ignored
        entirely (the `Accum` pointer may not even be valid).

    It is important to specify the correct accumulator type, as the
    `ProtoVector` implementation may use this information to transmit the
    accumulator across different processes.
*/
typedef int AccumType;

/** An operation suitable for use in `#ProtoVectorVtable.fold_map`.

    @param[in,out] f_ctx      The same `f_ctx` pointer provided to `fold_map`.
    @param[in,out] accum      The accumulated value so far.
    @param[in]     val        A value that needs to be merged into `accum`.
    @param[in]     offset     Offset of this slice relative to the beginning.
    @param[in,out] data       Partial slices of the vectors.
    @param[in]     num_elems  Number of elements in each slice.

    The function must update `accum` to the monoidal summation of:

      - `accum`,
      - `val`, and
      - all values derived from the subvector elements (starting at index 0
        and ending at `num_elems - 1`)

    The function may also perform an element-wise operation on the given
    subvectors.  It can modify the subvectors, in which case the changes will
    be reflected on the original vectors.

    For a given value of `f_ctx`, function `f` must be thread-safe.  That is,
    it must not access global variables without synchronization.
*/
typedef void (*FoldMapFn)(void *f_ctx,
                          Accum accum,
                          Accum val,
                          size_t offset,
                          double *restrict *data,
                          size_t num_elems);

/** A proto-vector allows element-wise operations and reductions on a vector,
    as well as creation and destruction.

    All functions within this struct require exclusive access to the `self`
    object.  Therefore, one must not call any of `ProtoVector`'s functions
    with the same `self` object while any another function is running.
*/
struct ProtoVectorVtable {
    /** Creates a new uninitialized vector with a fixed length determined by
        the `ProtoVector`. */
    Vector (*create)(void *self);

    /** Destroys a vector. */
    void (*destroy)(void *self, Vector vector);

    /** Applies an element-wise operation (map) to a set of vectors and then
        accumulates a combined result (fold) of a monoidal operation.

        This is a generic interface for performing many kinds of operations on
        vectors, possibly in parallel.

        @param[in,out] self
        A pointer to the associated `ProtoVector`.

        @param[in,out] accum
        A pointer to an array containing the monoidal identity.  This value
        shall be updated with the result at the end of the calculation.  If
        `accum_type` is zero, the pointer is ignored.

        @param[in] accum_type
        If positive, the accumulator shall be an array of `unsigned char` with
        length `accum_type`.  If negative, the accumulator shall be an array
        of `double` with length `-accum_type`.  If zero, `accum` is ignored.

        @param[in] f
        The operation that is being performed.  See `#FoldMapFn`.

        @param[in] offset
        The starting index.

        @param[in,out] vectors
        The vectors that are to be read and possibly modified.  If a vector is
        being modified, it must appear only once inside the array of
        `vectors`.

        @param[in] num_vectors
        The number of `vectors`.

    */
    void (*fold_map)(void *self,
                     Accum accum,
                     AccumType accum_type,
                     FoldMapFn f,
                     void *f_ctx,
                     size_t offset,
                     Vector *vectors,
                     size_t num_vectors);
};

/** Sum the values of a vector. */
double vector_sum(void *self,
                  const struct ProtoVectorVtable *self_vtable,
                  Vector vector);

#ifdef __cplusplus
}
#endif
#endif
