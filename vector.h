#ifndef G_LD1RHWBBXD615K1PEKWNYBSHDF5Q8
#define G_LD1RHWBBXD615K1PEKWNYBSHDF5Q8
/** @file
    A generic interface for elementwise and aggregate operations on vectors.
*/
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/** A vector is an opaque data type containing a finite number of `double`s.  */
typedef void SgVector;

/** An accumulator is an array of some kind. */
typedef void SgVectorAccum;

/** The type of an accumulator.

      - A positive integer indicates that the accumulator is an array of
        `unsigned char`, with length given by the absolute value of the type.

      - A negative integer indicates that the accumulator is an array of
        `double`, with length given by the absolute value of the type.

      - Zero indicates that the accumulator is empty and should be ignored
        entirely (the `Accum` pointer may not even be valid).

    It is important to specify the correct accumulator type, as the
    `SgVectorDriver` implementation may use this information to transmit the
    accumulator across different processes.
*/
typedef int SgVectorAccumType;

/** An operation suitable for use in `#SgVectorDriverVt.operate`.

    @param[in,out] f_ctx      The same `f_ctx` pointer provided to `operate`.
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
typedef void SgVectorOperation(void *f_ctx,
                               SgVectorAccum *accum,
                               const SgVectorAccum *val,
                               size_t offset,
                               double **data,
                               size_t num_elems);

/** Minimal structure of an `SgVectorDriverBase` object (“base class”).  An
    actual instance may hold additional information, but all drivers must
    store `SgVectorDriverBase` as their first element. */
struct SgVectorDriverBase {
    /** The length of all vectors created by this driver.  Once a driver has
        been created, this attribute must not be modified. */
    size_t len;
};

/** The minimal set of operations that a vector driver must support.

    A vector driver allows element-wise operations and reductions on a vector,
    as well as creation and destruction.

    All functions within this struct require exclusive access to the `self`
    object.  Therefore, one must not call any of `SgVectorDriver`'s functions
    with the same `self` object while any another function is running.
*/
struct SgVectorDriverVt {
    /** Creates a new uninitialized vector with a fixed length determined by
        the `SgVectorDriver`. */
    SgVector *(*try_new)(struct SgVectorDriverBase *self);

    /** Destroys a vector. */
    void (*del)(struct SgVectorDriverBase *self, SgVector *vector);

    /** Applies an element-wise operation (map) to a set of vectors and then
        accumulates a combined result (fold) of a monoidal operation.

        This is a generic interface for performing many kinds of operations on
        vectors, possibly in parallel.

        @param[in,out] self
        A pointer to the associated `SgVectorDriver`.

        @param[in,out] accum
        A pointer to an array containing the monoidal identity.  This value
        shall be updated with the result at the end of the calculation.  If
        `accum_type` is zero, the pointer is ignored.

        @param[in] accum_type
        If positive, the accumulator shall be an array of `unsigned char` with
        length `accum_type`.  If negative, the accumulator shall be an array
        of `double` with length `-accum_type`.  If zero, `accum` is ignored.

        @param[in] f
        The operation that is being performed.  See `#SgVectorOperation`.

        @param[in] offset
        The starting index.

        @param[in,out] vectors
        The vectors that are to be read and possibly modified.  If a vector is
        being modified, it must appear only once inside the array of
        `vectors`.

        @param[in] num_vectors
        The number of `vectors`.

    */
    void (*operate)(struct SgVectorDriverBase *self,
                    SgVectorAccum *accum,
                    SgVectorAccumType accum_type,
                    SgVectorOperation *f,
                    void *f_ctx,
                    size_t offset,
                    SgVector **vectors,
                    size_t num_vectors);
};

/** A vector driver.  This is just a driver pointer combined with its vtable
    pointer for convenience. */
struct SgVectorDriver {
    struct SgVectorDriverBase *data;
    struct SgVectorDriverVt *vtable;
};

/** Creates a vectorm, returning `NULL` if it fails. */
SgVector *sg_vector_try_new(struct SgVectorDriver drv);

/** Creates a vector. */
SgVector *sg_vector_new(struct SgVectorDriver drv);

/** Destroys a vector. */
void sg_vector_del(struct SgVectorDriver drv, SgVector *vector);

/** Gets the length of any vector created by this vector driver. */
size_t sg_vector_len(struct SgVectorDriver drv);

/** Perform a generic operation on some vectors.  See
    `#SgVectorDriverVt.operate`. */
void sg_vector_operate(struct SgVectorDriver drv,
                       SgVectorAccum *accum,
                       SgVectorAccumType accum_type,
                       SgVectorOperation *f,
                       void *f_ctx,
                       size_t offset,
                       SgVector **vectors,
                       size_t num_vectors);

/** Sums the values of a vector. */
double sg_vector_sum(struct SgVectorDriver drv, const SgVector *vector);

/** The underlying operation used by `sg_vector_sum`. */
extern SgVectorOperation vector_sum_operation;

struct SgBasicVectorDriver {
    struct SgVectorDriverBase base;
};

struct SgBasicVectorDriver sg_basic_vector_driver_new(size_t len);

struct SgVectorDriver sg_basic_vector_driver_get(struct SgBasicVectorDriver *);

extern struct SgVectorDriverVt SG_BASIC_VECTOR_DRIVER_VT;

#ifdef __cplusplus
}
#endif
#endif
