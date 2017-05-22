/* Compatibility header for `inline`. */
#ifdef SG_INLINE_DEFINED
# error "must not include begin header twice in succession"
#else
# if !defined inline && __STDC_VERSION__ < 199901L && !defined __cplusplus
#  if defined __GNUC__
#   define inline __inline__
#  elif defined _MSC_VER
#   define inline __inline
#  else
#   define inline static
#  endif
#  define SG_INLINE_DEFINED 1
# else
#  define SG_INLINE_DEFINED 0
# endif
#endif
