/* Compatibility header for `restrict`. */
#ifdef SG_RESTRICT_DEFINED
# error "must not include begin header twice in succession"
#else
# if !defined restrict && __STDC_VERSION__ < 199901L
#  if defined __GNUC__
#   define restrict __restrict__
#  elif defined _MSC_VER
#   define restrict __restrict
#  else
#   define restrict
#  endif
#  define SG_RESTRICT_DEFINED 1
# else
#  define SG_RESTRICT_DEFINED 0
# endif
#endif
