/* Compatibility header for `thread_local`. */
#ifdef SG_THREAD_LOCAL_DEFINED
# error "must not include begin header twice in succession"
#else
# if !defined thread_local && __cplusplus < 201103L
#  if __STDC_VERSION__ >= 201112L
#   define thread_local _Thread_local
#  elif defined __GNUC__
#   define thread_local __thread
#  elif defined _MSC_VER
#   define thread_local __declspec(thread)
#  else
#   error "missing support for thread-local variables"
#  endif
#  define SG_THREAD_LOCAL_DEFINED 1
# else
#  define SG_THREAD_LOCAL_DEFINED 0
# endif
#endif
