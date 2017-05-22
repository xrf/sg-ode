/* Compatibility header for `thread_local`. */
#if !defined SG_THREAD_LOCAL_DEFINED
# error "must not include end header without matching begin"
#elif SG_THREAD_LOCAL_DEFINED
# undef thread_local
#endif
#undef SG_THREAD_LOCAL_DEFINED
