/* Compatibility header for `inline`. */
#if !defined SG_INLINE_DEFINED
# error "must not include end header without matching begin"
#elif SG_INLINE_DEFINED
# undef inline
#endif
#undef SG_INLINE_DEFINED
