/* Compatibility header for `restrict`. */
#if !defined SG_RESTRICT_DEFINED
# error "must not include end header without matching begin"
#elif SG_RESTRICT_DEFINED
# undef restrict
#endif
#undef SG_RESTRICT_DEFINED
