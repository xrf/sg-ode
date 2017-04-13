#if defined _WIN32 || defined __CYGWIN__
# ifdef SG_BUILD
#  define SG_EXTERN __declspec(dllexport)
# else
#  define SG_EXTERN __declspec(dllimport)
# endif
#elif __GNUC__ >= 4
# define SG_EXTERN __attribute__ ((visibility ("default")))
#else
# define SG_EXTERN
#endif
