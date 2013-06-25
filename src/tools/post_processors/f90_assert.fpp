

#if ( defined(NDEBUG) || defined(DARWIN_ABSOFT_COMPILER_WORKAROUND) )
# define ASSERT(x) !! assert( x )
#else
# define ASSERT(x) if(.not.(x)) call f90_assert(__FILE__,__LINE__)
#endif

#define INSIST(x) if(.not.(x)) call f90_assert(__FILE__,__LINE__)
