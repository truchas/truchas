#ifdef NDEBUG
# define ASSERT(x) !! assert( x )
# define ASSERT_MSG(x,msg) !! assert( x )
#else
# define ASSERT(x) if(.not.(x)) call t2d_assert(__FILE__,__LINE__,'')
# define ASSERT_MSG(x,msg) if(.not.(x)) call t2d_assert(__FILE__,__LINE__,msg)
#endif

#define INSIST(x) if(.not.(x)) call t2d_assert(__FILE__,__LINE__,'')
#define INSIST_MSG(x,msg) if(.not.(x)) call t2d_assert(__FILE__,__LINE__,msg)
