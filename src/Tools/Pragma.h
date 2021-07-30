#ifndef PRAGMA
#define PRAGMA

#define to_string(a) #a

#if defined(__clang__)
    // # pragma clang loop unroll(n)
    #define UNROLL(n) \
    _Pragma( to_string(clang loop unroll(full)))
    #define UNROLL_S(n) \
    _Pragma( to_string(clang loop unroll(full)))
#elif defined (__FUJITSU)
    // #pragma loop fullunroll_pre_simd n
    #define UNROLL(n) \
    _Pragma( to_string(loop fullunroll_pre_simd n))
    #define UNROLL_S(n) \
    _Pragma( to_string(loop fullunroll_pre_simd n))
#elif defined(__GNUC__)
    // #pragma GCC unroll (n)
	#define UNROLL(n) \
    _Pragma( to_string(GCC unroll (n)))
	#define UNROLL_S(n) \
    _Pragma( to_string(GCC unroll (n)))
#else
    //#pragma unroll(5)
    #define UNROLL(n) \
    _Pragma( to_string( pragma unroll (n)))
    #define UNROLL_S(n) 
#endif


#endif
