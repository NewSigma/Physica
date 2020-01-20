#ifndef _Physica_C_Utils_CUH
#define _Physica_C_Utils_CUH

#define CHECK(call)                                                                                        \
{                                                                                                          \
    const cudaError_t e = call;                                                                            \
    if (e != cudaSuccess)                                                                                  \
    {                                                                                                      \
		fprintf(stderr, "%s: %s (%s: %d)", cudaGetErrorName(e), cudaGetErrorString(e), __FILE__, __LINE__);\
		exit(1);                                                                                           \
    }                                                                                                      \
}

#endif