#ifndef FORCUDA_H
#define FORCUDA_H

#if defined(__CUDACC__)
#define FORCUDA __host__ __device__
#else
#define FORCUDA
#endif

#endif // FORCUDA_H