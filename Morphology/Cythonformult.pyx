cimport cython
cimport numpy as np
from libc.math cimport sqrt

ctypedef np.float64_t DTYPE_F64_t
ctypedef np.int64_t DTYPE_INT64_t 
ctypedef np.float32_t DTYPE_F32_t

cpdef potentialenergycalc(DTYPE_F64_t[:,:] mmc_image, 
                          DTYPE_INT64_t[:,:] mmc_image_mask):
    
    cdef DTYPE_F64_t potential_energy = 0.0
    cdef DTYPE_INT64_t i1
    cdef DTYPE_INT64_t j1
    cdef DTYPE_INT64_t i2
    cdef DTYPE_INT64_t j2
    cdef DTYPE_INT64_t N = mmc_image.shape[0]
    cdef DTYPE_F64_t r

    for i1 in range(N):
        for j1 in range(N):
             if mmc_image_mask[i1,j1] == 0:
                for i2 in range(i1 + 1, N):
                    for j2 in range(j1 + 1, N):
                         if mmc_image_mask[i2,j2] == 0:
                            r = sqrt((j1 - j2)**2 + (i1 - i2)**2)
                            potential_energy += (mmc_image[i1, j1] * mmc_image[i2, j2]) / r

    return potential_energy
