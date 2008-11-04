#ifndef RELAXATION_H
#define RELAXATION_H

#include "linalg.h"

template<class I, class T, class F>
void gauss_seidel(const I Ap[], 
                  const I Aj[], 
                  const T Ax[],
                        T  x[],
                  const T  b[],
                  const I row_start,
                  const I row_stop,
                  const I row_step)
{
    for(I i = row_start; i != row_stop; i += row_step) {
        I start = Ap[i];
        I end   = Ap[i+1];
        T rsum = 0;
        T diag = 0;
        
        for(I jj = start; jj < end; jj++){
            I j = Aj[jj];
            if (i == j)
                diag  = Ax[jj];
            else
                rsum += Ax[jj]*x[j];
        }
        
        //TODO raise error? inform user?
        if (diag != 0){
            x[i] = (b[i] - rsum)/diag;
        }
    }
}


template<class I, class T, class F>
void block_gauss_seidel(const I Ap[], 
                        const I Aj[], 
                        const T Ax[],
                              T  x[],
                        const T  b[],
                        const I row_start,
                        const I row_stop,
                        const I row_step,
                        const I blocksize)
{
    const I B2 = blocksize * blocksize;
    
    for(I i = row_start; i != row_stop; i += row_step) {
        I start = Ap[i];
        I end   = Ap[i+1];

        for(I bi = 0; bi < blocksize; bi++){
            T rsum = 0;
            T diag = 0;

            for(I jj = start; jj < end; jj++){
                I j = Aj[jj];
                const T * block_row = Ax + B2*jj + blocksize*bi;
                const T * block_x   = x + blocksize * j;

                if (i == j){
                    //diagonal block
                    diag = block_row[bi];
                    for(I bj = 0; bj < bi; bj++){
                        rsum += block_row[bj] * block_x[bj];
                    }
                    for(I bj = bi+1; bj < blocksize; bj++){
                        rsum += block_row[bj] * block_x[bj];
                    }
                } else {
                    for(I bj = 0; bj < blocksize; bj++){
                        rsum += block_row[bj] * block_x[bj];
                    }
                }
            }

            //TODO raise error? inform user?
            if (diag != 0){
                x[blocksize*i + bi] = (b[blocksize*i + bi] - rsum)/diag;
            }
        }
    }
}

template<class I, class T, class F>
void jacobi(const I Ap[], 
            const I Aj[], 
            const T Ax[],
                  T  x[],
            const T  b[],
                  T temp[],
            const I row_start,
            const I row_stop,
            const I row_step,
            const T omega[])
{
    T one = 1.0;
    T omega2 = omega[0];

    for(I i = row_start; i != row_stop; i += row_step) {
        temp[i] = x[i];
    }
    
    for(I i = row_start; i != row_stop; i += row_step) {
        I start = Ap[i];
        I end   = Ap[i+1];
        T rsum = 0;
        T diag = 0;
        
        for(I jj = start; jj < end; jj++){
            I j = Aj[jj];
            if (i == j)
                diag  = Ax[jj];
            else
                rsum += Ax[jj]*temp[j];
        }
        
        //TODO raise error? inform user?
        if (diag != 0){ 
            x[i] = (one - omega2) * temp[i] + omega2 * ((b[i] - rsum)/diag);
        }
    }
}

//
// Guass Seidel Indexed will relax a specific index field Id
//
template<class I, class T, class F>
void gauss_seidel_indexed(const I Ap[], 
                          const I Aj[], 
                          const T Ax[],
                                T  x[],
                          const T  b[],
                          const I Id[],
                          const I row_start,
                          const I row_stop,
                          const I row_step)
{
  for(I i = row_start; i != row_stop; i += row_step) {
    I inew = Id[i];
    I start = Ap[inew];
    I end   = Ap[inew+1];
    T rsum  = 0;
    T diag  = 0;

    for(I jj = start; jj < end; ++jj){
      I j = Aj[jj];
      if (inew == j){
        diag = Ax[jj];
      }
      else{
        rsum += Ax[jj]*x[j];
      }
    }

    if (diag != 0){
      x[inew] = (b[inew] - rsum)/diag;
    }
  }
}

//
// weighted-Jacobi on the normal equations, i.e. Kaczmarz type iteration
//
template<class I, class T, class F>
void kaczmarz_jacobi(const I Ap[], 
                     const I Aj[], 
                     const T Ax[],
                           T  x[],
                     const T  b[],
                     const T Tx[],
                           T temp[],
                     const I row_start,
                     const I row_stop,
                     const I row_step,
                     const T omega[])
{
    //rename
    const T * delta = Tx;
    const T omega2 = omega[0];

    for(I i = row_start; i < row_stop; i+=row_step)
    {   temp[i] = 0.0; }

    for(I i = row_start; i < row_stop; i+=row_step)
    {
        I start = Ap[i];
        I end   = Ap[i+1];
        for(I j = start; j < end; j++)
        {   temp[Aj[j]] += omega2*conjugate(Ax[j])*delta[i]; }
    }

    for(I i = row_start; i < row_stop; i+=row_step)
    {   x[i] += temp[i]; }
}

//
// Gauss Seidel on the normal equations, i.e. Kaczmarz type iteration
//
template<class I, class T, class F>
void kaczmarz_gauss_seidel(const I Ap[], 
                     const I Aj[], 
                     const T Ax[],
                           T  x[],
                     const T  b[],
                     const I row_start,
                     const I row_stop,
                     const I row_step,
                     const T Tx[])
{
    //rename
    const T * D_inv = Tx;
    
    for(I i = row_start; i != row_stop; i+=row_step)
    {
        I start = Ap[i];
        I end   = Ap[i+1];
        
        //First calculate "delta", the scaled residual term
        T delta = 0.0;
        for(I j = start; j < end; j++)
        {   delta += Ax[j]*x[Aj[j]]; }
        delta = (b[i] - delta)*D_inv[i];

        for(I j = start; j < end; j++)
        {   x[Aj[j]] += conjugate(Ax[j])*delta; }
    }

}

#endif
