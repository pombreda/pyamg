#ifndef LINALG_H
#define LINALG_H

#include <math.h>
#include <limits>

/*******************************************************************
 * Overloaded routines for real arithmetic for int, float and double
 *******************************************************************/

/* Sign-of Function overloaded for int, float and double
 * signof(x) =  1 if x > 0
 * signof(x) = -1 if x < 0
 * signof(0) =  1 if x = 0
 */
inline int signof(int a) { return (a<0 ? -1 : 1); }
inline float signof(float a) { return (a<0.0 ? -1.0 : 1.0); }
inline double signof(double a) { return (a<0.0 ? -1.0 : 1.0); }



/*******************************************************************
 *         Overloaded routines for complex arithmetic for  
 *         pyamg's complex class, float and double
 *******************************************************************/

/*
 * Return the complex conjugate of a number
 */
inline float conjugate(const float& x)
    { return x; }
inline double conjugate(const double& x)
    { return x; }
inline npy_cfloat_wrapper conjugate(const npy_cfloat_wrapper& x)
    { return npy_cfloat_wrapper(x.real, -x.imag); }
inline npy_cdouble_wrapper conjugate(const npy_cdouble_wrapper& x)
    { return npy_cdouble_wrapper(x.real, -x.imag); }

/* 
 * Return the real part of a number
 */
inline float real(const float& x)
    { return x; }
inline double real(const double& x)
    { return x; }
inline float real(const npy_cfloat_wrapper& x)
    { return x.real; }
inline double real(const npy_cdouble_wrapper& x)
    { return x.real; }

/*
 * Return the imaginary part of a number
 */
inline float imag(const float& x)
    { return 0.0; }
inline double imag(const double& x)
    { return 0.0; }
inline float imag(const npy_cfloat_wrapper& x)
    { return x.imag; }
inline double imag(const npy_cdouble_wrapper& x)
    { return x.imag; }

/* 
 * Return the norm, i.e. the magnitude, of a single number
 */
inline float mynorm(const float& x)
    { return fabs(x); }
inline double mynorm(const double& x)
    { return fabs(x); }
inline float mynorm(const npy_cfloat_wrapper& x)
    { return sqrt(x.real*x.real + x.imag*x.imag); }
inline double mynorm(const npy_cdouble_wrapper& x)
    { return sqrt(x.real*x.real + x.imag*x.imag); }

/* 
 * Return the norm squared of a single number, i.e.  save a square root
 */
inline float mynormsq(const float& x)
    { return (x*x); }
inline double mynormsq(const double& x)
    { return (x*x); }
inline float mynormsq(const npy_cfloat_wrapper& x)
    { return (x.real*x.real + x.imag*x.imag); }
inline double mynormsq(const npy_cdouble_wrapper& x)
    { return (x.real*x.real + x.imag*x.imag); }



/*******************************************************************
 *              Dense Linear Algebra Routines
 *      templated for pyamg's complex class, float and double
 *******************************************************************/

/* dot(x, y, n)
 *
 * Parameters
 * ----------
 * x : {float|complex array}
 *      n-vector
 * y : {float|complex array}
 *      n-vector
 * n : {int}
 *      size of x and y
 *
 * Return
 * ------
 * conjuate(x).T y
 *
 */
template<class I, class T>
inline T dot_prod(const T x[], const T y[], const I n)
{
    T sum = 0.0;
    for( I i = 0; i < n; i++)
    {   sum += conjugate(x[i])*y[i]; }
    return sum;
}


/* norm(x, n)
 *
 * Parameters
 * ----------
 * x : {float|complex array}
 *      n-vector
 * n : {int}
 *      size of x and y
 * normx : {scalar}
 *      output value
 *
 * Return
 * ------
 * normx = sqrt( <x, x> )
 *
 */
template<class I, class T, class F>
inline void norm(const T x[], const I n, F &normx)
{
    normx = sqrt(real(dot_prod(x,x,n)));
}


/* axpy(x, y, alpha, n)
 *
 * Parameters
 * ----------
 * x : {float|complex array}
 *      n-vector
 * y : {float|complex array}
 *      n-vector
 * n : {int}
 *      size of x and y
 * alpha : {scalar}
 *      
 * Return
 * ------
 * x = x + alpha*y
 *
 */
template<class I, class T>
inline void axpy(T x[], const T y[], const T alpha, const I n)
{
    for( I i = 0; i < n; i++)
    {   x[i] += alpha*y[i]; }
}


/* Transpose Ax by overwriting Bx
 * 
 * Parameters
 * ----------
 * Ax : {float|complex array}
 *      m x n dense array
 * Bx : {float|complex array}
 *      m x n dense array
 * m,n : {int}
 *      Dimensions of Ax and Bx
 *
 * Return
 * ------
 * Bx is overwritten with the transpose of Ax
 *
 * Notes
 * -----
 * There is a fair amount of hardcoding to make this routine very 
 * fast for small (<10) square matrices, although it works for general 
 * m x n matrices.
 *
 */
template<class I, class T>
inline void transpose(const T Ax[], T Bx[], const I m, const I n)
{
    // Almost all uses of this function are for 
    // m==n, m,n<10.  Hence the attempts at speed.
    
    //Hard code the smallest examples for speed
    if( (m==1) && (n==1))
    {   Bx[0] = Ax[0]; }
    else if( (m==2) && (n==2))
    {
        Bx[0] = Ax[0];
        Bx[1] = Ax[2];
        Bx[2] = Ax[1];
        Bx[3] = Ax[3];
    }
    else if( (m==3) && (n==3))
    {
        Bx[0] = Ax[0];
        Bx[1] = Ax[3];
        Bx[2] = Ax[6];
        Bx[3] = Ax[1];
        Bx[4] = Ax[4];
        Bx[5] = Ax[7];
        Bx[6] = Ax[2];
        Bx[7] = Ax[5];
        Bx[8] = Ax[8];
    }
    // Do some hard coding for the other common examples
    else if ( (m == n) && (m < 11))
    {
        I j = 0;
        for(I i = 0; i < m*m; i+=m)
        {
            if(m == 4)
            {   Bx[i] = Ax[j]; Bx[i+1] = Ax[j+4]; Bx[i+2] = Ax[j+8]; Bx[i+3] = Ax[j+12]; }
            if(m == 5)
            {   Bx[i] = Ax[j]; Bx[i+1] = Ax[j+5]; Bx[i+2] = Ax[j+10]; Bx[i+3] = Ax[j+15]; 
                Bx[i+4] = Ax[j+20]; }
            if(m == 6)
            {   Bx[i] = Ax[j]; Bx[i+1] = Ax[j+6]; Bx[i+2] = Ax[j+12]; Bx[i+3] = Ax[j+18]; 
                Bx[i+4] = Ax[j+24];  Bx[i+5] = Ax[j+30]; }
            if(m == 7)
            {   Bx[i] = Ax[j]; Bx[i+1] = Ax[j+7]; Bx[i+2] = Ax[j+14]; Bx[i+3] = Ax[j+21]; 
                Bx[i+4] = Ax[j+28];  Bx[i+5] = Ax[j+35]; Bx[i+6] = Ax[j+42]; }
            if(m == 8)
            {   Bx[i] = Ax[j]; Bx[i+1] = Ax[j+8]; Bx[i+2] = Ax[j+16]; Bx[i+3] = Ax[j+24]; 
                Bx[i+4] = Ax[j+32];  Bx[i+5] = Ax[j+40]; Bx[i+6] = Ax[j+48]; Bx[i+7] = Ax[j+56]; }
            if(m == 9)
            {   Bx[i] = Ax[j]; Bx[i+1] = Ax[j+9]; Bx[i+2] = Ax[j+18]; Bx[i+3] = Ax[j+27]; 
                Bx[i+4] = Ax[j+36];  Bx[i+5] = Ax[j+45]; Bx[i+6] = Ax[j+54]; Bx[i+7] = Ax[j+63]; 
                Bx[i+8] = Ax[j+72];}
            if(m == 10)
            {   Bx[i] = Ax[j]; Bx[i+1] = Ax[j+10]; Bx[i+2] = Ax[j+20]; Bx[i+3] = Ax[j+30]; 
                Bx[i+4] = Ax[j+40];  Bx[i+5] = Ax[j+50]; Bx[i+6] = Ax[j+60]; Bx[i+7] = Ax[j+70]; 
                Bx[i+8] = Ax[j+80]; Bx[i+9] = Ax[j+90];}

            j++;
        }
    }
    // Finally, the general case
    else
    {
        I Bcounter = 0;
        for(I i = 0; i < n; i++)
        {
            I Acounter = i;
            for(I j = 0; j < m; j++)
            {
                //B[i,j] = A[j,i]
                Bx[Bcounter] = Ax[Acounter];
                Bcounter++;
                Acounter+=n;
            }
        }
    }

    return;
}


/* Calculate Ax*Bx = S
 *
 * Parameters
 * ----------
 * Ax : {float|complex array} 
 *      Stored in row major
 * Arows : {int}
 *      Number of rows of A
 * Acols : {int}
 *      Number of columns of A
 * Atrans : {char}
 *      Not Used
 * Bx : {float|complex array} 
 *      Stored in col major
 * Brows : {int}
 *      Number of rows of B
 * Bcols : {int}
 *      Number of columns of B
 * Btrans : {char}
 *      Not Used
 * Sx : {float|complex array} 
 *      Output array, Contents are overwitten
 * Srows : {int}
 *      Number of rows of S
 * Scols : {int}
 *      Number of columns of S
 * Strans : {char}
 *      'T' gives S in col major
 *      'F' gives S in row major
 *
 * Return
 * ------
 * Sx = Ax*Bx in column or row major, depending on Strans.
 * Contents of Sx are overwritten
 *
 * Notes
 * -----
 * Naively calculates S(i,j) = A(i,:) B(:,j) by looping over the rows of A
 * and the columns of B.
 *
 */
template<class I, class T>
void gemm(const T Ax[], const I Arows, const I Acols, const char Atrans, 
          const T Bx[], const I Brows, const I Bcols, const char Btrans, 
          T Sx[], const I Srows, const I Scols, const char Strans)
{
    //Add checks for dimensions, but leaving them out speeds things up
    //Add functionality for transposes

    if(Strans == 'T')
    {
        I s_counter = 0; I a_counter =0; I b_counter =0; I a_start = 0;
        for(I i = 0; i < Arows; i++)
        {
            s_counter = i;
            b_counter = 0; 
            for(I j = 0; j < Bcols; j++)
            {
                Sx[s_counter] = 0.0;
                a_counter = a_start;
                for(I k = 0; k < Brows; k++)
                {
                    //S[i,j] += Ax[i,k]*B[k,j]
                    Sx[s_counter] += Ax[a_counter]*Bx[b_counter];
                    a_counter++; b_counter++;
                }
                s_counter+=Scols;
            }
            a_start += Acols;
        }
    }
    else if(Strans == 'F')
    {
        I s_counter = 0; I a_counter =0; I b_counter =0; I a_start = 0;
        for(I i = 0; i < Arows; i++)
        {
            b_counter = 0; 
            for(I j = 0; j < Bcols; j++)
            {
                Sx[s_counter] = 0.0;
                a_counter = a_start;
                for(I k = 0; k < Brows; k++)
                {
                    //S[i,j] += A[i,k]*B[k,j]
                    Sx[s_counter] += Ax[a_counter]*Bx[b_counter];
                    a_counter++; b_counter++;
                }
                s_counter++;
            }
            a_start += Acols;
        }
    }
}


/*
 * Compute the SVD of a matrix, Ax, using the Jacobi method.
 * Compute Ax = U S V.H
 *
 * Parameters
 * ----------
 * Ax : {float|complex array}
 *      m x n dense matrix, stored in col major form
 * U : {float|complex array}
 *      m x n dense matrix initialized to 0.0 
 *      Passed in as Tx
 * V : {float|complex array}
 *      n x n dense matrix initialized to 0.0 
 *      Passed in as Bx
 * S : {float|complex array}
 *      n x 1 dense matrix initialized to 0.0 
 *      Passed in as Sx
 * m,n : {int} 
 *      Dimensions of Ax, m > n.
 *
 * Return
 * ------
 * Returns Ax = U S V.H
 * U, V, S are modified in place
 * 
 * V : {array}
 *      Orthogonal n x n matrix, V, stored in col major
 * U : {array}
 *      Orthogonal m x nmatrix, U, stored in col major
 * S : {array}
 *      Singular values
 * int : {int}
 *      Function return value, 
 *      -1:  error
 *      0:  successful
 *      1:  did not converge
 *
 * Notes
 * -----
 * The Jacobi method is used to compute the SVD.  Conceptually, 
 * the Jacobi method applies successive Jacobi rotations, Q_i to
 * the system, Q_i^H Ax.H Ax Q_i.  Despite the normal equations 
 * appearing here, the actual method can be quite accurate.  
 * However, the method is slower than Golub-Reinsch for all 
 * but very small matrices.
 *
 * References
 * ----------
 * De Rijk, "A One-Sided Jacobi Algorithm for computing the singular value 
 * decomposition on a vector computer", SIAM J Sci and Statistical Comp,
 * Vol 10, No 2, p 359-371, March 1989.
 *
 */

template<class I, class T, class F>
I svd_jacobi (const T Ax[], T Tx[], T Bx[], F Sx[], const I m, const I n)
{
    // Not implemented for m < n matrices
    if( m < n)
    {   return -1; }

    // Rename
    const T * A = Ax;
    T * U = Tx;
    T * V = Bx;
    F * S = Sx;
    
    // Hard code fast 1x1 SVD
    if ( (n==1) && (m==1) )
    {
        F normA = mynorm(A[0]);

        V[0] = 1.0;
        S[0] = normA; 
        if(normA == 0.0)
        {   U[0] = 1.0; }
        else
        {   U[0] = A[0]/normA; }
        
        return 0;
    }
  
    // Workspace
    I i, j, k;
    I nsq = n*n;
    F normx;

    // Initialize the rotation counter and the sweep counter.
    I count = 1;
    I sweep = 0;

    // Always do at least 12 sweeps
    I sweepmax = std::max(8*n, 12);

    F tolerance = 10.0*m*std::numeric_limits<F>::epsilon();

    // Set V to the identity matrix
    for(i = 0; i < nsq; i++)
    {   V[i] = 0.0;}
    for(i = 0; i < nsq; i+= (n+1) )
    {   V[i] = 1.0;}

    // Copy A to U, note that the stop address &(A[m*n]) 
    // should go one past the final element to be copied
    std::copy(&(A[0]), &(A[m*n]), &(U[0]));

    // Store the column error estimates in S, for use during the orthogonalization
    for (j = 0; j < n; j++)
    {
        // S[j] = eps*norm(A[:,j])
        norm(&(U[j*m]), m, normx);
        S[j] = std::numeric_limits<F>::epsilon()*normx;
    }
  
    // Orthogonalize U by plane rotations.
    while( (count > 0) && (sweep <= sweepmax) )
    {
        // Initialize rotation counter.
        count = n*(n - 1)/2;
        I jm = 0;
        I jn = 0;

        for(j=0; j < n - 1; j++)
        {
            I km = (j+1)*m;
            I kn = (j+1)*n;

            for (k = j + 1; k < n; k++)
            {
                F cos, abserr_a, abserr_b;
                T sin, neg_conj_sin;
                I sorted, orthog, noisya, noisyb;

                F a; norm(&(U[jm]), m, a);              // || U[:,j] ||
                F b; norm(&(U[km]), m, b);              // || U[:,k] ||
                T d = dot_prod(&(U[jm]), &(U[km]), m);  // <U[:,j], U[:,k]>
                F norm_d = mynorm(d);

                // test for columns j,k orthogonal, or dominant errors 
                abserr_a = S[j];
                abserr_b = S[k];

                sorted = (a >= b);
                orthog = (norm_d <= tolerance*a*b);
                
                // Test to see if col a or b has become noise
                noisya = (a < abserr_a);
                noisyb = (b < abserr_b);

                // no need for rotations
                if(sorted && (orthog || noisya || noisyb))
                {
                    // if count ever = 0, then everything is sorted and orthogonal 
                    // (or possibly just noise)
                    count--;
                }
                
                // swap cols ||     Handle 0 matrix case
                else if(!sorted   || ( (norm_d == 0.0) && (a==b)  ) )
                {
                    // Apply rotation matrix,
                    // [ 0.0  1.0 ]
                    // [-1.0  0.0 ]
                    // Basically, swap columns in U and V with one sign flip
                    
                    S[j] = abserr_b;
                    S[k] = abserr_a;

                    // apply rotation by right multiplication to U
                    I koffset = km;
                    for(I joffset = jm; joffset < (jm + m); joffset++)
                    {
                        // for i = 0:m-1
                        //   U[i,j] =   0.0*U[i,j] - 1.0*U[i,k]
                        //   U[i,k] =   1.0*U[i,j] + 0.0*U[i,k]
                        const T Uij = U[joffset];
                        const T Uik = U[koffset];
                        U[joffset] = -Uik; 
                        U[koffset] =  Uij;
                        koffset++;
                    }
    
                    // apply rotation by right multiplication to V
                    koffset = kn;
                    for(I joffset = jn; joffset < (jn + n); joffset++)
                    {
                        // for i = 0:n-1
                        //   V[i,j] =   0.0*V[i,j] - 1.0*V[i,k]
                        //   V[i,k] =   1.0*V[i,j] + 0.0*V[i,k]
                        const T Vij = V[joffset];
                        const T Vik = V[koffset];
                        V[joffset] = -Vik; 
                        V[koffset] =  Vij; 
                        koffset++;
                    }
                }

                // Carry out Jacobi Rotations to orgthogonalize column's j and k in U
                else
                {
                    // calculate rotation angles for 
                    // jacobi_rot = [cos          sin]
                    //              [-conj(sin)   cos]
                    F tau = (b*b - a*a)/(2.0*norm_d);
                    F t = signof(tau)/(fabs(tau) + sqrt(1.0 + tau*tau));
                    cos = 1.0/(sqrt(1.0 + t*t));
                    sin = d*(t*cos/norm_d);
                    neg_conj_sin = conjugate(sin)*-1.0;
                
                    F norm_sin = mynorm(sin);
                    S[j] = fabs(cos)*abserr_a + norm_sin*abserr_b;
                    S[k] =  norm_sin*abserr_a + fabs(cos)*abserr_b;

                    // apply rotation by right multiplication to U
                    I koffset = km;
                    for(I joffset = jm; joffset < (jm + m); joffset++)
                    {
                        // for i = 0:m-1
                        //   U[i,j] =   cos*U[i,j] + -conj(sin)*U[i,k]
                        //   U[i,k] =   sin*U[i,j] +       cos*U[i,k]
                        const T Uij = U[joffset];
                        const T Uik = U[koffset];
                        U[joffset] = Uij*cos + neg_conj_sin*Uik; 
                        U[koffset] = sin*Uij + Uik*cos;
                        koffset++;
                    }
    
                    // apply rotation by right multiplication to V
                    koffset = kn;
                    for(I joffset = jn; joffset < (jn + n); joffset++)
                    {
                        // for i = 0:n-1
                        //   V[i,j] =   cos*V[i,j] + -conj(sin)*V[i,k]
                        //   V[i,k] =   sin*V[i,j] +       cos*V[i,k]
                        const T Vij = V[joffset];
                        const T Vik = V[koffset];
                        V[joffset] = Vij*cos + neg_conj_sin*Vik; 
                        V[koffset] = sin*Vij + Vik*cos; 
                        koffset++;
                    }
                }

                km += m;
                kn += n;
            } // end k loop

            jm += m;
            jn += n;
        } // end j loop
        
        //Sweep completed.
        sweep++;

    }// end while loop

    // Orthogonalization complete. Compute singular values.
    F sigma_tol=0.0;
    I Uoffset = 0;
    I iszero = n;
    for (j = 0; j < n; j++)
    {
        F curr_norm;
        norm(&(U[Uoffset]), m, curr_norm);              // || U[:,j] ||
        
        if(j == 0)
        {
            F alpha = std::max(2000.0, 100.0*std::max(m,n)*std::numeric_limits<F>::epsilon()*curr_norm);
            sigma_tol = alpha*std::numeric_limits<F>::epsilon();
        }

        // Determine if singular value is zero
        if( curr_norm <= sigma_tol )
        {   
            iszero--;                               // detect all zero matrix
            S[j] = 0.0;                             // singular
            for(i = Uoffset; i < (Uoffset + m); i++)
            {   U[i] = 0.0; }                       // annihilate U[:,j]
        }
        else
        {
            S[j] = curr_norm;                       // non-singular
            for(i = Uoffset; i < (Uoffset + m); i++)
            {   U[i] = U[i]/curr_norm; }            // normalize column U[:,j]
        }

        Uoffset += m;
    }

    if(iszero == 0)
    {
        // Set U and V to the identity matrix
        for(i = 0; i < nsq; i++)
        {   V[i] = 0.0;}
        for(i = 0; i < nsq; i+= (n+1) )
        {   V[i] = 1.0;}
        
        // U is already 0.0
        for(i = 0; i < n*m; i+= (m+1) )
        {   U[i] = 1.0;}

        return 0;
    }

    if(count > 0)
    {
        // reached sweep limit, i.e. did not converge
        return 1;
    }

    return 0;
}
 

/*
 * Solve a system with the SVD, i.e. use a robust Moore-Penrose 
 * Pseudoinverse to multiply the RHS
 * 
 * Parameters
 * ----------
 * A : {float|complex array} 
 *      m x n dense column major array, m>n
 * m,n : {int}
 *      Dimensions of A, m > n
 * b : {float|complex array}
 *      RHS, m-vector
 * sing_vals : {float array}
 *      Holds the singular values upon return
 * work : {float|complex array}
 *      worksize array for temporary space for routine 
 * worksize : {int}
 *      must be > m*n + n
 *
 * Return
 * ------
 * A^{-1} b replaces b
 * sing_vals holds the singular values
 *
 * Notes
 * -----
 * forcing preallocation of sing_vals and work, allows for 
 * efficient multiple calls to this routine
 *
 */
template<class I, class T, class F>
void svd_solve( T Ax[], I m, I n, T b[], F sing_vals[], T work[], I work_size)
{
    I mn = m*n;
    // Rename
    T * U = &(work[0]);
    T * V = &(work[mn]);
    T * x = &(work[2*mn]);
    const char trans = 'F';
    
    // calculate SVD
    svd_jacobi(&(Ax[0]), &(U[0]), &(V[0]), &(sing_vals[0]), n, n);
        
    // Forming conjugate(U.T) in row major requires just
    // conjugating the current entries of U in col major
    for(I i = 0; i < m*n; i++) 
    {   U[i] = conjugate(U[i]); }

    // A^{-1} b = V*Sinv*U.H*b, in 3 steps
    // Step 1, U.H*b
    gemm(&(U[0]), n, n, trans, &(b[0]), n, 1, trans,  
         &(x[0]), n, 1, trans);

    // Setp 2, scale x by Sinv
    for(I j = 0; j < n; j++)
    {
        if(sing_vals[j] != 0.0)
        {   x[j] = x[j]/sing_vals[j]; }
        else
        {   x[j] = 0.0; }
    }

    // Step 3, multiply by V
    // transpose V so that it is in row major for gemm
    transpose(&(V[0]), &(U[0]), n, n);
    gemm(&(U[0]), n, n, trans, &(x[0]), n, 1, trans,  
         &(b[0]), n, 1, trans);
    
    return;
}

 /* Replace each block of A with a Moore-Penrose pseudoinverse of that block.
 * Routine is designed to invert many small matrices at once.
 * Parameters
 * ----------
 * Ax : {float|complex array}  
 *      (m, n, n) array, assumed to be "raveled" and in row major form
 * m,n : int
 *      dimensions of Ax
 * TransA : char
 *      'T' or 'F'.  Decides whether to transpose each nxn block
 *      of A before inverting.  If using Python array, should be 'T'.  
 * 
 * Return
 * ------
 * Ax : {array}
 *      Ax is modified in place with the pseduo-inverse replacing each
 *      block of Ax.  Ax is returned in row-major form for Python
 *
 * Notes
 * -----
 * This routine is designed to be called once for a large m.
 * Calling this routine repeatably would not be efficient.
 *
 * This function offers substantial speedup over native Python
 * code for many small matrices, e.g. 5x5 and 10x10.  Tests have
 * indicated that matrices larger than 27x27 are faster if done 
 * in native Python.
 *
 * Examples
 * --------
 * >>> from pyamg.amg_core import pinv_array
 * >>> from scipy import arange, ones, array, dot
 * >>> A = array([arange(1,5, dtype=float).reshape(2,2), ones((2,2),dtype=float)])
 * >>> Ac = A.copy()
 * >>> pinv_array(A, 2, 2, 'T')
 * >>> print "Multiplication By Inverse\n" + str(dot(A[0], Ac[0]))
 * >>> print "Multiplication by PseudoInverse\n" + str(dot(Ac[1], dot(A[1], Ac[1])))
 * >>>
 * >>> A = Ac.copy()
 * >>> pinv_array(A,2,2,'F')
 * >>> print "Changing flag to \'F\' results in different Inverse\n" + str(dot(A[0], Ac[0]))
 * >>> print "A holds the inverse of the transpose\n" + str(dot(A[0], Ac[0].T))
 *
 */
template<class I, class T, class F>
void pinv_array(T Ax[], const I m, const I n, const char TransA)
{
    I nsq = n*n;
    I Acounter = 0;
    T * Tran = new T[nsq];
    T * U = new T[nsq];
    T * V = new T[nsq];
    T * SinvUh = new T[nsq];
    F * S = new F[n];
    const char t = 'F';

    for(I i = 0; i < m; i++)
    {
        if(TransA == 'T')
        {   // transpose block of A so that it is in col major for SVD
            transpose(&(Ax[Acounter]), &(Tran[0]), n, n); 
            
            // calculate SVD
            svd_jacobi(&(Tran[0]), &(U[0]), &(V[0]), &(S[0]), n, n);
        }
        else
        {
            // calculate SVD
            svd_jacobi(&(Ax[Acounter]), &(U[0]), &(V[0]), &(S[0]), n, n);
        }

        // invert S
        for(I j = 0; j < n; j++)
        {
            if(S[j] != 0.0)
            {   S[j] = 1.0/S[j]; }
        }
        
        // Sinv*conjugate(U.T), stored in column major form
        I counter = 0;
        for(I j = 0; j < n; j++)        // col of Uh
        {
            I Uoffset = j;
            for(I k = 0; k < n; k++)    // row of Uh
            {
                //Sinv(j) * conj(U(j,k)) ==> SinvUh(k,j)
                SinvUh[counter] = conjugate(U[Uoffset])*S[k];
                counter++;
                Uoffset+=n;
            }
        }

        // transpose V so that it is in row major for gemm
        transpose(&(V[0]), &(Tran[0]), n, n);

        // A^{-1} = V*SinvUh
        gemm(&(Tran[0]), n, n, t, &(SinvUh[0]), n, n, t,  
             &(Ax[Acounter]), n, n, t);

        Acounter += nsq;
    }
    
    delete[] Tran;
    delete[] U;
    delete[] V;
    delete[] S;
    delete[] SinvUh;
    
    return;
}


#endif
