__all__ =['approximate_spectral_radius','infinity_norm','diag_sparse',
          'hstack_csr','vstack_csr']

import numpy
import scipy
from numpy import fromfile, ascontiguousarray, mat, int32
from scipy import ravel, arange, concatenate, tile, asarray, sqrt, diff, \
                  rand, zeros, ones, empty, asmatrix, dot
from scipy.linalg import norm, eigvals
from scipy.sparse import isspmatrix, isspmatrix_csr, isspmatrix_csc, \
        isspmatrix_bsr, csr_matrix, csc_matrix, bsr_matrix, coo_matrix
from scipy.sparse.sputils import upcast


def approximate_spectral_radius(A,tol=0.1,maxiter=10,symmetric=None):
    """approximate the spectral radius of a matrix

    *Parameters*:
        A : dense or sparse matrix 
            E.g. csr_matrix, csc_matrix, ndarray, etc.

        tol : {scalar}
            Tolerance of approximation
            Currently unused

        maxiter : {integer}
            Maximum number of iterations to perform

        symmetric : {None,boolean}
            True  - if A is symmetric
                    Lanczos iteration is used (more efficient)
            False - if A is non-symmetric
                    Arnoldi iteration is used (less efficient)
            None  - symmetry of A unknown
                    Method chosen automatically (default)
    *Returns*:
        An approximation to the spectral radius of A (scalar value)

    """
    #from scipy.sandbox.arpack import eigen
    #return norm(eigen(A, k=1, ncv=min(10,A.shape[0]), which='LM', tol=tol, return_eigenvectors=False))
   
    if not isspmatrix(A):
        A = asmatrix(A) #convert dense arrays to matrix type
    
    if A.shape[0] != A.shape[1]:
        raise ValueError,'expected square matrix'

    maxiter = min(A.shape[0],maxiter)

    #TODO make method adaptive

    numpy.random.seed(0)  #make results deterministic

    v0  = rand(A.shape[1],1)
    v0 /= norm(v0)

    H  = zeros((maxiter+1,maxiter))
    V = [v0]

    #save past estimates
    #estimates = []

    for j in range(maxiter):
        w = A * V[-1]
   
        if symmetric:
            if j >= 1:
                H[j-1,j] = beta
                w -= beta * V[-2]

            alpha = dot(ravel(w),ravel(V[-1]))
            H[j,j] = alpha
            w -= alpha * V[-1]
            
            beta = norm(w)
            if (H[j+1,j] < 1e-10): break
            
            w /= beta
            H[j+1,j] = beta

            V.append(w)
            V = V[-2:] #retain only last two vectors
        else:
            #orthogonalize against Vs
            for i,v in enumerate(V):
                H[i,j] = dot(ravel(w),ravel(v))
                w -= H[i,j]*v
            H[j+1,j] = norm(w)
            if (H[j+1,j] < 1e-10): break
            
            w /= H[j+1,j] 
            V.append(w)
   
            # if upper 2x2 block of Hessenberg matrix H is almost symmetric,
            # and the user has not explicitly specified symmetric=False,
            # then switch to symmetric Lanczos algorithm
            if symmetric is not False and j == 1:
                if abs(H[1,0] - H[0,1]) < 1e-12:
                    symmetric = True
                    V = V[1:]
                    H[1,0] = H[0,1]
                    beta = H[2,1]
    
    return max([norm(x) for x in eigvals(H[:j+1,:j+1])])      



def infinity_norm(A):
    """
    Infinity norm of a sparse matrix (maximum absolute row sum).  This serves
    as an upper bound on spectral radius.
    """

    if isspmatrix_csr(A) or isspmatrix_csc(A):
        #avoid copying index and ptr arrays
        abs_A = A.__class__((abs(A.data),A.indices,A.indptr),shape=A.shape)
        return (abs_A * ones(A.shape[1],dtype=A.dtype)).max()
    else:
        return (abs(A) * ones(A.shape[1],dtype=A.dtype)).max()

def diag_sparse(A):
    """
    If A is a sparse matrix (e.g. csr_matrix or csc_matrix)
       - return the diagonal of A as an array

    Otherwise
       - return a csr_matrix with A on the diagonal
    """

    #TODO integrate into SciPy?
    if isspmatrix(A):
        return A.diagonal()
    else:
        return csr_matrix((asarray(A),arange(len(A)),arange(len(A)+1)),(len(A),len(A)))

def scale_rows(A,v,copy=True):
    from scipy.sparse.sparsetools import csr_scale_rows, bsr_scale_rows

    v = ravel(v)

    if isspmatrix_csr(A) or isspmatrix_bsr(A):
        M,N = A.shape
        if M != len(v):
            raise ValueError,'scale vector has incompatible shape'

        if copy:
            A = A.copy()
            A.data = asarray(A.data,dtype=upcast(A.dtype,v.dtype))
        else:
            v = asarray(v,dtype=A.dtype)

        if isspmatrix_csr(A):
            csr_scale_rows(M, N, A.indptr, A.indices, A.data, v)
        else:
            R,C = A.blocksize
            bsr_scale_rows(M/R, N/C, R, C, A.indptr, A.indices, ravel(A.data), v)

        return A
    elif isspmatrix_csc(A):
        return scale_columns(A.T,v)
    else:
        return scale_rows(csr_matrix(A),v)
        
def scale_columns(A,v,copy=True):
    from scipy.sparse.sparsetools import csr_scale_columns, bsr_scale_columns

    v = ravel(v)

    if isspmatrix_csr(A) or isspmatrix_bsr(A):
        M,N = A.shape
        if N != len(v):
            raise ValueError,'scale vector has incompatible shape'

        if copy:
            A = A.copy()
            A.data = asarray(A.data,dtype=upcast(A.dtype,v.dtype))
        else:
            v = asarray(v,dtype=A.dtype)

        if isspmatrix_csr(A):
            csr_scale_columns(M, N, A.indptr, A.indices, A.data, v)
        else:
            R,C = A.blocksize
            bsr_scale_columns(M/R, N/C, R, C, A.indptr, A.indices, ravel(A.data), v)

        return A
    elif isspmatrix_csc(A):
        return scale_rows(A.T,v)
    else:
        return scale_rows(csr_matrix(A),v)

def symmetric_rescaling(A,copy=True):
    if isspmatrix_csr(A) or isspmatrix_csc(A) or isspmatrix_bsr(A):
        if A.shape[0] != A.shape[1]:
            raise ValueError,'expected square matrix'

        D = diag_sparse(A)
        mask = D == 0

        D_sqrt = sqrt(abs(D))
        D_sqrt_inv = 1.0/D_sqrt
        D_sqrt_inv[mask] = 0

        DAD = scale_rows(A,D_sqrt_inv,copy=copy)
        DAD = scale_columns(DAD,D_sqrt_inv,copy=False)

        return D_sqrt,D_sqrt_inv,DAD

    else:
        return symmetric_rescaling(csr_matrix(A))


def hstack_csr(A,B):
    if not isspmatrix(A) or not isspmatrix(B):
        raise TypeError,'expected sparse matrix'

    if A.shape[0] != B.shape[0]:
        raise ValueError,'row dimensions must agree'

    A = A.tocoo()
    B = B.tocoo()
    I = concatenate((A.row,B.row))
    J = concatenate((A.col,B.col+A.shape[1]))
    V = concatenate((A.data,B.data))
    return coo_matrix((V,(I,J)),shape=(A.shape[0],A.shape[1]+B.shape[1])).tocsr()

def vstack_csr(A,B):
    #TODO OPTIMIZE THIS
    if not isspmatrix(A) or not isspmatrix(B):
        raise TypeError,'expected sparse matrix'

    if A.shape[1] != B.shape[1]:
        raise ValueError,'column dimensions must agree'

    A = A.tocoo()
    B = B.tocoo()
    I = concatenate((A.row,B.row+A.shape[0]))
    J = concatenate((A.col,B.col))
    V = concatenate((A.data,B.data))
    return coo_matrix((V,(I,J)),shape=(A.shape[0]+B.shape[0],A.shape[1])).tocsr()


##############################################################################################
#					JBS Utils			 	 	     #
##############################################################################################

def UnAmal(A,blocksize):
	#Input:	 A:		Amalmagated matrix, assumed to be in CSR format
	#	 blocksize:	Block size of unamalgamted matrix
	#
	#Output: A_UnAmal:	BSR matrix that is essentially a Kronecker product of 
	#			A and ones((blocksize,blocksize))
	data = ones( (A.indices.shape[0], blocksize, blocksize) )
	return bsr_matrix((data, A.indices, A.indptr), shape=(blocksize*A.shape[0], blocksize*A.shape[1]) )

def read_coord(filename):
	#Input:		filename:	File of x,y,z coordinates.  Formatting of file must be
	#				<Begin File>
	#				number_of_coordinates
	#				x_1	y_1	z_1
	#				x_2	y_2	z_2
	#				...	...	...
	#				x_n	y_n	z_n
	#				<End File>
	#
	#Output:	X,Y,Z:		(number_of_coordinate  x  1) vectors containing the coordinates

	fid = open(filename)

	N      = int(fid.readline()) 

	XYZ = fromfile(fid, sep=' ').reshape(-1,3)  # X,Y,Z in columns

	X  = ascontiguousarray(XYZ[:,0], dtype='float')
	Y  = ascontiguousarray(XYZ[:,1], dtype='float')
	Z  = ascontiguousarray(XYZ[:,2], dtype='float')

	return  X,Y,Z


def Coord2RBM(numNodes, numPDEs, x, y, z):
	#Input:		numNodes:	Number of nodes
	#		numPDEs:	Number of dofs per node
	#		x,y,z:		Coordinate vectors
	#
	#Output:	rbm:		Matrix of size (numNodes*numPDEs) x (1 | 6) containing the 6 rigid body modes

	#check inputs
	if(numPDEs == 1):
		numcols = 1
	elif( (numPDEs == 3) or (numPDEs == 6) ):
		numcols = 6
	else:
		raise ValueError("Coord2RBM(...) only supports 1, 3 or 6 PDEs per spatial location, i.e. numPDEs = [1 | 3 | 6].  You've entered " \
				+ str(numPDEs) + "." )

	if( (max(x.shape) != numNodes) or (max(y.shape) != numNodes) or (max(z.shape) != numNodes) ):
		raise ValueError("Coord2RBM(...) requires coordinate vectors of equal length.  Length must be numNodes = " + str(numNodes)) 

	#if( (min(x.shape) != 1) or (min(y.shape) != 1) or (min(z.shape) != 1) ):
	#	raise ValueError("Coord2RBM(...) requires coordinate vectors that are (numNodes x 1) or (1 x numNodes).") 


	#preallocate rbm
	rbm = mat(zeros((numNodes*numPDEs, numcols)))
	
	for node in range(numNodes):
		dof = node*numPDEs

		if(numPDEs == 1):
			rbm[node] = 1.0 
	            
		if(numPDEs == 6): 
			for ii in range(3,6):		#lower half = [ 0 I ]
				for jj in range(0,6):
					if(ii == jj):
						rbm[dof+ii, jj] = 1.0 
					else: 
						rbm[dof+ii, jj] = 0.0

		if((numPDEs == 3) or (numPDEs == 6) ): 
			for ii in range(0,3):		#upper left = [ I ]
				for jj in range(0,3):
					if(ii == jj):
						rbm[dof+ii, jj] = 1.0 
					else: 
						rbm[dof+ii, jj] = 0.0

			for ii in range(0,3):		#upper right = [ Q ]
				for jj in range(3,6):
					if( ii == (jj-3) ):
						rbm[dof+ii, jj] = 0.0
					else:
						if( (ii+jj) == 4):
							rbm[dof+ii, jj] = z[node]
						elif( (ii+jj) == 5 ): 
							rbm[dof+ii, jj] = y[node]
						elif( (ii+jj) == 6 ): 
							rbm[dof+ii, jj] = x[node]
		             			else:
							rbm[dof+ii, jj] = 0.0
			ii = 0 
			jj = 5 
			rbm[dof+ii, jj] *= -1.0
	
			ii = 1 
			jj = 3 
			rbm[dof+ii, jj] *= -1.0
	
			ii = 2 
			jj = 4 
			rbm[dof+ii, jj] *= -1.0
	
	return rbm


############################################################################################
#			JBS --- Define BSR helper functions				   #
############################################################################################

def BSR_Get_Row(A, i):
#	Input:	A:	Matrix assumed to be in BSR format
#		i:	row number
#
#	Output:	z:	 Actual nonzero values for row i
#		colindx: Array of column indices for the nonzeros of row i
	
	blocksize = A.blocksize[0]
	BlockIndx = i/blocksize
	rowstart = A.indptr[BlockIndx]
	rowend = A.indptr[BlockIndx+1]
	localRowIndx = i%blocksize

	#Get z
	indys = A.data[rowstart:rowend, localRowIndx, :].nonzero()
	z = A.data[rowstart:rowend, localRowIndx, :][indys[0], indys[1]]


	colindx = zeros((1, z.__len__()), dtype=int32)
	counter = 0

	for j in range(rowstart, rowend):
		coloffset = blocksize*A.indices[j]
		indys = A.data[j,localRowIndx,:].nonzero()[0]
		increment = indys.shape[0]
		colindx[0,counter:(counter+increment)] = coloffset + indys
		counter += increment		

	return mat(z).T, colindx[0,:]

def BSR_Row_WriteScalar(A, i, x):
#	Input:	A:	Matrix assumed to be in BSR format
#		i:	row number
#		x:	scalar to overwrite nonzeros of row i in A
#
#	Output:	A:	x is a scalar and all nonzeros in row i of A 
#			have been overwritten with x.  If x is a vector,
#			the first length(x) nonzeros in row i of A have been
#			overwritten with entries from x
	
	blocksize = A.blocksize[0]
	BlockIndx = i/blocksize
	rowstart = A.indptr[BlockIndx]
	rowend = A.indptr[BlockIndx+1]
	localRowIndx = i%blocksize

	#for j in range(rowstart, rowend):
	#	indys = A.data[j,localRowIndx,:].nonzero()[0]
	#	increment = indys.shape[0]
	#	A.data[j,localRowIndx,indys] = x
	
	indys = A.data[rowstart:rowend, localRowIndx, :].nonzero()
	A.data[rowstart:rowend, localRowIndx, :][indys[0], indys[1]] = x


def BSR_Row_WriteVect(A, i, x):
#	Input:	A:	Matrix assumed to be in BSR format
#		i:	row number
#		x:	Array of values to overwrite nonzeros in row i of A
#
#	Output:	A:	The nonzeros in row i of A have been
#			overwritten with entries from x.  x must be same
#			length as nonzeros of row i.  This is guaranteed
#			when this routine is used with vectors derived form
#			Get_BSR_Row
	
	blocksize = A.blocksize[0]
	BlockIndx = i/blocksize
	rowstart = A.indptr[BlockIndx]
	rowend = A.indptr[BlockIndx+1]
	localRowIndx = i%blocksize
	
	#Numpy occasionally seems to be written by idiot small children.  This line fixes one
	#	of the idiotic things about the array/matrix setup.  Sometimes I really wish
	#	for the Matlab matrix "slicing" interface rather than this.
	x = x.__array__().reshape( (max(x.shape),) )

	#counter = 0
	#for j in range(rowstart, rowend):
	#	indys = A.data[j,localRowIndx,:].nonzero()[0]
	#	increment = min(indys.shape[0], blocksize)
	#	A.data[j,localRowIndx,indys] = x[counter:(counter+increment), 0]
	#	counter += increment

	indys = A.data[rowstart:rowend, localRowIndx, :].nonzero()
	A.data[rowstart:rowend, localRowIndx, :][indys[0], indys[1]] = x

def BSR_Get_Colindices(A):
#	Input:	A:	BSR matrix from whom you want absolute column indices
#
#	Output: colindices:	List of arrays that is (num_nodes=A.rows/A.blocksize[0]) long.
#				Each entry holds that node's column indices, assuming the block 
#				is perfectly dense
	
	RowsPerBlock = A.blocksize[0]
	ColsPerBlock = A.blocksize[1]
	nRows = A.shape[0]
	colindices = []

	i = 0
	while(i < nRows):
		BlockIndx = i/RowsPerBlock
		rowstart = A.indptr[BlockIndx]
		rowend = A.indptr[BlockIndx+1]
		length = rowend - rowstart
		indys = zeros(length*ColsPerBlock, dtype=int32)

		counter = 0
		for j in range(rowstart, rowend):
			coloffset = ColsPerBlock*A.indices[j]
			indys[counter:(counter+ColsPerBlock)] = range(coloffset, coloffset+ColsPerBlock)
			counter = counter + ColsPerBlock
	
		colindices.append(indys)
	
		i += RowsPerBlock

	return colindices

###################################################################################################







































