## Scipy.linalg functionality:

### Basics
  
- [x] inv(a[, overwrite_a, check_finite])	Compute the inverse of a matrix.
- [x] solve(a, b[, sym_pos, lower, overwrite_a, ...])	Solves the linear equation set a * x = b for the unknown x for square a matrix.
- [x] solve_banded(l_and_u, ab, b[, overwrite_ab, ...])	Solve the equation a x = b for x, assuming a is banded matrix.
- [x] solveh_banded(ab, b[, overwrite_ab, ...])	Solve equation a x = b.
- [ ] solve_circulant(c, b[, singular, tol, ...])	Solve C x = b for x, where C is a circulant matrix.
- [x] solve_triangular(a, b[, trans, lower, ...])	Solve the equation a x = b for x, assuming a is a triangular matrix.
- [ ] solve_toeplitz(c_or_cr, b[, check_finite])	Solve a Toeplitz system using Levinson Recursion
- [ ] matmul_toeplitz(c_or_cr, x[, check_finite, ...])  Efficient Toeplitz Matrix-Matrix Multiplication using FFT
- [x] det(a[, overwrite_a, check_finite])	Compute the determinant of a matrix
- [x] norm(a[, ord, axis, keepdims])	Matrix or vector norm.
- [x] lstsq(a, b[, cond, overwrite_a, ...])	Compute least-squares solution to equation Ax = b.
- [x] pinv(a[, cond, rcond, return_rank, check_finite])	Compute the (Moore-Penrose) pseudo-inverse of a matrix.
- [x] pinvh(a[, cond, rcond, lower, return_rank, ...])	Compute the (Moore-Penrose) pseudo-inverse of a Hermitian matrix.
- [x] kron(a, b)	Kronecker product.
- [ ] khatri_rao(a, b)  Khatri-rao product
- [x] tril(m[, k])	Make a copy of a matrix with elements above the k-th diagonal zeroed.
- [x] triu(m[, k])	Make a copy of a matrix with elements below the k-th diagonal zeroed.
- [ ] orthogonal_procrustes(A, B[, check_finite])	Compute the matrix solution of the orthogonal Procrustes problem.
- [x] matrix_balance(A[, permute, scale, ...])	Compute a diagonal similarity transformation for row/column balancing.
- [ ] subspace_angles(A, B) Compute the subspace angles between two matrices.
- [ ] bandwidth(a)  Return the lower and upper bandwidth of a 2D numeric array.
- [x] issymmetric(a[, atol, rtol])  Check if a square 2D array is symmetric.
- [x] ishermitian(a[, atol, rtol])  Check if a square 2D array is Hermitian.
- [x] LinAlgError	Generic Python-exception-derived object raised by linalg functions.
- [ ] LinAlgWarning The warning emitted when a linear algebra related operation is close to fail conditions of the algorithm or loss of accuracy is expected.

### Eigenvalue Problems

- [x] eig(a[, b, left, right, overwrite_a, ...])	Solve an ordinary or generalized eigenvalue problem of a square matrix.
- [x] eigvals(a[, b, overwrite_a, check_finite, ...])	Compute eigenvalues from an ordinary or generalized eigenvalue problem.
- [x] eigh(a[, b, lower, eigvals_only, ...])	Solve an ordinary or generalized eigenvalue problem for a complex Hermitian or real symmetric matrix.
- [x] eigvalsh(a[, b, lower, overwrite_a, ...])	Solve an ordinary or generalized eigenvalue problem for a complex Hermitian or real symmetric matrix.
- [ ] eig_banded(a_band[, lower, eigvals_only, ...])	Solve real symmetric or complex hermitian band matrix eigenvalue problem.
- [ ] eigvals_banded(a_band[, lower, ...])	Solve real symmetric or complex hermitian band matrix eigenvalue problem.
- [ ] eigh_tridiagonal(d, e[, eigvals_only, ...]) Solve eigenvalue problem for a real symmetric tridiagonal matrix.
- [ ] eigvalsh_tridiagonal(d, e[, select, ...]) Solve eigenvalue problem for a real symmetric tridiagonal matrix.

### Decompositions

- [x] lu(a[, permute_l, overwrite_a, check_finite])	Compute pivoted LU decomposition of a matrix.
- [x] lu_factor(a[, overwrite_a, check_finite])	Compute pivoted LU decomposition of a matrix.
- [x] lu_solve(lu_and_piv, b[, trans, ...])	Solve an equation system, a x = b, given the LU factorization of a
- [x] svd(a[, full_matrices, compute_uv, ...])	Singular Value Decomposition.
- [x] svdvals(a[, overwrite_a, check_finite])	Compute singular values of a matrix.
- [x] diagsvd(s, M, N)	Construct the sigma matrix in SVD from singular values and size M, N.
- [ ] orth(A)	Construct an orthonormal basis for the range of A using SVD
- [ ] null_space(A[, rcond])  Construct an orthonormal basis for the null space of A using SVD
- [ ] ldl(A[, lower, hermitian, overwrite_a, ...])  Computes the LDLt or Bunch-Kaufman factorization of a symmetric/ hermitian matrix.
- [x] cholesky(a[, lower, overwrite_a, check_finite])	Compute the Cholesky decomposition of a matrix.
- [ ] cholesky_banded(ab[, overwrite_ab, lower, ...])	Cholesky decompose a banded Hermitian positive-definite matrix
- [x] cho_factor(a[, lower, overwrite_a, check_finite])	Compute the Cholesky decomposition of a matrix, to use in cho_solve
- [x] cho_solve(c_and_lower, b[, overwrite_b, ...])	Solve the linear equations A x = b, given the Cholesky factorization of A.
- [ ] cho_solve_banded(cb_and_lower, b[, ...])	Solve the linear equations A x = b, given the Cholesky factorization of A.
- [ ] polar(a[, side])	Compute the polar decomposition.
- [x] qr(a[, overwrite_a, lwork, mode, pivoting, ...])	Compute QR decomposition of a matrix.
- [ ] qr_multiply(a, c[, mode, pivoting, ...])	Calculate the QR decomposition and multiply Q with a matrix.
- [ ] qr_update(Q, R, u, v[, overwrite_qruv, ...])	Rank-k QR update
- [ ] qr_delete(Q, R, k, int p=1[, which, ...])	QR downdate on row or column deletions
- [ ] qr_insert(Q, R, u, k[, which, rcond, ...])	QR update on row or column insertions
- [x] rq(a[, overwrite_a, lwork, mode, check_finite])	Compute RQ decomposition of a matrix.
- [x] qz(A, B[, output, lwork, sort, overwrite_a, ...])	QZ decomposition for generalized eigenvalues of a pair of matrices.
- [ ] ordqz(A, B[, sort, output, overwrite_a, ...])	QZ decomposition for a pair of matrices with reordering.
- [x] schur(a[, output, lwork, overwrite_a, sort, ...])	Compute Schur decomposition of a matrix.
- [ ] rsf2csf(T, Z[, check_finite])	Convert real Schur form to complex Schur form.
- [x] hessenberg(a[, calc_q, overwrite_a, ...])	Compute Hessenberg form of a matrix.
- [ ] cdf2rdf(w, v) Converts complex eigenvalues w and eigenvectors v to real eigenvalues in a block diagonal form wr and the associated real eigenvectors vr, such that.
- [ ] cossin(X[, p, q, separate, swap_sign, ...]) Compute the cosine-sine (CS) decomposition of an orthogonal/unitary matrix.
- [ ] scipy.linalg.interpolative – Interpolative matrix decompositions

### Matrix Functions

- [x] expm(A[, q])	Compute the matrix exponential using Pade approximation.
- [ ] logm(A[, disp])	Compute matrix logarithm.
- [x] cosm(A)	Compute the matrix cosine.
- [x] sinm(A)	Compute the matrix sine.
- [x] tanm(A)	Compute the matrix tangent.
- [x] coshm(A)	Compute the hyperbolic matrix cosine.
- [x] sinhm(A)	Compute the hyperbolic matrix sine.
- [x] tanhm(A)	Compute the hyperbolic matrix tangent.
- [ ] signm(A[, disp])	Matrix sign function.
- [ ] sqrtm(A[, disp, blocksize])	Matrix square root.
- [ ] funm(A, func[, disp])	Evaluate a matrix function specified by a callable.
- [ ] expm_frechet(A, E[, method, compute_expm, ...])	Frechet derivative of the matrix exponential of A in the direction E.
- [ ] expm_cond(A[, check_finite])	Relative condition number of the matrix exponential in the Frobenius norm.
- [ ] fractional_matrix_power(A, t)	Compute the fractional power of a matrix.

### Matrix Equation Solvers

- [ ] solve_sylvester(a, b, q)	Computes a solution (X) to the Sylvester equation AX+XB=QAX+XB=Q.
- [ ] solve_continuous_are(a, b, q, r[, e, s, ...])	Solves the continuous-time algebraic Riccati equation (CARE).
- [ ] solve_discrete_are(a, b, q, r[, e, s, balanced])	Solves the discrete-time algebraic Riccati equation (DARE).
- [ ] solve_continuous_lyapunov(a, q)	Solves the continuous Lyapunov equation AX+XAH=QAX+XAH=Q.
- [ ] solve_discrete_lyapunov(a, q[, method])	Solves the discrete Lyapunov equation AXAH−X+Q=0AXAH−X+Q=0.

### Sketches and Random Projections

- [ ] clarkson_woodruff_transform(input_matrix, ...)  Applies a Clarkson-Woodruff Transform/sketch to the input matrix.

### Special Matrices

- [x] block_diag(arrs)	Create a block diagonal matrix from provided arrays.
- [x] circulant(c)	Construct a circulant matrix.
- [x] companion(a)	Create a companion matrix.
- [x] dft(n[, scale])	Discrete Fourier transform matrix.
- [x] fiedler(a)  Returns a symmetric Fiedler matrix
- [ ] fiedler_companion(a)  Returns a Fiedler companion matrix
- [x] hadamard(n[, dtype])	Construct a Hadamard matrix.
- [x] hankel(c[, r])	Construct a Hankel matrix.
- [x] helmert(n[, full])	Create a Helmert matrix of order n.
- [x] hilbert(n)	Create a Hilbert matrix of order n.
- [ ] invhilbert(n[, exact])	Compute the inverse of the Hilbert matrix of order n.
- [x] leslie(f, s)	Create a Leslie matrix.
- [x] pascal(n[, kind, exact])	Returns the n x n Pascal matrix.
- [x] invpascal(n[, kind, exact])	Returns the inverse of the n x n Pascal matrix.
- [x] toeplitz(c[, r])	Construct a Toeplitz matrix.
- [x] tri(N[, M, k, dtype])	Construct (N, M) matrix filled with ones at and below the k-th diagonal.
