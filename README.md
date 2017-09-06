[![Build Status](https://travis-ci.org/konovod/linalg.svg?branch=master)](https://travis-ci.org/konovod/linalg)

# linalg
Linear algebra library in Crystal, uses LAPACKE.
- direct access to LAPACK methods
- convenient Matrix(T) class, supports T=Float32, Float64 and Complex.
- high-level interface similar to scipy.linalg or MATLAB.

Killing SciPy, one module at a time.

## Installation

1. Install LAPACKE (and dependencies - LAPACK and BLAS). `sudo apt install libopenblas-base liblapacke` for Ubuntu, `sudo pacman -S lapacke` (for better performance use `openblas-lapack` package from AUR) for Arch.


2. (for Ubuntu) it seems package doesn't create symlink, so use
- `sudo ln -s /usr/lib/liblapacke.so.3 /usr/lib/liblapacke.so`
- `sudo ln -s /usr/lib/openblas-base/libblas.so.3 /usr/lib/libcblas.so`

Add this to your application's `shard.yml`:

```yaml
dependencies:
  linalg:
    github: konovod/linalg
```

## Usage

```crystal
require "linalg"
```
Basic type aliases are
- Mat = Matrix(Float64)
- Mat32 = Matrix(Float32)
- MatComplex = Matrix(Complex)

Complex consisting of two Float32 isn't supported for now (it is easy, but I'm not sure if it's useful).

Types with prefix G (GMat, GMat32, GMatComplex) are for actually allocated matrices,
others are automatically converted to them when needed.

```crystal
#suggested to don't prefix LA:: everywhere
include LA

# create matrix from array of arrays (or tuple... everything Indexable)
m = GMat[
  [1, 2, 3],
  [4, 5, 6],
  [7, 8, 9],
  [10, 11, 12],
]

# or using block
m = GMat32.new(3, 4) { |i, j| i*3 + j + 1 }
# or using one of other ways, check "spec" directory

# there are also matrices with special values, memory for elements isn't allocated until they are changed.
# NOTE currently virtual matrices support is incomplete, so memory is allocated always
a = Mat.identity(3)
puts a # =>
# [1.0, 0.0, 0.0]
# [0.0, 1.0, 0.0]
# [0.0, 0.0, 1.0]

# do basic arithmetics
pp 2 * a - Mat.diag([2, 2, 2]) == Mat.zeros(3, 3) # => true

# basic algebra
a = Mat.rand(5, 5) + 2 * Mat.identity(5)
pp (a.inv * a - Mat.identity(5)).norm < 1e-6

b = Mat.rand(5, 1)
x = LA.solve(a, b) # or a.solve(b)
pp (a*x - b).norm < 1e-6

m = GMat[[-2, 4, 1], [2, -4, 1], [1, 1, 1]]
pp m.eigvals # => [-6.0, -1.0, 2.0]

# extract submatrices (memory isn't copied as they reference to basic matrix)
m = GMat[
  [1, 2, 3],
  [4, 5, 6],
  [7, 8, 9],
]
pp m.columns[2] # LA::SubMatrix(Float64) (3x1, None):
# [3.0]
# [6.0]
# [9.0]

x = m[1..1, 1..2]
pp x        # => [5.0, 6.0]
x[0, 0] = 0 # m[1,1] is now 0 (questionable feature? maybe should be ##[]! for modifiable submatrices and ##[] for CoW?)
y = x.clone # now y is a separate matrix
y[0, 0] = 1 # m[1,1] is still 0
pp m[1, 1]

```
other present features:

- svd (`Mat##svd` or `Mat##svdvals` for just values)
- lu decomposition (`Mat##lu`)
```crystal
# to just get P L U matrices
p, l, u = a.lu

# to get them in compact form and use for solving linear equations:
a = GMat32[
  [2, 4],
  [2, 8]
  ]

lu = a.lu_factor # lu is LUMatrix(T) - immutable object that can return it's content and solve systems
puts lu.solve(GMat32[[2], [4]])
```
- matrix rank determination (using SVD or QRP)
- linear least squares problem (`LA.solvels` to just get decision or `LA.lstsq` to also get rank and singular values (TODO - and residues))
- cholesky decomposition (`##cholesky`, `##cholesky!`, `##cho_solve`)
- `hessenberg` form
- `qr`, `rq`, `lq`, `ql` decompositions
- `schur` and `qz` (generalized schur) decomposition
- generalized eigenproblem (`eigs(a, b, ...)`)
- creating special matrices like `pascal` or `toeplitz` (check scipy.md for a full list)
- matrix exponent and trigonomertic functions
- matrix exponentiation (to integer powers only atm (TODO - fractional))


There is also concept of `Mat##flags` that represent properties of matrix (symmetric, positive definite etc), they are used to automatically select faster algorithms from LAPACK. Flags are partially enforced by runtime checks, with the possibility of user override. For example, if we say that `a.assume!(MatrixFlags::Symmetric)` then `a.transpose` or `a + Mat.diag(*a.size)` will also have this flag, so the LAPACK routines for symmetrical matrices will be used. In fact, `a.transpose` will return matrix clone as for symmetric matrices A=A'.

Supported flags:
```crystal
enum MatrixFlags
  Symmetric
  Hermitian
  PositiveDefinite
  Orthogonal
  UpperTriangular
  LowerTriangular
  Triangular      = UpperTriangular | LowerTriangular
```
NOTE for complex matrices `Orthogonal` flag means `Unitary`.

Main functions for flags are:
```crystal
  a.assume!(flag) # sets matrix flag without check, can lead to incorrect results if matrix do not have corresponding property.
  a.detect(flag) # checks if matrix has property, if yes sets the flag. Returns true if check positive
  a.detect # detect all possible flags
  a.flags # returns matrix flags
```
Most operations - matrix addition, multiplication, inversion, transposition and decompositions correctly update flags, but any direct access like `a[i,j] = 0` or `a.map!{|v| v+1}` resets flags to `None`, so use `a.detect` after them if you need to preserve flags (or `a.assume!(f)` if detection is too slow).

## Development

### Roadmap:

##### Important

- [x] saving/loading from files
- [ ] Virtual(lazy) matrices, ways to evade allocations during calculations
- [x] Matrix exponent and trigonometric
- [ ] other matrix functions
- [ ] Banded matrices
- [ ] Column-major storage (optional?)
- [ ] Other missing features from LAPACK (mostly selectable and orderable eigenvalues)
- [ ] Sparse matrices (perhaps out of scope/deserves separate shard)
- [ ] Other missing features from scipy.linalg (pseudoinverse, lyapunov/ricatti/sylvester equations, other things i don't know algorithms for)

##### Not so important

- [x] saving/loading to matlab-like string
- [ ] better pretty-print, with alignment and various precision
- [x] use blas for multiplication
- [ ] more flags support (inversion of diagonal matrix and other trivial cases)

## Contributing

1. Fork it ( https://github.com/konovod/linalg/fork )
2. Create your feature branch (git checkout -b my-new-feature)
3. Commit your changes (git commit -am 'Add some feature')
4. Push to the branch (git push origin my-new-feature)
5. Create a new Pull Request
