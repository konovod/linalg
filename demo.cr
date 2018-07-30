require "./src/linalg"
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

# do basic arithmetics
a = Mat.eye(3)
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
x[0, 0] = 0 # m[1,1] is now 0 (questionable feature? maybe should be #[]! for modifiable submatrices and #[] for CoW?)
y = x.clone # now y is a separate matrix
y[0, 0] = 1 # m[1,1] is still 0
pp m[1, 1]
