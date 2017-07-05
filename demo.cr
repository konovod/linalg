require "./src/linalg"
include Linalg

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
a = Mat.identity(3) # =>
# [1 0 0]
# [0 1 0]
# [0 0 1]

# do basic arithmetics
2 * a - Mat.diag([2, 2, 2]) == Mat.zeros(3, 3) # => true

# basic algebra
a = Mat.rand(5, 5) + 2 * Mat.identity(5)
(a.inv * a - Mat.identity(5)).abs < 1e-6

b = Mat.rand(5, 1)
x = Linalg.solve(a, b) # or a.solve(b)
(a.inv*x - b).abs < 1e-6

m = GMat[[-2, 4, 1], [2, -4, 1], [1, 1, 1]]
m.eigvals # => [-6, -1, 2]

# extract submatrices (memory isn't copied as they reference to basic matrix)
m = GMat[
  [1, 2, 3],
  [4, 5, 6],
  [7, 8, 9],
])
r = m.columns[2] # =>
# [3]
# [6]
# [7]
x = m[1..1, 1..2] # =>
# [5 6]
x[0, 0] = 0 # m[1,1] is now 0 (questionable feature? maybe should be ##[]! for modifiable submatrices and ##[] for CoW?)
y = x.clone # now y is a separate matrix
y[0, 0] = 1 # m[1,1] is still 0
