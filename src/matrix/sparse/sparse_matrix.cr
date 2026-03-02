# Abstract base class for sparse matrix formats.
#
# Provides common functionality for sparse matrix implementations (CSR, CSC, COO).
# Sparse matrices store only non-zero elements, saving memory for matrices with sparsity.

module LA::Sparse
  abstract class Matrix(T) < LA::Matrix(T)
    getter nrows : Int32 = 0
    getter ncolumns : Int32 = 0
    property flags : MatrixFlags = MatrixFlags::None

    # Number of non-zero elements in the matrix.
    abstract def nonzeros : Int32

    # Compares this sparse matrix with another for equality.
    #
    # Arguments:
    #   - other (Sparse::Matrix(T)) : Matrix to compare with.
    #
    # Returns:
    #   - Bool : `true` if matrices are equal.
    def ==(other : Sparse::Matrix(T))
      return false unless nrows == other.nrows && ncolumns == other.ncolumns
      if other.nonzeros + self.nonzeros < nrows * ncolumns // 4
        a, b = self, other
        a, b = b, a if a.nonzeros < b.nonzeros
        a.each_with_index(all: false) do |v, i, j|
          return false if b.unsafe_fetch(i, j) != v
        end
        b.each_with_index(all: false) do |v, i, j|
          return false if a.unsafe_fetch(i, j) != v
        end
        true
      else
        super(other)
      end
    end

    # Converts the sparse matrix to a dense GeneralMatrix.
    #
    # Returns:
    #   - GeneralMatrix(T) : Dense matrix representation.
    def to_general
      result = GeneralMatrix(T).new(nrows, ncolumns)
      self.each_with_index(all: false) do |v, i, j|
        result.unsafe_set(i, j, v)
      end
      result.flags = flags
      result
    end

    def to_dense
      to_general
    end

    def dup
      self.class.new(self)
    end

    def clone
      dup
    end

    def add(m : DenseMatrix, *, alpha = 1, beta = 1)
      m.add(self, alpha: beta, beta: alpha)
    end

    # In-place: keeps only elements on or above the k-th diagonal.
    #
    # Arguments:
    #   - k (Int32) : Diagonal offset. Default: 0.
    #
    # Returns:
    #   - self : The modified matrix.
    def triu!(k = 0)
      flags = self.flags
      select_index! { |i, j| i <= j - k }
      self.flags = self.flags.triu(k >= 0, square?)
      self
    end

    # In-place: keeps only elements on or below the k-th diagonal.
    #
    # Arguments:
    #   - k (Int32) : Diagonal offset. Default: 0.
    #
    # Returns:
    #   - self : The modified matrix.
    def tril!(k = 0)
      flags = self.flags
      select_index! { |i, j| i >= j - k }
      self.flags = self.flags.tril(k <= 0, square?)
      self
    end

    # Returns a new matrix with elements below the k-th diagonal zeroed.
    #
    # Arguments:
    #   - k (Int32) : Diagonal offset. Default: 0.
    #
    # Returns:
    #   - Sparse::Matrix(T) : New matrix with elements below diagonal zeroed.
    def tril(k = 0)
      result = select_with_index do |v, i, j|
        i >= j - k
      end
      result.flags = self.flags.tril(k <= 0, square?)
      result
    end

    # Returns a new matrix with elements above the k-th diagonal zeroed.
    #
    # Arguments:
    #   - k (Int32) : Diagonal offset. Default: 0.
    #
    # Returns:
    #   - Sparse::Matrix(T) : New matrix with elements above diagonal zeroed.
    def triu(k = 0)
      result = select_with_index do |v, i, j|
        i <= j - k
      end
      result.flags = self.flags.triu(k >= 0, square?)
      result
    end

    # Returns a new matrix keeping only elements where the block returns true.
    #
    # Arguments:
    #   - block : A block that receives the value and returns true to keep.
    #
    # Returns:
    #   - Sparse::Matrix(T) : New matrix with filtered elements.
    def select(& : T -> Bool)
      select_with_index { |v, i, j| yield(v) }
    end

    # Returns a new matrix keeping only elements at positions where the block returns true.
    #
    # Arguments:
    #   - block : A block that receives (i, j) and returns true to keep.
    #
    # Returns:
    #   - Sparse::Matrix(T) : New matrix with filtered elements.
    def select_index(& : (Int32, Int32) -> Bool)
      select_with_index { |v, i, j| yield(i, j) }
    end
  end
end
