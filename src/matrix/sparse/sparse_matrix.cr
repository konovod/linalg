# TODO - inline docs

module LA::Sparse
  abstract class Matrix(T) < LA::Matrix(T)
    getter nrows : Int32 = 0
    getter ncolumns : Int32 = 0
    property flags : MatrixFlags = MatrixFlags::None

    abstract def nonzeros : Int32

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

    def triu!(k = 0)
      flags = self.flags
      select_index! { |i, j| i <= j - k }
      self.flags = self.flags.triu(k >= 0, square?)
      self
    end

    def tril!(k = 0)
      flags = self.flags
      select_index! { |i, j| i >= j - k }
      self.flags = self.flags.tril(k <= 0, square?)
      self
    end

    def tril(k = 0)
      result = select_with_index do |v, i, j|
        i >= j - k
      end
      result.flags = self.flags.tril(k <= 0, square?)
      result
    end

    def triu(k = 0)
      result = select_with_index do |v, i, j|
        i <= j - k
      end
      result.flags = self.flags.triu(k >= 0, square?)
      result
    end

    def select(& : T -> Bool)
      select_with_index { |v, i, j| yield(v) }
    end

    def select_index(& : (Int32, Int32) -> Bool)
      select_with_index { |v, i, j| yield(i, j) }
    end

    def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
      result = T.additive_identity
      case kind
      in .frobenius?
        each { |v| result += v*v }
        Math.sqrt(result)
      in .one?
        sums = Slice(T).new(ncolumns, T.additive_identity)
        each_with_index do |v, row, col|
          sums[col] += v.abs
        end
        sums.max
      in .inf?
        sums = Slice(T).new(nrows, T.additive_identity)
        each_with_index do |v, row, col|
          sums[row] += v.abs
        end
        sums.max
      in .max_abs?
        each { |v| result = v.abs if result < v }
        result
      end
    end
  end
end
