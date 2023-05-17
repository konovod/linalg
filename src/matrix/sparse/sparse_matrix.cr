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

    # returns element-wise sum
    def +(m : Sparse::Matrix)
      result = clone.add!(T.new(1), m)
      result.flags = self.flags.sum(m.flags)
      result
    end

    def -(m : Sparse::Matrix)
      result = clone.add!(-T.new(1), m)
      result.flags = self.flags.sum(m.flags)
      result
    end

    def +(m : LA::Matrix)
      m.clone.add! self
    end

    def -(m : LA::Matrix)
      (-m).add! self
    end

    def add!(k : Number | Complex, m : LA::Matrix)
      raise ArgumentError.new "can't `add!` dense matrix to sparse"
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
  end
end
