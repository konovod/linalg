require "complex"
require "./matrix"
require "./general_matrix"

module Linalg
  module Matrix(T)
    def self.tri(nrows, ncolumns, k = 0)
      GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        i >= j - k ? 1 : 0
      end
    end

    def self.block_diag(*args)
      nrows = args.sum &.nrows
      ncolumns = args.sum &.ncolumns
      GeneralMatrix(T).new(nrows, ncolumns).tap do |result|
        row = 0
        column = 0
        args.each do |arg|
          result[row...row + arg.nrows, column...column + arg.ncolumns] = arg
          row += arg.nrows
          column += arg.ncolumns
        end
      end
    end

    def self.toeplitz(column : Indexable | Matrix, row : Indexable | Matrix | Nil = nil)
      row = row.to_a if row.is_a? Matrix
      column = column.to_a if column.is_a? Matrix
      if row
        GeneralMatrix(T).new(column.size, row.size) do |i, j|
          k = i - j
          if k >= 0
            column[k]
          else
            row[-k]
          end
        end
      else
        GeneralMatrix(T).new(column.size, column.size) do |i, j|
          k = i - j
          if k >= 0
            column[k]
          else
            column[-k].conj
          end
        end
      end
    end

    def self.circulant(c)
      GeneralMatrix(T).new(c.size, c.size) do |i, j|
        k = i - j
        c[(k + c.size) % c.size]
      end
    end

    def self.leslie(f, s)
      GeneralMatrix(T).new(s.size + 1, f.size).tap do |matrix|
        f.each_with_index { |fi, i| matrix.unsafe_set 0, i, T.new(fi) }
        s.each_with_index { |si, i| matrix.unsafe_set i + 1, i, T.new(si) }
      end
    end

    def self.companion(a)
      k = -1.0/a[0]
      GeneralMatrix(T).new(a.size - 1, a.size - 1).tap do |matrix|
        (a.size - 1).times { |i| matrix.unsafe_set 0, i, T.new(a[i + 1]*k) }
        (a.size - 2).times { |i| matrix.unsafe_set i + 1, i, T.new(1) }
      end
    end

    # TODO - faster implementation
    def self.hadamard(n)
      raise ArgumentError.new("size must be positive") unless n > 0
      raise ArgumentError.new("size must be power of two") unless n.popcount == 1
      return GeneralMatrix(T).new([[1]]) if n == 1
      return GeneralMatrix(T).new([[1, 1], [1, -1]]) if n == 2
      return hadamard(n/2).kron(hadamard(2))
    end

    def self.hankel(column : Indexable | Matrix, row : Indexable | Matrix | Nil = nil)
      row = row.to_a if row.is_a? Matrix
      column = column.to_a if column.is_a? Matrix
      if row
        GeneralMatrix(T).new(column.size, row.size) do |i, j|
          k = i + j
          if k < column.size
            column[k]
          else
            row[k - column.size + 1]
          end
        end
      else
        GeneralMatrix(T).new(column.size, column.size) do |i, j|
          k = i + j
          if k < column.size
            column[k]
          else
            0
          end
        end
      end
    end

    def self.helmert(n, full = false)
      if full
        result = GeneralMatrix(T).new(n, n)
      else
        result = GeneralMatrix(T).new(n - 1, n)
      end
      # first row
      if full
        result[0, 0...n] = T.new(Math.sqrt(1.0/n))
        rowdelta = 1
      else
        rowdelta = 0
      end
      # rest
      (n - 1).times do |i|
        x = i + 1
        v = T.new(Math.sqrt(1.0/(x + x*x)))
        result.unsafe_set i + rowdelta, i + 1, -v*x
        result[i + rowdelta, 0..i] = v
      end
      result
    end

    def self.hilbert(n)
      GeneralMatrix(T).new(n, n) do |i, j|
        T.new(1.0) / (i + j + 1)
      end
    end
  end

  enum DFTScale
    None
    N
    SqrtN
  end

  module Matrix(T)
    def self.dft(n, scale : DFTScale = DFTScale::None)
      {% raise "DFT matrix must be Complex" unless T == Complex %}
      j = Complex.new(0, 1)
      w = (-2*Math::PI*j / n).exp
      result = Matrix(T).ones(n, n).clone
      result.each_index do |i, j|
        next if i == 0 || j == 0
        if j == 1
          result.unsafe_set(i, j, w*result.unsafe_at(i - 1, j))
        else
          result.unsafe_set(i, j, result.unsafe_at(i, 1)*result.unsafe_at(i, j - 1))
        end
      end
      unless scale.none?
        case scale
        when .sqrt_n?
          scale = 1.0 / Math.sqrt(n)
        else
          scale = 1.0 / n
        end
        result.map! { |v| scale*v }
      end
      result
    end
  end
end

module Linalg
  enum PascalKind
    Upper
    Lower
    Symmetric
  end

  module Matrix(T)
    private def self.n_choose_k(n, k)
      (1..n).product(T.new(1.0)) / ((1..k).product(T.new(1.0)) * (1..n - k).product(T.new(1.0)))
    end

    def self.pascal(n, kind : PascalKind = PascalKind::Symmetric)
      case kind
      when .upper?
        zeros(n, n).tap do |m|
          m.each_index do |i, j|
            next if i > j
            if i == j || i == 0
              m.unsafe_set(i, j, T.new(1.0))
            else
              m.unsafe_set(i, j, m.unsafe_at(i - 1, j - 1) + m.unsafe_at(i, j - 1))
            end
          end
          m.assume!(MatrixFlags::UpperTriangular)
        end
      when .lower?
        zeros(n, n).tap do |m|
          m.each_index do |i, j|
            next if i < j
            if i == j || j == 0
              m.unsafe_set(i, j, T.new(1.0))
            else
              m.unsafe_set(i, j, m.unsafe_at(i - 1, j - 1) + m.unsafe_at(i - 1, j))
            end
          end
          m.assume!(MatrixFlags::LowerTriangular)
        end
      else
        zeros(n, n).tap do |m|
          m.each_index do |i, j|
            if i == 0 || j == 0
              m.unsafe_set(i, j, T.new(1.0))
            else
              m.unsafe_set(i, j, m.unsafe_at(i - 1, j) + m.unsafe_at(i, j - 1))
            end
          end
          m.assume!(MatrixFlags::Symmetric | MatrixFlags::Hermitian)
        end
      end
    end
  end
end
