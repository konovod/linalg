require "./matrix"
require "./submatrix"

module LA
  # banded matrix, heap-allocated
  class BandedMatrix(T) < Matrix(T)
    getter nrows : Int32
    getter ncolumns : Int32
    property flags = MatrixFlags::None
    getter upper_band : Int32
    getter lower_band : Int32
    getter raw_banded : Slice(T)

    @[AlwaysInline]
    private def band_len
      ncolumns
    end

    @[AlwaysInline]
    private def bands_size
      band_len*(@upper_band + @lower_band + 1)
    end

    private def ij2index(i, j)
      return nil unless {0, j - @upper_band}.max <= i <= {nrows - 1, j + @lower_band}.min
      ai = @upper_band + i - j
      ai*band_len + j
    end

    private def index2ij(index) : {Int32, Int32}?
      ai = index/band_len
      j = index % band_len
      i = ai + j - @upper_band
      return nil if i < 0 || i >= nrows
      {i, j}
    end

    def initialize(@nrows, @ncolumns, @upper_band, @lower_band = upper_band, @flags = MatrixFlags::None)
      check_type
      @raw_banded = Slice(T).new(bands_size, T.new(0))
    end

    def initialize(@nrows, @ncolumns, @upper_band, @lower_band = upper_band, @flags = MatrixFlags::None, &block)
      check_type
      @raw_banded = Slice(T).new(bands_size) do |index|
        if ij = index2ij(index)
          i, j = ij
          T.new(yield(i, j))
        else
          T.new(0)
        end
      end
    end

    def self.new(nrows, ncolumns, upper_band, values : Indexable)
      new(nrows, ncolumns, upper_band, upper_band, values)
    end

    def initialize(@nrows, @ncolumns, @upper_band, @lower_band, values : Indexable)
      check_type
      # TODO check values sizes
      @raw_banded = Slice(T).new(bands_size) do |index|
        band = index/band_len
        offset = band - @upper_band
        i = (index % band_len)
        if offset >= 0
          len = {nrows - offset, ncolumns}.min
        else
          len = {nrows, ncolumns + offset}.min
        end
        raise ArgumentError.new("#{band} item length is #{values[band].size}, expected #{len} ") if i == 0 && len != values[band].size
        i += offset if offset < 0
        if i >= 0 && i < len
          T.new(values[band][i])
        else
          T.new(0)
        end
      end
    end

    def self.new(matrix : BandedMatrix)
      new(matrix.nrows, matrix.ncolumns, matrix.upper_band, matrix.lower_band, matrix.flags).tap do |result|
        matrix.raw_banded.each_with_index do |v, i|
          result.raw_banded[i] = T.new(v)
        end
      end
    end

    def self.new(matrix : Matrix, tolerance = matrix.tolerance)
      upper_band = (1..matrix.ncolumns - 1).bsearch do |i|
        matrix.diag(i).all? { |v| v.abs <= tolerance }
      end
      upper_band = (upper_band || matrix.ncolumns) - 1

      lower_band = (1..matrix.nrows - 1).bsearch do |i|
        matrix.diag(-i).all? { |v| v.abs <= tolerance }
      end
      lower_band = (lower_band || matrix.nrows) - 1

      new(matrix.nrows, matrix.ncolumns, upper_band, lower_band, matrix.flags) do |i, j|
        matrix.unsafe_at(i, j)
      end
    end

    def unsafe_at(i, j)
      if index = ij2index(i, j)
        @raw_banded.unsafe_at(index)
      else
        T.new(0.0)
      end
    end

    def unsafe_set(i, j, value)
      clear_flags # TODO - not always?
      if index = ij2index(i, j)
        @raw_banded[index] = T.new(value)
      else
        # TODO - resize?
        raise IndexError.new
      end
    end

    def dup
      BandedMatrix(T).new(self)
    end

    def clone
      dup
    end

    def each_index(*, all = false, &block)
      if all
        super(all: false) { |i, j| yield(i, j) }
      else
        bands_size.times do |i|
          if rowcol = index2ij(i)
            yield(*rowcol)
          end
        end
      end
    end

    def map_with_index(&block)
      BandedMatrix(T).new(nrows, ncolumns, upper_band, lower_band) { |i, j| yield(unsafe_at(i, j), i, j) }
    end

    def to_unsafe
      {% if T == Complex %}
        @raw_banded.to_unsafe.as(LibCBLAS::ComplexDouble*)
      {% else %}
        @raw_banded.to_unsafe
      {% end %}
    end

    def ==(other : BandedMatrix(T))
      # optimization - don't check corners that are empty on both matrices
      return false unless nrows == other.nrows && ncolumns == other.ncolumns
      maxupper = {upper_band, other.upper_band}.max
      maxlower = {lower_band, other.lower_band}.max
      range = -maxupper..maxlower
      each_index(all: true) do |i, j|
        next unless range.includes? i - j
        return false if other.unsafe_at(i, j) != unsafe_at(i, j)
      end
      true
    end

    def transpose!
      return self if flags.symmetric?
      newraw = Slice(T).new(nrows*(@upper_band + @lower_band + 1), T.new(0))
      each_with_index do |v, i, j|
        # raise "#{i}, #{j}" unless {0, j - @lower_band}.max <= j <= {ncolumns - 1, i + @upper_band}.min
        ai = @lower_band + j - i
        newraw[ai*nrows + i] = v
      end
      @upper_band, @lower_band = @lower_band, @upper_band
      @nrows, @ncolumns = @ncolumns, @nrows
      @raw_banded = newraw
      self.flags = flags.transpose
      self
    end

    # def reshape!(anrows, ancolumns)
    # def tril!(k = 0)
    # def triu!(k = 0)
    # def to_real
    # def to_imag
    # def +(m : Matrix(T))
    # def -(m : Matrix(T))
    # def transpose
    # def conjtranspose
    # def kron(b : Matrix(T))
    # def tril(k = 0)
    # def triu(k = 0)
    # def self.rand(nrows, ncolumns, rng = Random::DEFAULT)
    # def self.diag(nrows, ncolumns, values)
    # def self.diag(nrows, ncolumns, &block)
    # def self.identity(n)
    # def cat(other : Matrix(T), dimension)
  end

  alias BMat = BandedMatrix(Float64)
  alias BMat32 = BandedMatrix(Float32)
  alias BMatComplex = BandedMatrix(Complex)
end
