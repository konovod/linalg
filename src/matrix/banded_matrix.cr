require "./matrix"
require "./submatrix"

module LA
  # banded matrix, heap-allocated
  class BandedMatrix(T) < Matrix(T)
    getter nrows : Int32
    getter ncolumns : Int32
    property flags : MatrixFlags = MatrixFlags::None
    getter upper_band : Int32
    getter lower_band : Int32
    protected getter raw_banded : Slice(T)

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
      ai = index//band_len
      j = index % band_len
      i = ai + j - @upper_band
      return nil if i < 0 || i >= nrows
      {i, j}
    end

    def initialize(@nrows, @ncolumns, @upper_band, @lower_band = upper_band, @flags = MatrixFlags::None)
      check_type
      @raw_banded = Slice(T).new(bands_size, T.new(0))
      resized_flags
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
      resized_flags
    end

    def self.new(nrows, ncolumns, upper_band, values : Indexable)
      new(nrows, ncolumns, upper_band, upper_band, values)
    end

    def initialize(@nrows, @ncolumns, @upper_band, @lower_band, values : Indexable)
      check_type
      # TODO check values sizes
      @raw_banded = Slice(T).new(bands_size) do |index|
        band = index//band_len
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
      resized_flags
    end

    def self.new(matrix : BandedMatrix)
      new(matrix.nrows, matrix.ncolumns, matrix.upper_band, matrix.lower_band, matrix.flags).tap do |result|
        matrix.raw_banded.each_with_index do |v, i|
          result.raw_banded[i] = T.new(v)
        end
      end
    end

    def self.diag(nrows, ncolumns, values)
      new(nrows, ncolumns, 0, 0, MatrixFlags.for_diag(nrows == ncolumns)) do |i, j|
        values[i]
      end
    end

    def self.new(matrix : Matrix, tolerance = matrix.tolerance)
      upper_band = matrix.ncolumns - 1
      (matrix.ncolumns - 1).to(1) do |i|
        break unless matrix.diag(i).all? { |v| v.abs <= tolerance }
        upper_band = i - 1
      end
      lower_band = matrix.nrows - 1
      (matrix.nrows - 1).to(1) do |i|
        break unless matrix.diag(-i).all? { |v| v.abs <= tolerance }
        lower_band = i - 1
      end
      new(matrix.nrows, matrix.ncolumns, upper_band, lower_band, matrix.flags) do |i, j|
        matrix.unsafe_fetch(i, j)
      end
    end

    def unsafe_fetch(i, j)
      if index = ij2index(i, j)
        @raw_banded.unsafe_fetch(index)
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
      BandedMatrix(T).new(nrows, ncolumns, upper_band, lower_band) { |i, j| yield(unsafe_fetch(i, j), i, j) }
    end

    def map_with_index_f64(&block)
      BandedMatrix(Float64).new(nrows, ncolumns, upper_band, lower_band) { |i, j| yield(unsafe_fetch(i, j), i, j) }
    end

    def map_with_index_complex(&block)
      BandedMatrix(Complex).new(nrows, ncolumns, upper_band, lower_band) { |i, j| yield(unsafe_fetch(i, j), i, j) }
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
        return false if other.unsafe_fetch(i, j) != unsafe_fetch(i, j)
      end
      true
    end

    def transpose!
      return self if flags.symmetric?
      newraw = Slice(T).new(band_len*(@upper_band + @lower_band + 1), T.new(0))
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

    private def resized_flags
      if @upper_band == 0 && @lower_band == 0
        @flags = MatrixFlags.for_diag(nrows == ncolumns)
      elsif @upper_band == 0
        @flags = MatrixFlags::LowerTriangular
      elsif @lower_band == 0
        @flags = MatrixFlags::UpperTriangular
      else
        @flags = MatrixFlags::None
      end
    end

    def set_bands(aupper, alower) : Nil
      return if aupper == @upper_band && alower == @lower_band
      raise ArgumentError.new "upper_band must be non-negative" unless aupper >= 0
      raise ArgumentError.new "lower_band must be non-negative" unless alower >= 0
      newraw = Slice(T).new(band_len*(aupper + alower + 1), T.new(0))
      each_with_index do |v, i, j|
        next unless {0, j - aupper}.max <= i <= {nrows - 1, j + alower}.min
        ai = aupper + i - j
        newraw[ai*band_len + j] = v
      end
      @raw_banded = newraw
      @upper_band, @lower_band = aupper, alower
      resized_flags
    end

    def upper_band=(value)
      set_bands(value, @lower_band)
    end

    def lower_band=(value)
      set_bands(@upper_band, value)
    end

    def tril!(k = 0)
      if k >= 0
        # just update upper_band
        self.upper_band = {@upper_band, k}.min
      else
        # set upper band and cleanup lower elements
        self.upper_band = 0
        each_with_index { |v, i, j| unsafe_set(i, j, T.new(0)) if i < j - k }
        resized_flags
      end
      self
    end

    def triu!(k = 0)
      if k <= 0
        # just update lower_band
        self.lower_band = {@lower_band, -k}.min
      else
        # set lower band and cleanup upper elements
        self.lower_band = 0
        each_with_index { |v, i, j| unsafe_set(i, j, T.new(0)) if i > j - k }
        resized_flags
      end
      self
    end

    # returns element-wise sum
    def +(m : BandedMatrix(T))
      assert_same_size(m)
      result = BandedMatrix(T).new(nrows, ncolumns, {upper_band, m.upper_band}.max, {lower_band, m.lower_band}.max, flags.sum(m.flags)) do |i, j|
        self.unsafe_fetch(i, j) + m.unsafe_fetch(i, j)
      end
    end

    # returns element-wise subtract
    def -(m : BandedMatrix(T))
      assert_same_size(m)
      result = BandedMatrix(T).new(nrows, ncolumns, {upper_band, m.upper_band}.max, {lower_band, m.lower_band}.max, flags.sum(m.flags)) do |i, j|
        self.unsafe_fetch(i, j) - m.unsafe_fetch(i, j)
      end
    end

    # returns transposed matrix
    def transpose
      return clone if flags.symmetric?
      BandedMatrix(T).new(ncolumns, nrows, lower_band, upper_band, flags.transpose) do |i, j|
        unsafe_fetch(j, i)
      end
    end

    # returns conjtransposed matrix
    def conjtranspose
      {% if T != Complex %}
        return transpose
      {% end %}
      return clone if flags.hermitian?
      BandedMatrix(T).new(ncolumns, nrows, lower_band, upper_band, flags.transpose) do |i, j|
        unsafe_fetch(j, i).conj
      end
    end

    def add!(k : Number, m : BandedMatrix)
      assert_same_size(m)
      set_bands({upper_band, m.upper_band}.max, {lower_band, m.lower_band}.max)
      oldflags = flags
      map_with_index! { |v, i, j| v + k*m.unsafe_fetch(i, j) }
      self.flags = oldflags.sum(m.flags.scale(k.is_a?(Complex) && k.imag != 0))
      self
    end

    def add!(k : Number, m : Matrix)
      raise ArgumentError.new "can't add! non-banded matrix"
    end

    def self.rand(nrows, ncolumns, upper_band : Int32, lower_band : Int32, rng : Random = Random::DEFAULT)
      new(nrows, ncolumns, upper_band, lower_band) { |i, j| rng.rand }
    end

    def self.rand(nrows, ncolumns, upper_band : Int32, rng : Random = Random::DEFAULT)
      rand(nrows, ncolumns, upper_band, upper_band, rng)
    end

    def tril(k = 0)
      if k >= 0
        # just update upper_band
        return clone if k >= upper_band
        self.class.new(nrows, ncolumns, k, lower_band) do |i, j|
          unsafe_fetch(i, j)
        end
      else
        # set upper band and cleanup lower elements
        self.class.new(nrows, ncolumns, 0, lower_band) do |i, j|
          i >= j - k ? unsafe_fetch(i, j) : 0
        end
      end
    end

    def triu(k = 0)
      if k <= 0
        # just update lower_band
        return clone if -k >= lower_band
        self.class.new(nrows, ncolumns, upper_band, -k) do |i, j|
          unsafe_fetch(i, j)
        end
      else
        # set lower band and cleanup upper elements
        self.class.new(nrows, ncolumns, upper_band, 0) do |i, j|
          i <= j - k ? unsafe_fetch(i, j) : 0
        end
      end
    end

    # def kron(b : Matrix(T))
    # def cat(other : Matrix(T), dimension)

    # returns matrix norm
    def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
      # TODO - check if not square
      let = case kind
            when .frobenius?
              'F'
            when .one?
              'o'
            when .inf?
              'I'
            else
              'M'
            end.ord.to_u8

      worksize = kind.inf? ? nrows : 0
      lapack_util(langb, worksize, let, @nrows, @lower_band, @upper_band, matrix(self), @lower_band + @upper_band + 1)
    end

    def det(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      if flags.triangular?
        return diag.product
      end
      lru = overwrite_a ? self : self.clone
      lru.upper_band = @lower_band + @upper_band
      ipiv = Slice(Int32).new(nrows)
      pp! @lower_band, @upper_band, lru, lru.lower_band, lru.upper_band, lru.raw_banded
      lapack(gbtrf, nrows, nrows, @lower_band, @upper_band, lru, 2*@lower_band + @upper_band + 1, ipiv)
      lru.clear_flags
      lru.diag.product
    end
  end

  alias BMat = BandedMatrix(Float64)
  alias BMat32 = BandedMatrix(Float32)
  alias BMatComplex = BandedMatrix(Complex)
end
