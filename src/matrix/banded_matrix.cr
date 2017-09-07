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

    def initialize(@nrows, @ncolumns, @upper_band, @lower_band, values : Indexable)
      check_type
      @raw_banded = Slice(T).new(bands_size) do |index|
        i = index % band_len
        band = index/band_len
        T.new(values[band][i])
      end
    end

    # def self.from(matrix)
    #   new(matrix.nrows, matrix.ncolumns, matrix.raw)
    # end

    # def initialize(@nrows, @ncolumns, values : Indexable, @flags = MatrixFlags::None)
    #   check_type
    #   @raw = Slice(T).new(nrows*ncolumns) { |i| T.new(values[i]) }
    # end

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

    # def dup
    #   GeneralMatrix(T).new(@nrows, @ncolumns, @raw, @flags)
    # end

    # def clone
    # end

    # def to_unsafe
    # end
  end

  alias BMat = BandedMatrix(Float64)
  alias BMat32 = BandedMatrix(Float32)
  alias BMatComplex = BandedMatrix(Complex)
end
