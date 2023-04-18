require "./matrix"
require "./submatrix"

module LA
  abstract class Matrix(T)
    private macro def_indexable(name, offset, size)
      # Indexable(SubMatrix(T)) that allows iterating over {{name.id}}
      struct ::LA::Utils::{{name.id.capitalize}}(T)
        include Indexable(SubMatrix(T))
        # Creates {{name.id.capitalize}}(T) from matrix `base`
        protected def initialize(@base : Matrix(T))
        end
        # :nodoc:
        def size
          @base.n{{name.id}}
        end
        # :nodoc:
        def unsafe_fetch(index)
          # TODO unsafe_new for submatrix?
          SubMatrix(T).new(@base, {{offset}}, {{size}})
        end
      end

      # Returns Indexable(SubMatrix(T)) that allows iterating over {{name.id}}
      def {{name.id}}
        {{name.id.capitalize}}(T).new(self)
      end
    end

    def_indexable(columns, {0, index}, {@base.nrows, 1})
    def_indexable(rows, {index, 0}, {1, @base.ncolumns})

    # TODO - more macro magic?

    struct ::LA::Utils::Columns(T)
      # Returns `SubMatrix` that consist of given columns
      #
      # Example:
      # ```
      # m = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
      # m.columns[1...3].should eq m[0..2, 1..2]
      # ```
      def [](range : Range)
        start = range.begin || 0
        start += @base.ncolumns if start < 0

        aend = range.end || @base.ncolumns
        aend += @base.ncolumns if aend < 0

        if range.end && !range.excludes_end? # compensate if user wanted the endpoint only if not nil
          aend += 1
        end

        size = aend - start

        SubMatrix(T).new(@base, {0, start}, {@base.nrows, size})
      end
    end

    struct ::LA::Utils::Rows(T)
      # Returns `SubMatrix` that consist of given rows
      #
      # Example:
      # ```
      # m = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
      # m.rows[1..2].should eq m[1..2, 0..3]
      # ```
      def [](range : Range)
        start = range.begin || 0
        start += @base.nrows if start < 0

        aend = range.end || @base.nrows
        aend += @base.nrows if aend < 0

        if range.end && !range.excludes_end? # compensate if user wanted the endpoint only if not nil
          aend += 1
        end

        size = aend - start

        SubMatrix(T).new(@base, {start, 0}, {size, @base.ncolumns})
      end
    end

    # :nodoc:
    struct Diagonal(T)
      include Indexable(T)

      # Creates `Indexable(T)` that allows iterating over `offset`-th diagonal of `base` matrix
      protected def initialize(@base : Matrix(T), @offset = 0)
        raise ArgumentError.new("Offset #{offset} is too big (matrix size #{@base.nrows}x#{@base.ncolumns})") if size <= 0
      end

      def size
        if @offset >= 0
          {@base.nrows, @base.ncolumns - @offset}.min
        else
          {@base.nrows + @offset, @base.ncolumns}.min
        end
      end

      def unsafe_fetch(index)
        if @offset >= 0
          @base.unsafe_fetch(index, index + @offset)
        else
          @base.unsafe_fetch(index - @offset, index)
        end
      end
    end

    # Returns `Indexable(T)` that allow iterating over k-th diagonal of matrix
    #
    # Example:
    # ```
    # m = GMat32.new([[-1, 2, 3, 4],
    #                 [5, -6, 7, 8],
    #                 [9, 10, -11, 12]])
    # m.diag(0).to_a.should eq [-1, -6, -11]
    # m.diag(1).to_a.should eq [2, 7, 12]
    # m.diag(2).to_a.should eq [3, 8]
    # m.diag(3).to_a.should eq [4]
    # expect_raises(ArgumentError) { m.diag(4) }
    # m.diag(-1).to_a.should eq [5, 10]
    # m.diag(-2).to_a.should eq [9]
    # expect_raises(ArgumentError) { m.diag(-3) }
    # expect_raises(ArgumentError) { m.diag(-4) }
    # ```
    def diag(offset = 0)
      Diagonal(T).new(self, offset)
    end
  end
end
