require "./matrix"
require "./submatrix"

module LA
  module Matrix(T)
    private macro def_indexable(name, offset, size)
      struct {{name.id.capitalize}}(T)
        include Indexable(SubMatrix(T))
        protected def initialize(@base : Matrix(T))
        end
        def size
          @base.n{{name.id}}
        end
        # unsafe_new for submatrix?
        def unsafe_at(i)
          SubMatrix(T).new(@base, {{offset}}, {{size}})
        end
      end
      def {{name.id}}
        {{name.id.capitalize}}(T).new(self)
      end
    end

    def_indexable(columns, {0, i}, {@base.nrows, 1})
    def_indexable(rows, {i, 0}, {1, @base.ncolumns})

    # TODO - more macro magic?
    struct Columns(T)
      def [](range : Range)
        start = range.begin + (range.begin < 0 ? @base.ncolumns : 0)
        size = range.end - start + (range.end < 0 ? @base.ncolumns : 0)
        SubMatrix(T).new(@base, {0, start}, {@base.nrows, size})
      end
    end

    struct Rows(T)
      def [](range : Range)
        size = range.size
        size += @base.nrows if range.end < 0
        SubMatrix(T).new(@base, {range.begin, 0}, {size, @base.ncolumns})
      end
    end

    struct Diagonal(T)
      include Indexable(T)

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

      # unsafe_new for submatrix?
      def unsafe_at(i)
        if @offset >= 0
          @base.unsafe_at(i, i + @offset)
        else
          @base.unsafe_at(i - @offset, i)
        end
      end
    end

    def diag(offset = 0)
      Diagonal(T).new(self, offset)
    end
  end
end
