require "./matrix"
require "./submatrix"

module Linalg
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
        SubMatrix(T).new(@base, {0, range.begin}, {@base.nrows, range.size})
      end
    end

    struct Rows(T)
      def [](range : Range)
        SubMatrix(T).new(@base, {range.begin, 0}, {range.size, @base.ncolumns})
      end
    end
  end
end
