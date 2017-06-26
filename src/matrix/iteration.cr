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
  end
end
