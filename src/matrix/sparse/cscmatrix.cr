require "./sparse_matrix.cr"

# TODO - inline docs

module LA::Sparse
  macro csxmatrix(name, arows, acolumns, index1, index2, norm)
    class {{name}}Matrix(T) < Matrix(T)
        protected getter raw_{{acolumns}} : Array(Int32)
        protected getter raw_{{arows}} : Array(Int32)
        protected getter raw_values : Array(T)

        def initialize(@nrows, @ncolumns, capacity = 0)
        @raw_{{arows}} = Array(Int32).new(n{{arows}} + 1, 0)
        @raw_{{acolumns}} = Array(Int32).new(capacity)
        @raw_values = Array(T).new(capacity)
        @flags = MatrixFlags.for_diag(@n{{arows}} == @n{{acolumns}})
        end

        def initialize(@nrows, @ncolumns, raw_{{arows}} : Array(Int32), raw_{{acolumns}} : Array(Int32), raw_values : Array(T), @flags = MatrixFlags::None, *, dont_clone : Bool = false)
        if raw_{{arows}}.size != @n{{arows}} + 1
            raise ArgumentError.new("Can't construct #{self.class} from arrays of different size: {{arows}}.size(#{raw_{{arows}}.size}) != n{{arows}}+1 (#{@n{{arows}} + 1}")
        end
        if raw_{{acolumns}}.size != raw_values.size
            raise ArgumentError.new("Can't construct #{self.class} from arrays of different size: {{acolumns}}.size(#{raw_{{acolumns}}.size}) != values.size(#{raw_values.size})")
        end
        if dont_clone
            @raw_{{arows}} = raw_{{arows}}
            @raw_{{acolumns}} = raw_{{acolumns}}
            @raw_values = raw_values
        else
            @raw_{{arows}} = raw_{{arows}}.dup
            @raw_{{acolumns}} = raw_{{acolumns}}.dup
            @raw_values = raw_values.dup
        end
        end

        def nonzeros : Int32
          @raw_values.size
        end

        def self.new(matrix : self)
          new(matrix.nrows, matrix.ncolumns, matrix.raw_{{arows}}, matrix.raw_{{acolumns}}, matrix.raw_values, flags: matrix.flags)
        end

        def self.new(matrix : {{name}}Matrix)
          new(matrix.nrows, matrix.ncolumns, matrix.raw_{{arows}}.dup, matrix.raw_{{acolumns}}.dup, matrix.raw_values.map { |v| T.new(v) }, dont_clone: true, flags: matrix.flags)
        end

        def self.new(matrix : LA::Matrix)
        if matrix.is_a? Sparse::Matrix
            nonzeros = matrix.nonzeros
        else
            nonzeros = matrix.count { |v| !v.zero? }
        end
        result = new(matrix.nrows, matrix.ncolumns, nonzeros)
        index = 0
        result.n{{arows}}.times do |i|
            result.n{{acolumns}}.times do |j|
            v = matrix.unsafe_fetch({{index1}}, {{index2}})
            next if v.zero?
            result.raw_values << T.new(v)
            result.raw_{{acolumns}} << j
            index += 1
            end
            result.raw_{{arows}}[i + 1] = index
        end
        result.flags = matrix.flags
        result
        end

        private def ij2index(i, j) : Int32?
        row_start = @raw_{{arows}}[i]
        row_end = @raw_{{arows}}[i + 1]
        index = (row_start...row_end).bsearch { |n| @raw_{{acolumns}}[n] >= j }
        return nil unless index
        @raw_{{acolumns}}[index] == j ? index : nil
        end

        def unsafe_fetch(i, j)
        if index = ij2index({{index1}}, {{index2}})
            @raw_values.unsafe_fetch(index)
        else
            T.new(0.0)
        end
        end

        def unsafe_set(i, j, value)
        if index = ij2index({{index1}}, {{index2}})
            @raw_values.unsafe_put(index, T.new(value))
        else
            raise "cannot add values to #{self.class}"
        end
        end

        # Returns diagonal matrix of given size with diagonal elements taken from array `values`
        #
        # Raises if `values.size > {n{{arows}}, n{{acolumns}}}.min`
        def self.diag(nrows, ncolumns, values)
        raise ArgumentError.new("Too much elements (#{values.size}) for diag matrix [#{nrows}x#{ncolumns}]") if values.size > {nrows, ncolumns}.min
        n = values.size
        new(nrows, ncolumns,
            Array(Int32).new(n{{arows}} + 1) { |i| i <= n ? i : n },
            Array(Int32).new(n) { |i| i },
            Array(T).new(n) { |i| T.new(values[i]) },
            MatrixFlags.for_diag(n{{arows}} == n{{acolumns}}),
            dont_clone: true
        )
        end

        # Returns diagonal matrix of given size with diagonal elements equal to block value
        def self.diag(nrows, ncolumns, &block)
        n = {n{{arows}}, n{{acolumns}}}.min
        new(nrows, ncolumns,
            Array(Int32).new(n{{arows}} + 1) { |i| i <= n ? i : n },
            Array(Int32).new(n) { |i| i },
            Array(T).new(n) { |i| T.new(yield(i)) },
            MatrixFlags.for_diag(n{{arows}} == n{{acolumns}}),
            dont_clone: true
        )
        end

        def each_index(*, all = false, &block)
        if all
            super(all: true) { |i, j| yield(i, j) }
        else
            @n{{arows}}.times do |i|
            row_start = @raw_{{arows}}[i]
            row_end = @raw_{{arows}}[i + 1]
            (row_start...row_end).each do |index|
                j = @raw_{{acolumns}}[index]
                yield({{index1}}, {{index2}})
            end
            end
        end
        end

        def each_with_index(*, all = false, &block)
        if all
            super(all: true) { |v, i, j| yield(v, i, j) }
        else
            @n{{arows}}.times do |i|
            row_start = @raw_{{arows}}[i]
            row_end = @raw_{{arows}}[i + 1]
            (row_start...row_end).each do |index|
                j = @raw_{{acolumns}}[index]
                yield(@raw_values[index], {{index1}}, {{index2}})
            end
          end
        end
        end

        def map_with_index!(&block)
          @n{{arows}}.times do |i|
            row_start = @raw_{{arows}}[i]
            row_end = @raw_{{arows}}[i + 1]
            (row_start...row_end).each do |index|
              j = @raw_{{acolumns}}[index]
              @raw_values[index] = T.new(yield(@raw_values[index], {{index1}}, {{index2}}))
            end
          end
        end

        def map_with_index(&block)
        row = 0
        values = @raw_values.map_with_index do |v, index|
          while index >= @raw_{{arows}}[row + 1]
            row += 1
          end
          i,j = row, @raw_{{acolumns}}[index]
          newv = yield(v, {{index1}}, {{index2}})
          newv
        end
        {{name}}Matrix(T).new(@nrows, @ncolumns, @raw_{{arows}}.dup, @raw_{{acolumns}}.dup, values, dont_clone: true)
        # clone.tap &.map_with_index! { |v, i, j| yield(v, i, j) }
        end

        def map_with_index_f64(&block)
          row = 0
          values = @raw_values.map_with_index do |v, index|
            while index >= @raw_{{arows}}[row + 1]
              row += 1
            end
            i,j = row, @raw_{{acolumns}}[index]
            newv = yield(v, {{index1}}, {{index2}})
            Float64.new(newv)
          end
          {{name}}Matrix(Float64).new(@nrows, @ncolumns, @raw_{{arows}}.dup, @raw_{{acolumns}}.dup, values, dont_clone: true)
        end

        def map_with_index_complex(&block)
          row = 0
          values = @raw_values.map_with_index do |v, index|
            while index >= @raw_{{arows}}[row + 1]
              row += 1
            end
            i,j = row, @raw_{{acolumns}}[index]
            newv = yield(v, {{index1}}, {{index2}})
            Complex.new(newv)
          end
          {{name}}Matrix(Complex).new(@nrows, @ncolumns, @raw_{{arows}}.dup, @raw_{{acolumns}}.dup, values, dont_clone: true)
        end

        # def transpose! TODO

        def transpose
        return clone if flags.symmetric?
        result = self.class.new(ncolumns, nrows, raw_{{arows}}: Array(Int32).new(n{{acolumns}} + 1, 0), raw_{{acolumns}}: Array(Int32).new(nonzeros, 0), raw_values: Array(T).new(nonzeros, T.new(0)))
        t_row_pos = Slice(Int32).new(n{{acolumns}}, 0)
        each_index do |i,j|
            t_row_pos[{{index2}}] += 1
        end
        n{{acolumns}}.times do |col|
            result.raw_{{arows}}[col + 1] = result.raw_{{arows}}[col] + t_row_pos[col]
        end
        t_row_pos.fill(0)
        each_with_index do |v, i, j|
            index = t_row_pos[{{index2}}] + result.raw_{{arows}}[{{index2}}]
            result.raw_values[index] = v
            result.raw_{{acolumns}}[index] = {{index1}}
            t_row_pos[{{index2}}] += 1
        end
        result.flags = self.flags.transpose
        result
        end

        def conjtranspose
        \{% if T != Complex %}
            return transpose
        \{% end %}
        return clone if flags.hermitian?
        transpose.map!(&.conj).tap { |r| r.flags = self.flags.transpose }
        end

        def add(m : self, *, alpha = 1, beta = 1)
        assert_same_size(m)
        result = self.class.new(nrows, ncolumns, nonzeros + m.nonzeros)
        row = Slice(T?).new(n{{acolumns}}, value: nil.as(T?))
        cols = Array(Int32).new(n{{acolumns}})
        n{{arows}}.times do |i|
            row.fill(nil.as(T?))
            cols.clear
            self.scatter(row, alpha, i, cols)
            m.scatter(row, beta, i, cols)
            cols.sort!
            cols.each do |col|
            v = row[col]
            next unless v
            next if v.zero?
            result.raw_{{acolumns}} << col
            result.raw_values << v
            end
            result.raw_{{arows}}[i + 1] = result.raw_values.size
        end
        result.flags = flags.add(m.flags, alpha, beta)
        result
        end

        def clear
        @raw_{{arows}}.fill(0)
        @raw_values.clear
        @raw_{{acolumns}}.clear
        end

        # def self.rand(n{{arows}}, n{{acolumns}}, *, nonzero_elements, rng : Random = Random::DEFAULT)
        # def self.rand(n{{arows}}, n{{acolumns}}, *, fill_factor, rng : Random = Random::DEFAULT)
        def select!(& : T -> Bool)
          select_with_index! { |v, i, j| yield(v) }
        end

        def select_index!(& : (Int32, Int32) -> Bool)
          select_with_index! { |v, i, j| yield(i, j) }
        end

        def select_with_index!(& : (T, Int32, Int32) -> Bool)
        new_index = 0
        delta = 0
        row = 0
        n = @raw_values.size
        n.times do |index|
            while index >= @raw_{{arows}}[row + 1]
            @raw_{{arows}}[row + 1] -= delta
            row += 1
            end
            col = @raw_{{acolumns}}[index]
            v = @raw_values[index]
            i, j = row, col
            stay = yield(v, {{index1}}, {{index2}})
            if stay
            @raw_{{acolumns}}[new_index] = col
            @raw_values[new_index] = v
            new_index += 1
            else
            delta += 1
            end
        end
        (row + 1..@raw_{{arows}}.size - 1).each do |r|
            @raw_{{arows}}[r] = new_index
        end
        @raw_{{acolumns}}.truncate(0, new_index)
        @raw_values.truncate(0, new_index)
        end

        def select_with_index(& : (T, Int32, Int32) -> Bool)
        result = self.class.new(n{{arows}}, n{{acolumns}}, nonzeros)
        row = 0
        n = @raw_values.size
        n.times do |index|
            while index >= @raw_{{arows}}[row + 1]
            result.raw_{{arows}}[row + 1] = result.raw_values.size
            row += 1
            end
            col = @raw_{{acolumns}}[index]
            v = @raw_values[index]
            i,j = row, col
            stay = yield(v, {{index1}}, {{index2}})
            if stay
            result.raw_{{acolumns}} << col
            result.raw_values << v
            end
        end
        total = result.raw_values.size
        (row + 1..@raw_{{arows}}.size - 1).each do |r|
            result.raw_{{arows}}[r] = total
        end
        result
        end

        def resize!(arows, acolumns)
        return if arows == self.nrows && acolumns == self.ncolumns
        if arows < self.nrows || acolumns < self.ncolumns
            select_index! do |i, j|
            i < arows && j < acolumns
            end
        end
        if a{{arows}} < @n{{arows}}
            @raw_{{arows}}.truncate(0, a{{arows}} + 1)
        else
            n = nonzeros
            (a{{arows}} - @n{{arows}}).times { @raw_{{arows}} << n }
        end
        @n{{arows}} = a{{arows}}
        @n{{acolumns}} = a{{acolumns}}
        clear_flags
        self
        end

        protected def scatter(vector : Slice(T?), scale : T, row : Int32, cols : Array(Int32))
        row_start = raw_{{arows}}[row]
        row_end = raw_{{arows}}[row + 1]
        (row_start...row_end).each do |index|
            col = raw_{{acolumns}}[index]
            v = raw_values[index]*scale
            old = vector[col]
            if old.nil?
            vector[col] = v
            cols << col
            else
            vector[col] = old + v
            end
        end
        end

        def *(b : self)
        if ncolumns != b.nrows
            raise ArgumentError.new("matrix size should match (#{shape_str} * #{b.shape_str}")
        end
        c = self.class.new(nrows, b.ncolumns, nonzeros + b.nonzeros)
        i, j = self, b
        a, b = {{index1}}, {{index2}}
        row_ci = Slice(T?).new(b.n{{acolumns}}, value: nil.as(T?))
        ci_cols = Array(Int32).new(b.n{{acolumns}})
        a.n{{arows}}.times do |i|
            row_ci.fill(nil.as(T?))
            ci_cols.clear
            row_start = a.raw_{{arows}}[i]
            row_end = a.raw_{{arows}}[i + 1]
            (row_start...row_end).each do |index|
            k = a.raw_{{acolumns}}[index]
            aik = a.raw_values[index]
            # aik x bk* = ci*
            b.scatter(row_ci, aik, k, ci_cols)
            end
            # copy row to matrix c
            ci_cols.sort!
            ci_cols.each do |col|
            v = row_ci[col]
            next unless v
            next if v.zero?
            c.raw_{{acolumns}} << col
            c.raw_values << v
            end
            c.raw_{{arows}}[i + 1] = c.raw_values.size
        end
        c.flags = a.flags.mult(b.flags)
        c
        end

        def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
        result = T.additive_identity
        case kind
        when .{{norm}}?
            n{{arows}}.times do |i|
            sum = T.additive_identity
            row_start = @raw_{{arows}}[i]
            row_end = @raw_{{arows}}[i + 1]
            (row_start...row_end).each do |index|
                sum += @raw_values[index].abs
            end
            result = sum if result < sum
            end
            result
        else
            super(kind)
        end
        end
    end
  end

  csxmatrix(CSR, rows, columns, i, j, inf)
  csxmatrix(CSC, columns, rows, j, i, one)
end
