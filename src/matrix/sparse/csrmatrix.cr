require "./sparse_matrix.cr"

# TODO - inline docs

module LA::Sparse
  # class CSRMatrix(T) < Matrix(T)
  #   protected getter raw_columns : Array(Int32)
  #   protected getter raw_rows : Array(Int32)
  #   protected getter raw_values : Array(T)

  #   def initialize(@nrows, @ncolumns, capacity = 0)
  #     @raw_rows = Array(Int32).new(nrows + 1, 0)
  #     @raw_columns = Array(Int32).new(capacity)
  #     @raw_values = Array(T).new(capacity)
  #     @flags = MatrixFlags.for_diag(@nrows == @ncolumns)
  #   end

  #   def initialize(@nrows, @ncolumns, raw_rows : Array(Int32), raw_columns : Array(Int32), raw_values : Array(T), @flags = MatrixFlags::None, *, dont_clone : Bool = false)
  #     if raw_rows.size != @nrows + 1
  #       raise ArgumentError.new("Can't construct CSR Matrix from arrays of different size: rows.size(#{raw_rows.size}) != nrows+1 (#{@nrows + 1}")
  #     end
  #     if raw_columns.size != raw_values.size
  #       raise ArgumentError.new("Can't construct CSR Matrix from arrays of different size: columns.size(#{raw_columns.size}) != values.size(#{raw_values.size})")
  #     end
  #     if dont_clone
  #       @raw_rows = raw_rows
  #       @raw_columns = raw_columns
  #       @raw_values = raw_values
  #     else
  #       @raw_rows = raw_rows.dup
  #       @raw_columns = raw_columns.dup
  #       @raw_values = raw_values.dup
  #     end
  #   end

  #   def nonzeros : Int32
  #     @raw_values.size
  #   end

  #   def self.new(matrix : CSRMatrix(T))
  #     new(matrix.nrows, matrix.ncolumns, matrix.raw_rows, matrix.raw_columns, matrix.raw_values, flags: matrix.flags)
  #   end

  #   def self.new(matrix : CSRMatrix)
  #     new(matrix.nrows, matrix.ncolumns, matrix.raw_rows.dup, matrix.raw_columns.dup, matrix.raw_values.map { |v| T.new(v) }, dont_clone: true, flags: matrix.flags)
  #   end

  #   def self.new(matrix : LA::Matrix)
  #     if matrix.is_a? Sparse::Matrix
  #       nonzeros = matrix.nonzeros
  #     else
  #       nonzeros = matrix.count { |v| !v.zero? }
  #     end
  #     result = new(matrix.nrows, matrix.ncolumns, nonzeros)
  #     index = 0
  #     result.nrows.times do |row|
  #       result.ncolumns.times do |column|
  #         v = matrix.unsafe_fetch(row, column)
  #         next if v.zero?
  #         result.raw_values << T.new(v)
  #         result.raw_columns << column
  #         index += 1
  #       end
  #       result.raw_rows[row + 1] = index
  #     end
  #     result.flags = matrix.flags
  #     result
  #   end

  #   private def ij2index(i, j) : Int32?
  #     row_start = @raw_rows[i]
  #     row_end = @raw_rows[i + 1]
  #     index = (row_start...row_end).bsearch { |n| @raw_columns[n] >= j }
  #     return nil unless index
  #     @raw_columns[index] == j ? index : nil
  #   end

  #   def unsafe_fetch(i, j)
  #     if index = ij2index(i, j)
  #       @raw_values.unsafe_fetch(index)
  #     else
  #       T.new(0.0)
  #     end
  #   end

  #   def unsafe_set(i, j, value)
  #     if index = ij2index(i, j)
  #       @raw_values.unsafe_put(index, T.new(value))
  #     else
  #       raise "cannot add values to CSR matrix"
  #     end
  #   end

  #   # Returns diagonal matrix of given size with diagonal elements taken from array `values`
  #   #
  #   # Raises if `values.size > {nrows, ncolumns}.min`
  #   def self.diag(nrows, ncolumns, values)
  #     raise ArgumentError.new("Too much elements (#{values.size}) for diag matrix [#{nrows}x#{ncolumns}]") if values.size > {nrows, ncolumns}.min
  #     n = values.size
  #     new(nrows, ncolumns,
  #       Array(Int32).new(nrows + 1) { |i| i <= n ? i : n },
  #       Array(Int32).new(n) { |i| i },
  #       Array(T).new(n) { |i| T.new(values[i]) },
  #       MatrixFlags.for_diag(nrows == ncolumns),
  #       dont_clone: true
  #     )
  #   end

  #   # Returns diagonal matrix of given size with diagonal elements equal to block value
  #   def self.diag(nrows, ncolumns, &block)
  #     n = {nrows, ncolumns}.min
  #     new(nrows, ncolumns,
  #       Array(Int32).new(nrows + 1) { |i| i <= n ? i : n },
  #       Array(Int32).new(n) { |i| i },
  #       Array(T).new(n) { |i| T.new(yield(i)) },
  #       MatrixFlags.for_diag(nrows == ncolumns),
  #       dont_clone: true
  #     )
  #   end

  #   def each_index(*, all = false, &block)
  #     if all
  #       super(all: true) { |i, j| yield(i, j) }
  #     else
  #       @nrows.times do |i|
  #         row_start = @raw_rows[i]
  #         row_end = @raw_rows[i + 1]
  #         (row_start...row_end).each do |index|
  #           yield(i, @raw_columns[index])
  #         end
  #       end
  #     end
  #   end

  #   def each_with_index(*, all = false, &block)
  #     if all
  #       super(all: true) { |v, i, j| yield(v, i, j) }
  #     else
  #       @nrows.times do |i|
  #         row_start = @raw_rows[i]
  #         row_end = @raw_rows[i + 1]
  #         (row_start...row_end).each do |index|
  #           yield(@raw_values[index], i, @raw_columns[index])
  #         end
  #       end
  #     end
  #   end

  #   def map_with_index!(&block)
  #     @nrows.times do |i|
  #       row_start = @raw_rows[i]
  #       row_end = @raw_rows[i + 1]
  #       (row_start...row_end).each do |index|
  #         @raw_values[index] = T.new(yield(@raw_values[index], i, @raw_columns[index]))
  #       end
  #     end
  #   end

  #   def map_with_index(&block)
  #     row = 0
  #     values = @raw_values.map_with_index do |v, i|
  #       while i >= @raw_rows[row + 1]
  #         row += 1
  #       end
  #       newv = yield(v, row, @raw_columns[i])
  #       newv
  #     end
  #     CSRMatrix(T).new(@nrows, @ncolumns, @raw_rows.dup, @raw_columns.dup, values, dont_clone: true)
  #     # clone.tap &.map_with_index! { |v, i, j| yield(v, i, j) }
  #   end

  #   # def transpose! TODO

  #   def transpose
  #     return clone if flags.symmetric?
  #     result = self.class.new(ncolumns, nrows, raw_rows: Array(Int32).new(ncolumns + 1, 0), raw_columns: Array(Int32).new(nonzeros, 0), raw_values: Array(T).new(nonzeros, T.new(0)))
  #     t_row_pos = Slice(Int32).new(ncolumns, 0)
  #     each_index do |row, col|
  #       t_row_pos[col] += 1
  #     end
  #     ncolumns.times do |col|
  #       result.raw_rows[col + 1] = result.raw_rows[col] + t_row_pos[col]
  #     end
  #     t_row_pos.fill(0)
  #     each_with_index do |v, row, col|
  #       index = t_row_pos[col] + result.raw_rows[col]
  #       result.raw_values[index] = v
  #       result.raw_columns[index] = row
  #       t_row_pos[col] += 1
  #     end
  #     result.flags = self.flags.transpose
  #     result
  #   end

  #   def conjtranspose
  #     {% if T != Complex %}
  #       return transpose
  #     {% end %}
  #     return clone if flags.hermitian?
  #     transpose.map!(&.conj).tap { |r| r.flags = self.flags.transpose }
  #   end

  #   def add(m : CSRMatrix(T), *, alpha = 1, beta = 1)
  #     assert_same_size(m)
  #     result = CSRMatrix(T).new(nrows, ncolumns, nonzeros + m.nonzeros)
  #     row = Slice(T?).new(ncolumns, value: nil.as(T?))
  #     cols = Array(Int32).new(ncolumns)
  #     nrows.times do |i|
  #       row.fill(nil.as(T?))
  #       cols.clear
  #       self.scatter(row, alpha, i, cols)
  #       m.scatter(row, beta, i, cols)
  #       cols.sort!
  #       cols.each do |col|
  #         v = row[col]
  #         next unless v
  #         next if v.zero?
  #         result.raw_columns << col
  #         result.raw_values << v
  #       end
  #       result.raw_rows[i + 1] = result.raw_values.size
  #     end
  #     result.flags = flags.add(m.flags, alpha, beta)
  #     result
  #   end

  #   def clear
  #     @raw_rows.fill(0)
  #     @raw_values.clear
  #     @raw_columns.clear
  #   end

  #   # def self.rand(nrows, ncolumns, *, nonzero_elements, rng : Random? = nil)
  #   # def self.rand(nrows, ncolumns, *, fill_factor, rng : Random? = nil)
  #   def select!(& : T -> Bool)
  #     select_with_index! { |v, i, j| yield(v) }
  #   end

  #   def select_index!(& : (Int32, Int32) -> Bool)
  #     select_with_index! { |v, i, j| yield(i, j) }
  #   end

  #   def select_with_index!(& : (T, Int32, Int32) -> Bool)
  #     new_index = 0
  #     delta = 0
  #     row = 0
  #     n = @raw_values.size
  #     n.times do |index|
  #       while index >= @raw_rows[row + 1]
  #         @raw_rows[row + 1] -= delta
  #         row += 1
  #       end
  #       col = @raw_columns[index]
  #       v = @raw_values[index]
  #       stay = yield(v, row, col)
  #       if stay
  #         @raw_columns[new_index] = col
  #         @raw_values[new_index] = v
  #         new_index += 1
  #       else
  #         delta += 1
  #       end
  #     end
  #     (row + 1..@raw_rows.size - 1).each do |r|
  #       @raw_rows[r] = new_index
  #     end
  #     @raw_columns.truncate(0, new_index)
  #     @raw_values.truncate(0, new_index)
  #   end

  #   def select_with_index(& : (T, Int32, Int32) -> Bool)
  #     result = self.class.new(nrows, ncolumns, nonzeros)
  #     row = 0
  #     n = @raw_values.size
  #     n.times do |index|
  #       while index >= @raw_rows[row + 1]
  #         result.raw_rows[row + 1] = result.raw_values.size
  #         row += 1
  #       end
  #       col = @raw_columns[index]
  #       v = @raw_values[index]
  #       stay = yield(v, row, col)
  #       if stay
  #         result.raw_columns << col
  #         result.raw_values << v
  #       end
  #     end
  #     total = result.raw_values.size
  #     (row + 1..@raw_rows.size - 1).each do |r|
  #       result.raw_rows[r] = total
  #     end
  #     result
  #   end

  #   def resize!(anrows, ancolumns)
  #     return if anrows == self.nrows && ancolumns == self.ncolumns
  #     if anrows < self.nrows || ancolumns < self.ncolumns
  #       select_index! do |i, j|
  #         i < anrows && j < ancolumns
  #       end
  #     end
  #     if anrows < @nrows
  #       @raw_rows.truncate(0, anrows + 1)
  #     else
  #       n = nonzeros
  #       (anrows - @nrows).times { @raw_rows << n }
  #     end
  #     @nrows = anrows
  #     @ncolumns = ancolumns
  #     clear_flags
  #     self
  #   end

  #   protected def scatter(vector : Slice(T?), scale : T, row : Int32, cols : Array(Int32))
  #     row_start = raw_rows[row]
  #     row_end = raw_rows[row + 1]
  #     (row_start...row_end).each do |index|
  #       col = raw_columns[index]
  #       v = raw_values[index]*scale
  #       old = vector[col]
  #       if old.nil?
  #         vector[col] = v
  #         cols << col
  #       else
  #         vector[col] = old + v
  #       end
  #     end
  #   end

  #   def *(b : self)
  #     if ncolumns != b.nrows
  #       raise ArgumentError.new("matrix size should match (#{shape_str} * #{b.shape_str}")
  #     end
  #     c = self.class.new(nrows, b.ncolumns, nonzeros + b.nonzeros)
  #     row_ci = Slice(T?).new(b.ncolumns, value: nil.as(T?))
  #     ci_cols = Array(Int32).new(b.ncolumns)
  #     nrows.times do |i|
  #       row_ci.fill(nil.as(T?))
  #       ci_cols.clear
  #       row_start = @raw_rows[i]
  #       row_end = @raw_rows[i + 1]
  #       (row_start...row_end).each do |index|
  #         k = @raw_columns[index]
  #         aik = @raw_values[index]
  #         # aik x bk* = ci*
  #         b.scatter(row_ci, aik, k, ci_cols)
  #       end
  #       # copy row to matrix c
  #       ci_cols.each do |col|
  #         v = row_ci[col]
  #         next unless v
  #         next if v.zero?
  #         c.raw_columns << col
  #         c.raw_values << v
  #       end
  #       c.raw_rows[i + 1] = c.raw_values.size
  #     end
  #     c.flags = self.flags.mult(b.flags)
  #     c
  #   end

  #   def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
  #     result = T.additive_identity
  #     case kind
  #     when .inf?
  #       nrows.times do |i|
  #         rowsum = T.additive_identity
  #         row_start = @raw_rows[i]
  #         row_end = @raw_rows[i + 1]
  #         (row_start...row_end).each do |index|
  #           rowsum += @raw_values[index].abs
  #         end
  #         result = rowsum if result < rowsum
  #       end
  #       result
  #     else
  #       super(kind)
  #     end
  #   end
  # end
end
