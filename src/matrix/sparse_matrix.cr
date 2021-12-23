require "./matrix"
require "./submatrix"

module LA
  module Sparse
    abstract class Matrix(T) < LA::Matrix(T)
      getter nrows : Int32 = 0
      getter ncolumns : Int32 = 0
      property flags : MatrixFlags = MatrixFlags::None

      abstract def nonzeros : Int32
    end

    class COOMatrix(T) < Matrix(T)
      protected getter raw_columns : Array(Int32)
      protected getter raw_rows : Array(Int32)
      protected getter raw_values : Array(T)
      protected getter dictionary = {} of Tuple(Int32, Int32) => Int32

      def initialize(@nrows, @ncolumns, capacity = 0)
        @raw_rows = Array(Int32).new(capacity)
        @raw_columns = Array(Int32).new(capacity)
        @raw_values = Array(T).new(capacity)
      end

      def initialize(@nrows, @ncolumns, rows : Array(Int32), columns : Array(Int32), values : Array(T), @flags = MatrixFlags::None, *, dont_clone : Bool = false, dont_check : Bool = false, dictionary = nil)
        if rows.size != columns.size
          raise ArgumentError.new("Can't construct COO Matrix from arrays of different size: rows.size=#{rows.size} != columns.size=#{columns.size}")
        end
        if rows.size != values.size
          raise ArgumentError.new("Can't construct COO Matrix from arrays of different size: rows.size=#{rows.size} != values.size=#{values.size}")
        end
        if dont_clone
          @raw_rows = rows
          @raw_columns = columns
          @raw_values = values
        else
          @raw_rows = rows.dup
          @raw_columns = columns.dup
          @raw_values = values.dup
        end
        if dictionary
          @dictionary = dictionary.dup
        else
          @dictionary = {} of Tuple(Int32, Int32) => Int32
          nonzeros.times do |i|
            row = @raw_rows[i]
            unless (0...@nrows).includes?(row)
              raise ArgumentError.new("Can't construct COO Matrix from arrays: rows[#{i}]=#{@raw_rows[i]} is outside (#{0}..#{@nrows})")
            end
            column = @raw_columns[i]
            unless (0...@ncolumns).includes?(column)
              raise ArgumentError.new("Can't construct COO Matrix from arrays: columns[#{i}]=#{@raw_columns[i]} is outside (#{0}..#{@ncolumns})")
            end
            @dictionary[{row, column}] = i
          end
        end
      end

      def nonzeros : Int32
        @raw_rows.size
      end

      def self.new(matrix : COOMatrix(T))
        new(matrix.nrows, matrix.ncolumns, matrix.raw_rows, matrix.raw_columns, matrix.raw_values, dictionary: matrix.dictionary)
      end

      def self.new(matrix : COOMatrix)
        new(matrix.nrows, matrix.ncolumns, matrix.raw_rows.dup, matrix.raw_columns.dup, matrix.raw_values.map { |v| T.new(v) }, dont_clone: true, dictionary: matrix.dictionary)
      end

      def self.new(matrix : LA::Matrix)
        if matrix.is_a? Sparse::Matrix
          result = new(matrix.nrows, matrix.ncolumns, matrix.nonzeros)
        else
          result = new(matrix.nrows, matrix.ncolumns)
        end
        matrix.each_with_index(all: false) do |v, i, j|
          result.push_element(i, j, v)
        end
        result
      end

      def unsafe_fetch(i, j)
        if index = @dictionary[{i, j}]?
          @raw_values[index]
        else
          T.new(0.0)
        end
      end

      protected def push_element(i, j, v)
        @raw_rows << i
        @raw_columns << j
        @raw_values << T.new(v)
        @dictionary[{i, j}] = @raw_rows.size - 1
      end

      def add_element(i, j, v)
        clear_flags # TODO - not always?
        if index = @dictionary[{i, j}]?
          @raw_values[index] += T.new(value)
        else
          push_element(i, j, v)
        end
      end

      def unsafe_set(i, j, value)
        clear_flags # TODO - not always?
        if index = @dictionary[{i, j}]?
          @raw_values[index] = T.new(value)
        else
          # TODO - delete elements?
          push_element(i, j, value)
        end
      end

      def dup
        COOMatrix(T).new(self)
      end

      def clone
        dup
      end
    end
  end
end
