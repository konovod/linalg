require "./matrix"
require "./submatrix"

module LA
  module Sparse
    abstract class Matrix(T) < LA::Matrix(T)
      getter nrows : Int32 = 0
      getter ncolumns : Int32 = 0
      property flags : MatrixFlags = MatrixFlags::None
    end

    class COOMatrix(T) < Matrix(T)
      protected getter raw_columns : Array(Int32)
      protected getter raw_rows : Array(Int32)
      protected getter raw_values : Array(T)
      protected getter dictionary = {} of Tuple(Int32, Int32) => Int32

      def initialize(@nrows, @ncolumns)
        @raw_rows = Array(Int32).new
        @raw_columns = Array(Int32).new
        @raw_values = Array(T).new
      end

      def initialize(@nrows, @ncolumns, rows : Array(Int32), columns : Array(Int32), values : Array(T), @flags = MatrixFlags::None, *, dont_clone : Bool = false, dictionary = nil)
        raise "Can't construct COO Matrix from arrays of different size: rows.size=#{rows.size} != columns.size=#{columns.size}" if rows.size != columns.size
        raise "Can't construct COO Matrix from arrays of different size: rows.size=#{rows.size} != values.size=#{values.size}" if rows.size != values.size
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
          n = @raw_rows.size
          n.times do |i|
            @dictionary[{@raw_rows[i], @raw_columns[i]}] = i
          end
        end
      end

      def self.new(matrix : COOMatrix(T))
        new(matrix.nrows, matrix.ncolumns, matrix.raw_rows, matrix.raw_columns, matrix.raw_values, dictionary: matrix.dictionary)
      end

      def self.new(matrix : COOMatrix)
        new(matrix.nrows, matrix.ncolumns, matrix.raw_rows.dup, matrix.raw_columns.dup, matrix.raw_values.map { |v| T.new(v) }, dont_clone: true, dictionary: matrix.dictionary)
      end

      def self.new(matrix : LA::Matrix)
        result = new(matrix.nrows, matrix.ncolumns)
        matrix.each_with_index(all: false) do |v, i, j|
          result.unsafe_set(i, j, v)
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

      def unsafe_set(i, j, value)
        clear_flags # TODO - not always?
        if index = @dictionary[{i, j}]?
          @raw_values[index] = T.new(value)
        else
          # TODO - delete elements?
          @raw_rows << i
          @raw_columns << i
          @raw_values << T.new(value)
          @dictionary[{i, j}] = @raw_rows.size - 1
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
