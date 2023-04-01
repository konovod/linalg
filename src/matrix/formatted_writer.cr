require "./general_matrix"
require "./matrix"

abstract class LA::Matrix(T)
  # to_custom(io, "[", ",", "],[", "]")
  # Converts a matrix to string with custom format.
  # Example:
  # ```
  # a = LA::GMat[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
  # str = String.build do |io|
  #   a.to_custom(io, prefix: "(", columns_separator: ",", rows_separator: "|", postfix: ")")
  # end
  # str # => "(1,2,3|4,5,6|7,8,9)"
  # ```
  def to_custom(io, prefix, columns_separator, rows_separator, postfix)
    each_with_index(all: true) do |v, r, c|
      if c > 0
        io << columns_separator
      elsif r > 0
        io << rows_separator
      else
        io << prefix
      end
      io << v
    end
    io << postfix
  end

  # to_custom(io, "[", ",", "],[", "]")
  # Converts a matrix to string with custom format.
  # Example:
  # ```
  # a = LA::GMat[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
  # str = a.to_custom(prefix: "(", columns_separator: ",", rows_separator: "|", postfix: ")")
  # str # => "(1,2,3|4,5,6|7,8,9)"
  # ```
  def to_custom(prefix, columns_separator, rows_separator, postfix)
    String.build do |io|
      to_custom(io, prefix, columns_separator, rows_separator, postfix)
    end
  end

  # Converts a matrix to matlab format
  def to_matlab(io)
    to_custom(io, "[", ", ", "; ", "]")
  end

  # Converts a matrix to matlab format string
  # Example:
  # ```
  # LA::GMat[[1, 2, 3], [4, 5, 6], [7, 8, 9]].to_matlab # => "[1,2,3; 4,5,6; 7,8,9]"
  # ```
  def to_matlab
    String.build do |io|
      to_matlab io
    end
  end

  # Save a matrix to CSV (comma separated values) file
  # Example:
  # ```
  # LA::GMat.rand(30, 30).save_csv("./test.csv")
  # a = LA::GMat.load_csv("./test.csv")
  # ```
  def save_csv(filename)
    File.open(filename, "w") do |f|
      to_custom(f, "", ",", "\n", "\n")
    end
  end
end
