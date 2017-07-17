require "./general_matrix"
require "./matrix"

module LA::Matrix(T)
  # to_custom(io, "[", ",", "],[", "]")
  def to_custom(io, prefix, columns_separator, rows_separator, postfix)
    each_with_index do |v, r, c|
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

  def to_matlab(io)
    to_custom(io, "[", ", ", "; ", "]")
  end

  def to_matlab
    String.build do |io|
      to_matlab io
    end
  end

  def save_csv(filename)
    File.open(filename, "w") do |f|
      to_custom(f, "", ",", "\n", "\n")
    end
  end
end

# LA::Mat.rand(30, 30).save_csv("./test.csv")
