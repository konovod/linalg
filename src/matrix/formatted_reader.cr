require "./general_matrix"
require "./matrix"

module LA::Utils
  # :nodoc:
  # reads io until one of the delimiters found
  # returns nil if none of the delimiters found
  # otherwise returns tuple consisting of delimiter index and substring
  def self.prepare_multi_gets(delimiters)
    scanner = delimiters.map_with_index do |x, i|
      {size: x.bytesize, index: i, data: x, last: x.byte_at(x.bytesize - 1)}
    end.to_a
    scanner.sort_by { |x| -x[:size] }
    scanner
  end

  # :nodoc:
  def self.multi_gets(io, scanner) : {Int32, String}?
    # The 'hard' case: we read until we match the last byte,
    # and then compare backwards
    total_bytes = 0

    buffer = String::Builder.new
    while true
      return nil unless byte = io.read_utf8_byte
      buffer.write_byte(byte)
      total_bytes += 1
      scanner.each do |s|
        if s[:last] == byte &&
           (buffer.bytesize >= s[:size]) &&
           (buffer.buffer + total_bytes - s[:size]).memcmp(s[:data].to_unsafe, s[:size]) == 0
          buffer.back(s[:size])
          return {s[:index], buffer.to_s}
        end
      end
    end
    return nil
  end
end

class LA::Utils::MatrixParseError < Exception
end

abstract class LA::Matrix(T)
  # Converts a string of given format to matrix
  # Example:
  # ```
  # str = "( 1,2,3 | 4,5,6 | 7,8,9 )"
  # LA::GMat.from_custom(IO::Memory.new(str), prefix: "(", columns_separator: ",", rows_separator: "|", postfix: ")")
  # ```
  def self.from_custom(io, prefix, columns_separator, rows_separator, postfix)
    io.gets(prefix)
    postfix_empty = postfix == ""
    arr = postfix_empty ? {columns_separator, rows_separator} : {columns_separator, rows_separator, postfix}
    scan = Utils.prepare_multi_gets(arr)
    row = 0
    column = 0
    last_column = 0
    data = Array(T).new
    loop do
      token = Utils.multi_gets(io, scan)
      if token
        data << to_my_type(token[1])
        case token[0]
        when 0 # columns_separator
          raise MatrixParseError.new("row #{row} too long") if column >= last_column && row > 0
          column += 1
        when 1 # rows_separator
          if row == 0
            last_column = column
          else
            raise MatrixParseError.new("row #{row} too short") if column != last_column
          end
          row += 1
          column = 0
        else # postfix
          break
        end
      else
        if postfix_empty
          row -= 1
          break # go to matrix construction
        else
          raise MatrixParseError.new("postfix not found")
        end
      end
    end
    GeneralMatrix(T).new(row + 1, last_column + 1, data)
  end

  # Converts a string of given format to matrix
  # Example:
  # ```
  # str = "( 1,2,3 | 4,5,6 | 7,8,9 )"
  # LA::GMat.from_custom(str, "(", ",", "|", ")")
  # ```
  def self.from_custom(str : String, prefix, columns_separator, rows_separator, postfix)
    from_custom(IO::Memory.new(str), prefix, columns_separator, rows_separator, postfix)
  end

  # Converts a string from matlab format to matrix
  # Example:
  # ```
  # str = "[1,2,3; 4,5,6; 7,8,9]"
  # LA::GMat.from_matlab(str)
  # ```
  def self.from_matlab(s)
    from_custom(IO::Memory.new(s), "[", ",", ";", "]")
  end

  # Loads a matrix from CSV (comma separated values) file
  # Example:
  # ```
  # LA::GMat.rand(30, 30).save_csv("./test.csv")
  # a = LA::GMat.load_csv("./test.csv")
  # ```
  def self.load_csv(filename)
    File.open(filename, "r") do |f|
      return from_custom(f, "", ",", "\n", "")
    end
  end
end
