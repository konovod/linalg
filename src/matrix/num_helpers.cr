require "complex"

struct Complex
  def self.new(value)
    case value
    when Complex
      new(value.real, value.imag)
    else
      new(value, 0.0)
    end
  end
end
