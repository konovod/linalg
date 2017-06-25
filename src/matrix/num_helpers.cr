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

  def *(m : Matrix(Complex))
    m*self
  end
end

abstract struct Number
  def *(m : Matrix)
    m*self
  end

  def conj
    self
  end
end
