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

  def *(m : LA::Matrix(Complex))
    m*self
  end
end

abstract struct Number
  def *(m : LA::Matrix)
    m*self
  end

  def conj
    self
  end
end

module Enumerable(T)
  def product(initial : Complex, &block)
    reduce(initial) { |memo, e| memo * (yield e) }
  end

  def product(initial : Complex)
    product initial, &.itself
  end
end
