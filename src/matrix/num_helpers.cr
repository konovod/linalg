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

  def +(m : LA::Matrix)
    m + self
  end

  def -(m : LA::Matrix)
    (-m) + self
  end

  def sinh
    (self.exp - (-self).exp) / 2
  end

  def cosh
    (self.exp + (-self).exp) / 2
  end

  def sqrt
    m = self.abs
    s = imag.sign
    Complex.new(Math.sqrt((m + real)/2), s*Math.sqrt((m - real)/2))
  end
end

abstract struct Number
  def *(m : LA::Matrix)
    m*self
  end

  def +(m : LA::Matrix)
    m + self
  end

  def -(m : LA::Matrix)
    (-m) + self
  end

  def conj
    self
  end
end

module Math
  {% for op in %i(exp cos sin tan cosh sinh sqrt) %}
      def {{op.id}}(x : Complex)
        x.{{op.id}}
      end
      def self.{{op.id}}(x : Complex)
        x.{{op.id}}
      end
  {% end %}

  def exp(x : Complex)
    x.exp
  end

  def exp(x : Complex)
    x.exp
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
