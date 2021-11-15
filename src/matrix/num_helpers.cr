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
    (Math.exp(self) - Math.exp(-self)) / 2
  end

  def cosh
    (Math.exp(self) + Math.exp(-self)) / 2
  end

  def chop
    self.imag.zero? ? self.real : nil
  end

  def self.multiplicative_identity
    new(1.0, 0.0)
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
  {% for op in %i(cos sin tan cosh sinh) %}
      def {{op.id}}(x : Complex)
        x.{{op.id}}
      end
      def self.{{op.id}}(x : Complex)
        x.{{op.id}}
      end
  {% end %}
end

module Enumerable(T)
  def product(initial : Complex, &block)
    reduce(initial) { |memo, e| memo * (yield e) }
  end

  def product(initial : Complex)
    product initial, &.itself
  end
end
