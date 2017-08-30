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

struct Float32
  # Smallest finite value
  MIN_FINITE = 0xff7fffff_u32.unsafe_as(Float32) # -3.40282347e+38_f32
  # Largest finite value
  MAX_FINITE = 0x7f7fffff_u32.unsafe_as(Float32) # 3.40282347e+38_f32
end

struct Float64
  # Smallest finite value
  MIN_FINITE = 0xffefffffffffffff_u64.unsafe_as(Float64) # -1.7976931348623157e+308_f64
  # Largest finite value
  MAX_FINITE = 0x7fefffffffffffff_u64.unsafe_as(Float64) # 1.7976931348623157e+308_f64
end
