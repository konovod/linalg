require "complex"

struct Complex
  # Multiply scalar to matrix
  def *(m : LA::Matrix(Complex))
    m*self
  end

  # Adds scalar to matrix
  def +(m : LA::Matrix)
    m + self
  end

  # Substract matrix from scalar
  def -(m : LA::Matrix)
    (-m) + self
  end

  # Sinh for complex numbers
  def sinh
    (Math.exp(self) - Math.exp(-self)) / 2
  end

  # Cosh for complex numbers
  def cosh
    (Math.exp(self) + Math.exp(-self)) / 2
  end

  # Convert to real if imaginary paret is zero
  #
  # Returns nil if part is not zero
  def chop
    self.imag.zero? ? self.real : nil
  end
end

abstract struct Number
  # Multiply scalar to matrix
  def *(m : LA::Matrix)
    m*self
  end

  # Adds scalar to matrix
  def +(m : LA::Matrix)
    m + self
  end

  # Substract matrix from scalar
  def -(m : LA::Matrix)
    (-m) + self
  end

  # Returns self
  #
  # Because complex conjurgate for real numbers is number itself.
  # This method is useful to streamline work with complex and real numbers
  # def conj
  #   self
  # end
end

module Math
  {% for op in %i(cos sin tan cosh sinh) %}
      # :nodoc:
      def {{op.id}}(x : Complex)
        x.{{op.id}}
      end
      # {{op.id}} for complex numbers
      def self.{{op.id}}(x : Complex)
        x.{{op.id}}
      end
  {% end %}
end

module Enumerable(T)
  # Same as `#product` from stdlib, but for complex numbers
  def product(initial : Complex, &block)
    reduce(initial) { |memo, e| memo * (yield e) }
  end

  # Same as `#product` from stdlib, but for complex numbers
  def product(initial : Complex)
    product initial, &.itself
  end
end
