module LA
  {% for op in %i(expm cosm sinm tanm coshm sinhm tanhm) %}
     # See `Matrix#{{op.id}}`
     def self.{{op.id}}(mat)
       mat.{{op.id}}
     end
     # See `Matrix#{{op.id}}`
     def {{op.id}}(mat)
       mat.{{op.id}}
     end
  {% end %}

  abstract class Matrix(T)
    # Matrix exponential `exp(A)`
    # See `GeneralMatrix#expm` for details
    def expm
      self.to_general.expm
    end

    # Matrix cosine `cos(A)`, computed via complex exponential.
    def cosm
      # optimization idea for noncomplex matrix is from scipy
      {% if T == Complex %}
        0.5*((1.i*self).expm + (-1.i*self).expm)
      {% else %}
        GMatComplex.new(1.i*self).expm.to_real
      {% end %}
    end

    # Matrix sine `sin(A)`, computed via complex exponential.
    def sinm
      {% if T == Complex %}
        -0.5.i*((1.i*self).expm - (-1.i*self).expm)
      {% else %}
        GMatComplex.new(1.i*self).expm.to_imag
      {% end %}
    end

    # Matrix tangent `tan(A) = sinm(A) * cosm(A)⁻¹`
    def tanm
      self.sinm * self.cosm.inv!
    end

    # Matrix hyperbolic cosine `cosh(A) = (exp(A) + exp(-A))/2`
    def coshm
      0.5*(self.expm + (-self).expm)
    end

    # Matrix hyperbolic sine `sinh(A) = (exp(A) - exp(-A))/2`
    def sinhm
      0.5*(self.expm - (-self).expm)
    end

    # Matrix hyperbolic tangent `tanh(A) = sinhm(A) * coshm(A)⁻¹`
    def tanhm
      self.sinhm * self.coshm.inv!
    end
  end
end
