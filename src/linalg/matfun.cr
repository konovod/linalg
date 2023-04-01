# TODO - inline docs

module LA
  {% for op in %i(expm cosm sinm tanm coshm sinhm tanhm) %}
     def self.{{op.id}}(mat)
       mat.{{op.id}}
     end
     def {{op.id}}(mat)
       mat.{{op.id}}
     end
  {% end %}

  abstract class Matrix(T)
    # optimization idea for noncomplex matrix is from scipy
    def cosm
      {% if T == Complex %}
        0.5*((1.i*self).expm + (-1.i*self).expm)
      {% else %}
        GMatComplex.new(1.i*self).expm.to_real
      {% end %}
    end

    def sinm
      {% if T == Complex %}
        -0.5.i*((1.i*self).expm - (-1.i*self).expm)
      {% else %}
        GMatComplex.new(1.i*self).expm.to_imag
      {% end %}
    end

    def tanm
      self.sinm * self.cosm.inv!
    end

    def coshm
      0.5*(self.expm + (-self).expm)
    end

    def sinhm
      0.5*(self.expm - (-self).expm)
    end

    def tanhm
      self.sinhm * self.coshm.inv!
    end
  end
end
