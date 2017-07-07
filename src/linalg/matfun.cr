module Linalg
  def self.expm(mat)
    mat.expm
  end

  def self.sinm(mat)
    mat.sinm
  end

  def self.cosm(mat)
    mat.cosm
  end

  module Matrix(T)
    # optimization idea for noncomplex matrix is from scipy
    def cosm
      {% if T == Complex %}
        0.5*(expm(1.i*self)-expm(-1.i*self))
      {% else %}
        expm(GMatComplex.new(1.i*self)).to_real
      {% end %}
    end

    def sinm
      {% if T == Complex %}
        -0.5*(expm(1.i*self)-expm(-1.i*self))
      {% else %}
        expm(GMatComplex.new(1.i*self)).to_imag
      {% end %}
    end
  end
end
