module LA
  abstract class Matrix(T)
    def lu(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = nrows
      n = ncolumns
      k = {nrows, ncolumns}.min
      ipiv = Slice(Int32).new(m)
      lapack(getrf, nrows, ncolumns, matrix(a), nrows, matrix(ipiv))
      a.clear_flags
      # apply all transformation of piv to "own" piv
      piv = (1..m).to_a
      k.times do |i|
        tmp = piv[i]
        piv[i] = piv[ipiv[i] - 1]
        piv[ipiv[i] - 1] = tmp
      end
      p = GeneralMatrix(T).new(m, m)
      m.times do |i|
        p[piv[i] - 1, i] = T.new(1)
      end
      l = GeneralMatrix(T).new(m, k) do |i, j|
        case i <=> j
        when 0
          T.new(1)
        when 1
          a[i, j]
        else
          T.new(0)
        end
      end
      u = GeneralMatrix(T).new(k, n) do |i, j|
        if i <= j
          a[i, j]
        else
          T.new(0)
        end
      end
      {p, l, u}
    end

    def lu_factor!
      raise ArgumentError.new("matrix must be square") unless square?
      ipiv = Slice(Int32).new(nrows)
      lapack(getrf, nrows, ncolumns, matrix(self), nrows, matrix(ipiv))
      clear_flags
      LUMatrix(T).new(self, ipiv)
    end

    def lu_factor
      clone.lu_factor!
    end
  end

  enum LUTranspose
    None          = 0
    Transpose
    ConjTranspose
  end

  struct LUMatrix(T)
    @a : Matrix(T)
    @ipiv : Slice(Int32)

    # TODO - more macro magic?
    macro lapack(name, *args, **worksizes)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
       {% for arg, index in args %}
         {% if !(arg.stringify =~ /^matrix\(.*\)$/) %}
           %var{index} = {{arg}}
         {% end %}
       {% end %}
       info = 0
       LibLAPACK.{{typ}}{{name}}_(
         {% for arg, index in args %}
           {% if !(arg.stringify =~ /^matrix\(.*\)$/) %}
             pointerof(%var{index}),
           {% else %}
             {{arg.stringify.gsub(/^matrix\((.*)\)$/, "(\\1)").id}},
           {% end %}
         {% end %}
         pointerof(info))
      raise LinAlgError.new("LAPACK.{{typ}}{{name}} returned #{info}") if info != 0
    end

    def initialize(@a, @ipiv)
    end

    def size
      @a.nrows
    end

    def solve(b, transpose = LUTranspose::None, *, overwrite_b = false)
      raise ArgumentError.new("nrows of a and b must match") unless @a.nrows == b.nrows
      trans = case transpose
              when .none?           then 'N'
              when .transpose?      then 'T'
              when .conj_transpose? then 'C'
              else                       'N'
              end.ord.to_u8
      x = overwrite_b ? b : b.clone
      lapack(getrs, trans, size, b.ncolumns, matrix(@a), size, matrix(@ipiv), matrix(x), x.nrows)
      x.clear_flags
      x
    end
  end
end
