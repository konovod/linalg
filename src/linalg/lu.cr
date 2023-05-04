module LA
  class GeneralMatrix(T) < Matrix(T)
    # Compute pivoted LU decomposition of a matrix
    #
    # If you call `p,l,u = a.lu` for a matrix m*n `a` then
    #
    #  - `p` will be m*m permutation matrix
    #  - `l` will be m*k lower triangular matrix with unit diagonal. K = min(M, N)
    #  - `u` will be k*n upper triangular matrix
    #  - `(p*l*u)` will be equal to `a` (within calculation tolerance)
    #
    # if overwrite_a is true, a will be overriden in process of calculation
    #
    # Uses `getrf` LAPACK routine
    def lu(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = nrows
      n = ncolumns
      k = {nrows, ncolumns}.min
      ipiv = Slice(Int32).new(m)
      lapack(getrf, nrows, ncolumns, a, nrows, ipiv)
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

    # Compute pivoted LU decomposition of a matrix and returns it in a compact form,
    # useful for solving linear equation systems.
    #
    # Overrides source matrix in a process of calculation
    #
    # Uses `getrf` LAPACK routine
    def lu_factor! : LUMatrix(T)
      raise ArgumentError.new("matrix must be square") unless square?
      ipiv = Slice(Int32).new(nrows)
      lapack(getrf, nrows, ncolumns, self, nrows, ipiv)
      clear_flags
      LUMatrix(T).new(self, ipiv)
    end

    # Compute pivoted LU decomposition of a matrix and returns it in a compact form,
    # useful for solving linear equation systems.
    #
    # Uses `getrf` LAPACK routine
    def lu_factor : LUMatrix(T)
      clone.lu_factor!
    end
  end

  enum Enums::LUTranspose
    None          = 0
    Transpose
    ConjTranspose
  end

  # Struct holding lu-decomposition of a matrix
  #
  # Can be used to solve linear equations
  struct Utils::LUMatrix(T)
    # L and U matrices packed in one matrix
    getter a : Matrix(T)
    # indices of permutations
    getter ipiv : Slice(Int32)

    private macro lapack(name, *args)
    # TODO - more macro magic?

      {%
        lapack_funcs = {
          "getrs" => {4 => LapackHelper::ARG_MATRIX, 6 => LapackHelper::ARG_MATRIX, 7 => LapackHelper::ARG_MATRIX},
        }
      %}
      {% func_data = lapack_funcs[name.stringify] %}

      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
         {% for arg, index in args %}
           {% argtype = func_data[index + 1] %}
           {% if argtype == LapackHelper::ARG_MATRIX %}
           {% else %}
           %var{index} = {{arg}}
           {% end %}
         {% end %}
       %info = 0
       LibLAPACK.{{typ}}{{name}}(
       {% for arg, index in args %}
       {% argtype = func_data[index + 1] %}
       {% if argtype == LapackHelper::ARG_MATRIX %}
         {{arg}},
       {% else %}
        pointerof(%var{index}),
       {% end %}
       {% end %}
         pointerof(%info))
      raise LinAlgError.new("LAPACK.{{typ}}{{name}} returned #{%info}") if %info != 0
    end

    def initialize(@a, @ipiv)
    end

    # Returns size of a system
    def size
      @a.nrows
    end

    # Solves `A*x = b` equations system with given b and returns `x`
    #
    # `transpose` allows to solve `A.t * x = b` and `A.conjt * x = b` systems instead
    #
    # if overwrite_b is true, b will be overriden in process of calculation
    #
    # Uses `getrs` LAPACK routine
    def solve(b, transpose = LUTranspose::None, *, overwrite_b = false)
      raise ArgumentError.new("nrows of a and b must match") unless @a.nrows == b.nrows
      trans = case transpose
              in .none?           then 'N'
              in .transpose?      then 'T'
              in .conj_transpose? then 'C'
              end.ord.to_u8
      x = overwrite_b ? b : b.clone
      lapack(getrs, trans, size, b.ncolumns, @a, size, @ipiv, x, x.nrows)
      x.clear_flags
      x
    end
  end
end
