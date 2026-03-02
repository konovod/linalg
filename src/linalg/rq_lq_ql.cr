# Implements RQ, LQ, and QL decompositions using macros.

private macro decomposition(letters, trix)
  def {{letters}}(*, overwrite_a = false)
    a = overwrite_a ? self : clone
    tau = {{letters}}_initial(a)
    r = a.clone
    r.{{trix}}!
    m = a.nrows
    n = a.ncolumns
    k = {m,n}.min
    lapack(org{{letters}}, m, n, k, a, a.nrows, tau)
    a.assume! MatrixFlags::Orthogonal
    {% if letters.id.chars[0] == 'q' %}
      {a, r}
    {% else %}
      {r, a}
    {% end %}
  end

  # :nodoc:
  private def {{letters}}_initial(a)
    m = a.nrows
    n = a.ncolumns
    tau = GeneralMatrix(T).new(1, n)
    lapack(ge{{letters}}f, m, n, a, a.nrows, tau)
    a.clear_flags
    tau
  end

  # :nodoc:
  def {{letters}}_r(*, overwrite_a = false)
    a = overwrite_a ? self : clone
    tau = {{letters}}_initial(a)
    a.{{trix}}!
    a
  end

end

private macro decomposition_wrap(letters)
  # See `GeneralMatrix#{{letters}}`
  def {{letters}}
    to_general.{{letters}}(overwrite_a: true)
  end
  # See `GeneralMatrix#{{letters}}_r`
  def {{letters}}_r
    to_general.{{letters}}_r(overwrite_a: true)
  end
end

module LA
  class Matrix(T)
    decomposition_wrap(rq)
    decomposition_wrap(lq)
    decomposition_wrap(ql)
  end

  class GeneralMatrix(T) < Matrix(T)
    # Computes RQ decomposition.
    #
    # Arguments:
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), GeneralMatrix(T)) : Q and R factors.
    #
    # LAPACK Routines:
    #   - `gerqf` (RQ factorization)
    #   - `orgrq` (generate Q from elementary reflectors)
    decomposition(rq, triu)

    # Computes LQ decomposition.
    #
    # Arguments:
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), GeneralMatrix(T)) : L and Q factors.
    #
    # LAPACK Routines:
    #   - `gelqf` (LQ factorization)
    #   - `orglq` (generate Q from elementary reflectors)
    decomposition(lq, tril)

    # Computes QL decomposition.
    #
    # Arguments:
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), GeneralMatrix(T)) : Q and L factors.
    #
    # LAPACK Routines:
    #   - `geqlf` (QL factorization)
    #   - `orgql` (generate Q from elementary reflectors)
    decomposition(ql, tril)
  end
end
