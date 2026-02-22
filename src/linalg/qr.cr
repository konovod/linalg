# TODO - inline docs

module LA
  class Matrix(T)
    def qr_raw(*, pivoting = false)
      to_general.qr_raw(overwrite_a: true)
    end

    def qr_r(*, pivoting = false)
      to_general.qr_r(overwrite_a: true)
    end

    def qr(*, pivoting = false)
      to_general.qr(overwrite_a: true)
    end
  end

  class GeneralMatrix(T) < Matrix(T)
    # Internal helper for QR factorization.
    #
    # Performs the LAPACK QR factorization and returns TAU and pivot information.
    #
    # Arguments:
    #   - a (GeneralMatrix(T)) : Matrix to factorize (will be modified).
    #   - pivoting (Bool) : If `true`, uses QR with column pivoting.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), Array(Int32)) : TAU array and pivot indices.
    #
    # LAPACK Routines:
    #   - `geqrf` (QR without pivoting)
    #   - `geqp3` (QR with column pivoting)
    private def qr_initial(a, pivoting)
      m = a.nrows
      n = a.ncolumns
      tau = GeneralMatrix(T).new(1, n)
      if pivoting
        jpvt = Array(Int32).new(n, 0)
        lapack(geqp3, m, n, a, a.nrows, jpvt, tau, worksize: [2*n])
        a.clear_flags
      else
        jpvt = Array(Int32).new(0)
        lapack(geqrf, m, n, a, a.nrows, tau)
        a.clear_flags
      end
      {tau, jpvt}
    end

    def qr_raw(*, overwrite_a = false, pivoting = false)
      a = overwrite_a ? self : clone
      tau, pvt = qr_initial(a, pivoting)
      {a, tau, pvt}
    end

    def qr_r(*, overwrite_a = false, pivoting = false)
      a = overwrite_a ? self : clone
      tau, pvt = qr_initial(a, pivoting)
      a.triu!
      {a, pvt}
    end

    # Computes the complete QR decomposition with orthogonal Q.
    #
    # Arguments:
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #   - pivoting (Bool) : If `true`, performs QR with column pivoting. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), GeneralMatrix(T), Array(Int32)) :
    #       Orthogonal matrix Q (m×m or m×n depending on dimensions)
    #       Upper triangular matrix R (n×n)
    #       Pivot indices (empty if pivoting = false)
    #
    # Properties:
    #   - For an m×n matrix A with m ≥ n: A = Q * R
    #   - For m < n: A = Q * R where Q is m×m and R is m×n
    #   - When pivoting is enabled: A * P = Q * R where P is the permutation matrix
    #     defined by pivot indices
    #
    # Side Effects:
    #   - Marks Q as `Orthogonal` on successful decomposition
    #
    # LAPACK Routines:
    #   - `geqrf`/`geqp3` (QR factorization)
    #   - `orgqr` (generate Q from elementary reflectors)
    def qr(*, overwrite_a = false, pivoting = false)
      a = overwrite_a ? self : clone
      tau, pvt = qr_initial(a, pivoting)
      r = a.clone
      r.triu!
      m = a.nrows
      n = a.ncolumns
      k = {m, n}.min
      lapack(orgqr, m, n, k, a, a.nrows, tau)
      a.assume! MatrixFlags::Orthogonal
      {a, r, pvt}
    end
  end
end
