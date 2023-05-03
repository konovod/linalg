# TODO - inline docs

module LA
  class BandedMatrix(T) < Matrix(T)
    # returns matrix norm
    def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
      # TODO - check if not square
      let = case kind
            in .frobenius?
              'F'
            in .one?
              'o'
            in .inf?
              'I'
            in .max_abs?
              'M'
            end.ord.to_u8

      worksize = kind.inf? ? nrows : 0
      lapack_util(langb, worksize, let, @nrows, @lower_band, @upper_band, matrix(self), @lower_band + @upper_band + 1)
    end

    def det(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      if flags.triangular?
        return diag.product
      end
      lru = overwrite_a ? self : self.clone
      lru.upper_band = @lower_band + @upper_band
      ipiv = Slice(Int32).new(nrows)
      lapack(gbtrf, nrows, nrows, @lower_band, @upper_band, lru, 2*@lower_band + @upper_band + 1, ipiv)
      lru.clear_flags
      lru.diag.product
    end

    def solve(b : GeneralMatrix(T), *, overwrite_a = false, overwrite_b = false)
      raise ArgumentError.new("nrows of a and b must match") unless nrows == b.nrows
      raise ArgumentError.new("a must be square") unless square?
      a = overwrite_b ? self : self.clone
      x = overwrite_b ? b : b.clone
      n = nrows
      ku = upper_band
      kl = lower_band

      if flags.triangular?
        kd = flags.lower_triangular? ? kl : ku
        if flags.lower_triangular?
          self.upper_band = 0
        else
          self.lower_band = 0
        end
        lapack(tbtrs, uplo, 'N'.ord.to_u8, 'N'.ord.to_u8, n, kd, b.ncolumns, a, kd + 1, x, n)
      elsif flags.positive_definite?
        lapack(pbsv, 'U'.ord.to_u8, n, upper_band, b.ncolumns, a, upper_band + lower_band + 1, x, n)
      else
        a.upper_band = kl + ku
        ipiv = Slice(Int32).new(n)
        lapack(gbsv, n, kl, ku, b.ncolumns, a, 2*kl + ku + 1, ipiv, x, b.nrows)
      end

      a.clear_flags
      x.clear_flags
      x
    end
  end
end
