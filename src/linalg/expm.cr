module LA
  # this is a direct conversion to Crystal of matlab code

  # Coefficients of leading terms in the backward error functions h_{2m+1}.
  private Coeff = {1.0/100800, 1.0/10059033600, 1.0/4487938430976000,
                   1.0/5914384781877411840000.0, 1.0/113250775606021113483283660800000000.0}

  private M_VALS = {3, 5, 7, 9, 13}
  # theta_m for m=1:13.
  private THETA = {
    # 3.650024139523051e-008
    # 5.317232856892575e-004
    1.495585217958292e-002, # m_vals = 3
    # 8.536352760102745e-002
    2.539398330063230e-001, # m_vals = 5
    # 5.414660951208968e-001
    9.504178996162932e-001, # m_vals = 7
    # 1.473163964234804e+000
    2.097847961257068e+000, # m_vals = 9
    # 2.811644121620263e+000
    # 3.602330066265032e+000
    # 4.458935413036850e+000
    4.250000000000000e+000,
  } # m_vals = 13

  module Matrix(T)
    def expm(*, schur_fact = false)
      # EXPM_NEW  Matrix exponential.
      #  EXPM_NEW(A) is the matrix exponential of A computed using
      #  an improved scaling and squaring algorithm with a Pade approximation.
      #  It exploits triangularity (if any) of A.
      #  EXPM_NEW(A,1) uses an initial transformation to complex Schur form.
      #  EXPM_NEW(A,2) uses an initial transformation to real Schur form if A
      #  is real.
      #
      #  Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
      #  Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
      #  970-989, 2009.
      #  Awad H. Al-Mohy and Nicholas J. Higham, April 20, 2010.
      raise ArgumentError.new("Matrix must be square for expm") unless square?
      return Matrix(T).identity(nrows) if norm(MatrixNorm::One) == 0
      if flags.diagonal?
        return Matrix(T).diag(nrows, nrows) do |i|
          v = unsafe_at(i, i)
          {% if T == Complex %} v.exp {% else %} Math.exp(v) {% end %}
        end
      end
      return clone.expm2_by_2! if nrows == 2
      schur_fact = false if flags.triangular?
      if schur_fact
        a, q = self.schur
      else
        a = self
        q = Matrix(T).zeros(0, 0) # to prevent nilability of q
      end
      n = a.nrows
      s = 0

      a2 = a*a
      eta1 = {a2.normAm(2)**(1.0/4), a2.normAm(3)**(1.0/6)}.max
      t = eval_alpha(a, 1)
      if eta1 <= THETA[0] && t == 0
        f = a.padeApproximantOfDegree(M_VALS[0], a2)
        raise "" unless q
        return schur_fact ? q*f*q.conjt : f
      end
      #
      a4 = a2*a2
      eta2 = {a4.normAm(1)**(1.0/4), a2.normAm(3)**(1.0/6)}.max
      t = eval_alpha(a, 2)
      if eta2 <= THETA[1] && t == 0
        f = a.padeApproximantOfDegree(M_VALS[1], a2, a4)
        raise "" unless q
        return schur_fact ? q*f*q.conjt : f
      end
      #
      a6 = a2*a4
      eta3 = {a6.normAm(1)**(1.0/6), a4.normAm(2)**(1.0/8)}.max
      h = [0.0, 0.0, 0.0, 0.0]
      #
      (2..3).each do |i|
        if eta3 <= THETA[i]
          h = eval_alpha(a, i + 1)
          if h == 0
            f = a.padeApproximantOfDegree(M_VALS[i], a2, a4, a6)
            return schur_fact ? q*f*q.conjt : f
          end
        end
      end
      #
      eta4 = {a4.normAm(2)**(1.0/8), a2.normAm(5)**(1.0/10)}.max
      eta5 = {eta3, eta4}.min
      s = {(Math.log2(eta5/THETA.last)).ceil, 0}.max # Zero must be here
      t = eval_alpha(a/2 ** s, 5)
      s = s + t
      a = a/2 ** s; a2 = a2/2 ** (2*s); a4 = a4/2 ** (4*s); a6 = a6/2 ** (6*s) # Scaling
      #
      f = a.padeApproximantOfDegree(M_VALS.last, a2, a4, a6)
      lower = a.flags.lower_triangular? && !a.flags.upper_triangular?
      if lower
        a.t!
        f.t!
      end
      if a.flags.triangular? || schur_fact
        f = expm_sqtri(a, f, s)
        f.transpose! if lower
      else
        s.to_i.times { f = f*f }
      end
      return schur_fact ? q*f*q.conjt : f
    end

    protected def eval_alpha(a, k)
      eps = real_type_const(EPSILON)
      u = eps/2
      alpha = Coeff[k - 1]*a.map(&.abs).normAm(2*M_VALS[k - 1] + 1, overwrite_a: true) / a.norm(MatrixNorm::One)
      t = {(Math.log2(alpha/u)/(2*M_VALS[k - 1])).ceil, 0}.max
    end

    protected def padeApproximantOfDegree(m, a2, a4 = nil, a6 = nil)
      # PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
      #   F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
      #   Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
      #   Series are evaluated in decreasing order of powers, which is
      #   in approx. increasing order of maximum norms of the terms.

      c = getPadeCoefficients(m)
      n = self.nrows

      # Evaluate Pade approximant.
      case m
      when 3, 5, 7, 9
        m2 = (m + 3) / 2 # ((m + 1)/2.0).ceil
        apowers = Array(Matrix(T)).new(m2)
        apowers << Matrix(T).eye(n)
        apowers << a2
        if a4
          apowers << a4
        end
        if a6
          apowers << a6
        end
        (apowers.size...m2).each do |j|
          apowers << apowers.last*a2
        end
        u = Matrix(T).zeros(n, n)
        v = Matrix(T).zeros(n, n)
        #
        (2..m + 1).reverse_each.step(2).each do |j|
          u += c[j - 1] * apowers[j/2 - 1]
        end
        u = self*u
        (1..m).reverse_each.step(2).each do |j|
          v += c[j - 1] * apowers[(j + 1)/2 - 1]
        end
      when 13
        raise "" unless a4
        raise "" unless a6
        # For optimal evaluation need different formula for m >= 12.
        u = self * (a6*(c[14 - 1]*a6 + c[12 - 1]*a4 + c[10 - 1]*a2) + c[8 - 1]*a6 + c[6 - 1]*a4 + c[4 - 1]*a2 + c[2 - 1]*Matrix(T).eye(n))

        v = a6*(c[13 - 1]*a6 + c[11 - 1]*a4 + c[9 - 1]*a2) + c[7 - 1]*a6 + c[5 - 1]*a4 + c[3 - 1]*a2 + c[1 - 1]*Matrix(T).eye(n)
      else
        raise ""
      end
      return (v - u).inv!*(u + v)
    end

    protected def getPadeCoefficients(m)
      # GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
      #    C = GETPADECOEFFICIENTS returns coefficients of numerator
      #    of [M/M] Pade approximant, where M = 3,5,7,9,13.
      case m
      when 3
        {120.0, 60.0, 12.0, 1.0}
      when 5
        {30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0}
      when 7
        {17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0}
      when 9
        {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0,
         2162160.0, 110880.0, 3960.0, 90.0, 1.0}
      when 13
        {64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
         1187353796428800.0, 129060195264000.0, 10559470521600.0,
         670442572800.0, 33522128640.0, 1323241920.0,
         40840800.0, 960960.0, 16380.0, 182.0, 1.0}
      else raise ""
      end
    end

    protected def normAm(m, *, overwrite_a = false)
      # NORMAM   Estimate of 1-norm of power of matrix.
      #    NORMAM(A,m) estimates norm(A^m,1).
      #    If A has nonnegative elements the estimate is exact.
      #    [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
      #    matrix-vector products computed involving A or A^*.
      #    Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
      #    Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
      #    970-989, 2009.

      #    Awad H. Al-Mohy and Nicholas J. Higham, April 19, 2010.
      n = nrows
      a = overwrite_a ? self : clone

      if a.raw.all? { |v| {% if T == Complex %} v.real >= 0{% else %} v >= 0{% end %} }
        e = Matrix(T).ones(n, 1)
        a.transpose!
        m.times { e = a * e }
        return e.norm(MatrixNorm::Inf)
      else
        # TODO - normest?
        # [c,v,w,it] = normest1(@afun_power);
        # mv = it(2)*2*m; % Since t = 2.
        return a.map!(&.abs).normAm(m, overwrite_a: overwrite_a)
      end
    end

    #   function Z = afun_power(flag,X)
    #        %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

    #        if isequal(flag,'dim')
    #           Z = n;
    #        elseif isequal(flag,'real')
    #           Z = isreal(A);
    #        else

    #           [p,q] = size(X);
    #           if p ~= n, error('Dimension mismatch'), end

    #           if isequal(flag,'notransp')
    #              for i = 1:m, X = A*X; end
    #           elseif isequal(flag,'transp')
    #              for i = 1:m, X = A'*X; end
    #           end

    #           Z = X;

    #        end

    #   end
    # end

    def expm_sqtri(t, f, s)
      # EXPM_SQTRI   Squaring phase of scaling and squaring method.
      #   X = EXPM_SQTRI(T/2^s,F,s) carries out the squaring phase
      #   of the scaling and squaring method for an upper quasitriangular T,
      #   given T/2^s and a Pade approximant F to e^{T/2^s}.
      #   It corrects the diagonal blocks blocks at each step.
      #   This M-file exploits Code Fragment 2.1 and Code Fragment 2.2 of the
      #   reference below.
      #   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
      #   Algorithm for the Matrix Exponential,SIAM J. Matrix Anal. Appl. 31(3):
      #   970-989, 2009.
      #   Awad H. Al-Mohy and Nicholas J. Higham, April 19, 2010.
      k = 1
      # To turn off exact superdiagonal computation force "istriangular = 0".
      istriangular = t.flags.upper_triangular?
      raise "not implemented" unless istriangular
      # c = diag(-1).to_a.map { |x| x.abs > 0 ? 1 : 0 } # sum(c) = number of 2by2 full blocks
      #    NUMBLK blocks with i'th block in rows/cols INDX{i}.
      # numblk = nrows - c.sum # The number of blocks
      #    indx = cell(numblk,1);
      #    if c(end) == 0
      #        indx{end} = n; c = [c ; 0];
      #    end
      #    for j = 1:numblk
      #        if c(k)
      #            indx{j} = k:k+1; k = k+2;
      #        else
      #            indx{j} = k; k = k+1;
      #        end
      #    end
      i = 0
      while i <= s + 0.001
        f = f*f if i > 0
        if istriangular
          #         Compute diagonal and first superdiagonal.
          (1..f.nrows).step(2) do |j|
            if j < f.nrows - 1
              f[j - 1..j, j - 1..j] = (2 ** i * t[j - 1..j, j - 1..j]).expmT2by2!
            else
              f[j, j] = Math.exp(2 ** i * t[j, j])
            end
          end
        else
          #         Quasitriangular case: compute (block) diagonal only.
          #        for j = 1:numblk
          #            F(indx{j},indx{j}) = expm2_by_2( 2^i * T(indx{j},indx{j}) );
        end
        i += 1
      end
      return f
    end

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    protected def expm2_by_2!
      # EXPM2_BY_2  Exponential for a general 2-by-2 matrix A.
      case nrows
      when 1
        unsafe_set(0, 0, Math.exp(unsafe_at(0, 0)))
      when 2
        flags.upper_triangular? ? expmT2by2! : expm2by2full!
      else
        raise "BUG"
      end
      self
    end

    private def sinch(x)
      x == 0 ? T.new(1.0) : Math.sinh(x)/x
    end

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    private def expm2by2full!
      # EXPM2BY2FULL   Exponential of 2-by-2 full matrix.
      a = unsafe_at(0, 0)
      b = unsafe_at(0, 1)
      c = unsafe_at(1, 0)
      d = unsafe_at(1, 1)

      delta = Math.sqrt((a - d)*(a - d) + 4*b*c)
      k = Math.exp((a + d)/2)
      unsafe_set(0, 0, k*(Math.cosh(delta/2) + (a - d)/2*sinch(delta/2)))
      unsafe_set(0, 1, k*b*sinch(delta/2))
      unsafe_set(1, 0, k*c*sinch(delta/2))
      unsafe_set(1, 1, k*(Math.cosh(delta/2) + (d - a)/2*sinch(delta/2)))
      self
    end

    # %%%%%%%%%%%%%%%%%%%%%%%%
    protected def expmT2by2!
      # EXPMT2BY2    Exponential of 2-by-2 upper triangular matrix.
      # EXPMT2BY2(A) is the exponential of the 2-by-2 upper triangular matrix A.

      # Modified from FUNM (EXPM2by2).

      a1 = unsafe_at(0, 0)
      a2 = unsafe_at(1, 1)
      c = unsafe_at(0, 1)
      ave = (a1 + a2)/2
      {% if T == Complex %}ave = ave.real{% end %}
      df = (a1 - a2).abs/2
      if Math.max(ave, df) < Math.log(real_type_const(MAX_FINITE))
        # Formula fine unless it overflows.
        x12 = c*Math.exp((a1 + a2)/2) * sinch((a2 - a1)/2)
      else
        # Formula can suffer cancellation.
        x12 = c*(Math.exp(a2) - Math.exp(a1))/(a2 - a1)
      end

      unsafe_set(0, 0, Math.exp(a1))
      unsafe_set(0, 1, x12)
      unsafe_set(1, 0, 0.0)
      unsafe_set(1, 1, Math.exp(a2))
      self
    end
  end
end
