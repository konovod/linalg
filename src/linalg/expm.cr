module Linalg
  # this is a direct conversion to Crystal of matlab code

  # Coefficients of leading terms in the backward error functions h_{2m+1}.
  private Coeff = [1.0/100800, 1.0/10059033600, 1.0/4487938430976000,
                   1.0/5914384781877411840000.0, 1.0/113250775606021113483283660800000000.0]

  private M_VALS = {3, 5, 7, 9, 13}
  # theta_m for m=1:13.
  private THETA = [
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
  ] # m_vals = 13

  def expm_new(a, schur_fact = false)
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
    raise ArgumentError.new("Matrix must be square for expm") unless a.square?

    schur_fact = false if a.flags.triangular?
    if schur_fact
      a, q = a.schur
    end
    #
    #
    n = a.nrows
    have_A4 = false # prior evaluation of A4
    have_A6 = false # prior evaluation of A6

    s = 0

    a2 = a*a
    eta1 = {normAm(a2, 2)**(1.0/4), normAm(a2, 3)**(1.0/6)}.max
    t = eval_alpha(a, 1)
    if eta1 <= THETA[0] && t == 0
      f = padeApproximantOfDegree(M_VALS[0])
      return schur_fact ? q*f*q.conjt : f
    end
    #
    a4 = a2*a2; have_a4 = 1
    eta2 = {normAm(a4, 1)**(1.0/4), normAm(a2, 3)**(1.0/6)}.max
    t = eval_alpha(a, 2)
    if eta2 <= THETA[1] && t == 0
      f = padeApproximantOfDegree(M_VALS[1])
      return schur_fact ? q*f*q.conjt : f
    end
    #
    a6 = a2*a4; have_a6 = 1
    eta3 = {normAm(a6, 1)**(1.0/6), normAm(a4, 2)**(1.0/8)}.max
    h = [0.0, 0.0, 0.0, 0.0]
    #
    (2..3).each do |i|
      if eta3 <= THETA[i]
        h = eval_alpha(a, i + 1)
        if h == 0
          f = padeApproximantOfDegree(M_VALS[i])
          return schur_fact ? q*f*q.conjt : f
        end
      end
    end
    #
    eta4 = {normAm(a4, 2)**(1.0/8), normAm(a2, 5)**(1.0/10)}.max
    eta5 = {eta3, eta4}.min
    s = {ceil(Math.log2(eta5/THETA.last)), 0}.max # Zero must be here
    t = eval_alpha(A/2 ^ s, 5)
    s = s + t
    a = a/2 ^ s; a2 = a2/2 ^ (2*s); a4 = a4/2 ^ (4*s); a6 = a6/2 ^ (6*s) # Scaling
    #
    f = padeApproximantOfDegree(M_VALS.last)
    if a.flags.lower_triangular?
      a = a.t
      f.t!
    end
    #
    if a.flags.triangular? || schur_fact
      f = expm_sqtri(a, f, s)
      f.transpose! if a.flags.lower_triangular?
    else
      s.times { f = f*f }
    end
    return schur_fact ? q*f*q.conjt : f
  end

  private def eval_alpha(a, k)
    u = eps/2
    alpha = Coeff[k - 1]*normAm(a.map(&.abs), 2*M_VALS[k - 1] + 1)/a.norm(MatrixNorm::One)
    t = {ceil(Math.log2(alpha/u)/(2*M_VALS[k - 1])), 0}.max
  end

  private def padeApproximantOfDegree(m)
    # PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
    #   F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
    #   Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
    #   Series are evaluated in decreasing order of powers, which is
    #   in approx. increasing order of maximum norms of the terms.

    c = getPadeCoefficients

    # Evaluate Pade approximant.
    case m
    when 3, 5, 7, 9
      apowers = Array(typeof(a)).new(ceil((m + 1)/2))
      apowers << eye(n)
      apowers << a2
      if have_A4
        apowers << a4
      end
      if have_A6
        apowers << a6
      end
      apowers.size...ceil((m + 1)/2).each do |j|
        apowers << apowers.last*apowers[1]
      end
      #
      #         U = zeros(n); V = zeros(n);
      #
      #         for j = m+1:-2:2
      #             U = U + c(j)*Apowers{j/2};
      #         end
      #         U = A*U;
      #         for j = m:-2:1
      #             V = V + c(j)*Apowers{(j+1)/2};
      #         end
      #
      #     case 13
      #
      #         % For optimal evaluation need different formula for m >= 12.
      #         U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) ...
      #                  + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n) );
      #
      #         V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) ...
      #             + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n);
      #
    end
    # F = (-U+V)\(U+V);

  end

  def getPadeCoefficients
    # GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
    #    C = GETPADECOEFFICIENTS returns coefficients of numerator
    #    of [M/M] Pade approximant, where M = 3,5,7,9,13.
    # switch m
    #     case 3
    #         c = [120, 60, 12, 1];
    #     case 5
    #         c = [30240, 15120, 3360, 420, 30, 1];
    #     case 7
    #         c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
    #     case 9
    #         c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
    #              2162160, 110880, 3960, 90, 1];
    #     case 13
    #         c = [64764752532480000, 32382376266240000, 7771770303897600, ...
    #              1187353796428800,  129060195264000,   10559470521600, ...
    #              670442572800,      33522128640,       1323241920,...
    #              40840800,          960960,            16380,  182,  1];
    # end
  end
end
