abstract class LA::Matrix(T)
  # Taken from https://github.com/Exilor/matrix/
  def **(other : Int)
    raise ArgumentError.new("matrix must be square") unless square?
    m = self
    if other == 0
      self.class.identity(nrows)
    elsif other < 0
      repeated_square_power(other.abs).inv!
    elsif other == 1
      clone
    else
      repeated_square_power(other)
    end
  end

  private def repeated_square_power(n)
    result = self.class.identity(nrows)
    square = self

    while n > 0
      result *= square if (n & 1) == 1
      square *= square
      n >>= 1
    end

    result
  end

  #
  # function [X,nsq,m] = powerm_pade(A,p)
  # %POWERM_PADE The Schur-Pade algorithm for an arbitrary matrix power.
  # %   X = POWERM_PADE(A,P) computes the P'th power X of the matrix A,
  # %   for arbitrary real P and A with no nonpositive real eigenvalues,
  # %   by the Schur-Pade algorithm.
  # %   [X,NSQ,M] = POWERM_PADE(A, P) returns the number NSQ of matrix
  # %   square roots computed and the degree M of the Pade approximant used.
  # %   If A is singular or has any eigenvalues on the negative real axis,
  # %   a warning message is printed.
  #
  # %   This code is intended for double precision.
  #
  # %   See also LOGM, EXPM, FUNM.
  #
  # %   Reference: N. J. Higham and L. Lin, A Schur--Pad\'e Algorithm for
  # %   Fractional Powers of a Matrix.
  # %   MIMS EPrint 2010.91, The University of Manchester, October 2010,
  # %   revised February 2011.
  # %   Name of corresponding algorithm in that paper: Algorithm 5.1/SPade.
  #
  # %   Nicholas J. Higham and Lijing Lin, February 22, 2011.
  #
  # if ~isfloat(A) || ndims(A) ~= 2 || diff(size(A))
  #    error('MATLAB:powerm_pade:InputDim',...
  #          'First input must be a floating point, square matrix.');
  # end
  #
  # maxsqrt = 64;
  # n = length(A);
  # nsq = 0; m = 0;
  #
  # if p == 0, X = eye(n); return, end
  # if n == 1, X = A^p; return, end
  #
  # pint = floor(p);
  # pfrac = p - pint;       % pfrac >= 0
  # if pfrac == 0
  #    if p > 0             % positive integer
  #       X = mpower(A,p);  % Built-in MATLAB function.  Same as A^p.
  #    elseif p < 0         % negative integer
  #       X = inv(A)^(-p);
  #    end
  #    return
  # end
  # if pint == -1           % -1 < p < 0
  #     pfrac = p; pint = 0;
  # end
  #
  # % First form complex Schur form (if A not already upper triangular).
  # if isequal(A,triu(A))
  #    T = A; Q = eye(n);
  # else
  #    [Q,T] = schur(A,'complex');
  # end
  # diagT = diag(T);
  #
  # if ~all(diagT)
  #     warning('MATLAB:powerm_pade:zeroEig', ...
  #             'Matrix power may not exist for singular matrix.')
  # end
  # if any( imag(diagT) == 0 & real(diagT) <= 0 )
  #     warning('MATLAB:powerm_pade:nonPosRealEig', ...
  #            ['Principal matrix power is not defined for A with\n', ...
  #             '         nonpositive real eigenvalues. A non-principal matrix\n', ...
  #             '         power is returned.'])
  # end
  #
  # if isequal(T,diag(diagT)) % Handle special case of diagonal T.
  #    X = Q * diag(diagT.^p) * Q';
  #    return
  # end
  #
  # [X,nsq,m] = powerm_triang(T,pfrac,maxsqrt);
  #
  # % Experimentally, more accurate results are obtained if we transform
  # % back now and use original A for integer power rather than compute
  # % integer power with T and transform back at the end.
  #
  # X = Q*X*Q';
  # if pint > 0
  #     X = mpower(A,pint) * X;
  # elseif pint < 0
  #     X = inv(A)^(-pint) * X;
  # end
  #
  # if isreal(A) && norm(imag(X),1) <= 10*n*eps*norm(X,1)
  #    X = real(X);
  # end
  # end
  #
  # function [X,nsq,m] = powerm_triang(T,p,maxsqrt)
  # %POWERM_TRIANG   Power of triangular matrix by Pade-based algorithm.
  # %   X = POWERM_TRIANG(T,P,MAXSQRT) computes the P'th power X of
  # %   the upper triangular matrix T, for an arbitrary real number P in the
  # %   interval (-1,1),and T with no nonpositive real eigenvalues, by a
  # %   Pade-based algorithm. At most MAXSQRT matrix square roots are computed.
  # %   [X NSQ,M] = POWERM_TRIANG(T,P) returns the number NSQ of
  # %   square roots computed and the degree M of the Pade approximant used.
  #
  # if nargin < 3, maxsqrt = 100; end
  # if ~isequal(T,triu(T))
  #    error('MATLAB:powerm_triang:InputMat',...
  #          'First input argument must be an upper triangular matrix.');
  # end
  # n = length(T); nsq = 0;
  # if n == 1
  #    X = T^p; % compute the scalar power
  #    m = 0; return
  # elseif n == 2
  #    X = powerm2by2(T,p);  % Compute the power directly.
  #    m = 0; return
  # end
  #
  # T_old = T;
  # nsq = 0; q = 0;
  #
  # xvals = [ % Max norm(X) for degree m Pade approximant to (I-X)^p.
  #           1.512666672122460e-005    % m = 1
  #           2.236550782529778e-003    % m = 2
  #           1.882832775783885e-002    % m = 3 being used
  #           6.036100693089764e-002    % m = 4 being used
  #           1.239372725584911e-001    % m = 5 being used
  #           1.998030690604271e-001    % m = 6 being used
  #           2.787629930862099e-001    % m = 7 being used
  #           3.547373395551596e-001    % m = 8
  #           4.245558801949280e-001    % m = 9
  #           4.870185637611313e-001    % m = 10
  #           5.420549053918690e-001    % m = 11
  #           5.901583155235642e-001    % m = 12
  #           6.320530128774397e-001    % m = 13
  #           6.685149002867240e-001    % m = 14
  #           7.002836650662422e-001    % m = 15
  #           7.280253837034645e-001    % m = 16
  #           9.152924199170567e-001    % m = 32
  #           9.764341682154458e-001 ]; % m = 64
  #
  # while 1
  #
  #     normdiff = norm(T-eye(n),1);
  #
  #     if normdiff <= xvals(7)
  #
  #        q = q+1;
  #        j1 = find(normdiff <= xvals(3:7));
  #        j1 = j1(1) + 2;
  #        j2 = find(normdiff/2 <= xvals(3:7));
  #        j2 = j2(1) + 2;
  #        if j1-j2 < 2 || q == 2, m = j1; break, end
  #
  #     end
  #
  #     if nsq == maxsqrt, m = 16; break, end
  #     T = sqrtm_tri(T); nsq = nsq + 1;
  #
  # end
  #
  # X = powerm_cf(eye(n)-T,p,m);
  #
  # % Squaring phase, with directly computed diagonal and superdiagonal.
  #
  # for s = 0:nsq
  #     if s ~= 0, X = X*X; end  % Squaring
  #     for i = 1:n-1
  #         Tii = T_old(i:i+1,i:i+1);
  #         Si = powerm2by2(Tii,p/(2^(nsq-s)));
  #         X(i:i+1,i:i+1) = Si;
  #     end
  # end
  # end
  #
  # function R = sqrtm_tri(T)
  # %SQRTM_TRI   Upper triangular square root of upper triangular matrix.
  #
  # n = length(T);
  # R = zeros(n);
  # for j=1:n
  #     R(j,j) = sqrt(T(j,j));
  #     for i=j-1:-1:1
  #         R(i,j) = (T(i,j) - R(i,i+1:j-1)*R(i+1:j-1,j))/(R(i,i) + R(j,j));
  #     end
  # end
  # end
  #
  # function S = powerm_cf(Y,p,m)
  # %POWERM_CF   Evaluate Pade approximant bottom up.
  # %   POWERM_CF(Y,p,m) computes the [m/m] Pade approximant of the
  # %   matrix power (I-Y)^p by evaluating a continued fraction
  # %   representation in bottom-up fashion.
  #
  # k = 2*m;
  # n = length(Y);
  #
  # S = coeff(p,k)*Y;
  # for j=k-1:-1:1  % bottom-up
  #     S = coeff(p,j) * ( (eye(n) + S)\Y );
  # end
  # S = eye(n) + S;
  #
  # 	function c = coeff(p,i)
  # 	if i == 1, c = -p; return, end
  # 	jj = i/2;
  # 	if jj == round(jj)
  # 	   c = (-jj + p) / (2*(2*jj-1));
  # 	else
  # 	   jj = floor(jj);
  # 	   c = (-jj - p) / (2*(2*jj+1));
  # 	end
  # 	end
  #
  # end
  #
  # function X = powerm2by2(A,p)
  # %POWERM2BY2    Power of 2-by-2 upper triangular matrix.
  # %   POWERM2BY2(A,p) is the pth power of the 2-by-2 upper
  # %   triangular matrix A, where p is an arbitrary real number.
  #
  # a1 = A(1,1);
  # a2 = A(2,2);
  #
  # a1p = a1^p;
  # a2p = a2^p;
  #
  # loga1 = log(a1);
  # loga2 = log(a2);
  #
  # X = diag([a1p a2p]);
  #
  # if a1 == a2
  #
  #    X(1,2) = p*A(1,2)*a1^(p-1);
  #
  # elseif abs(a1) < 0.5*abs(a2) || abs(a2) < 0.5*abs(a1)
  #
  #    X(1,2) =  A(1,2) * (a2p - a1p) / (a2 - a1);
  #
  # else % Close eigenvalues.
  #
  #    w = atanh((a2-a1)/(a2+a1)) + 1i*pi*unwinding(loga2-loga1);
  #    dd = 2 * exp(p*(loga1+loga2)/2) * sinh(p*w) / (a2-a1);
  #    X(1,2) = A(1,2)*dd;
  #
  # end
  #
  # 	function u = unwinding(z,k)
  # 	%UNWINDING    Unwinding number.
  # 	%   UNWINDING(Z,K) is the K'th derivative of the
  # 	%   unwinding number of the complex number Z.
  # 	%   Default: k = 0.
  #
  # 	if nargin == 1 || k == 0
  # 	   u = ceil( (imag(z) - pi)/(2*pi) );
  # 	else
  # 	   u = 0;
  # 	end
  # 	end
  #
  # end

end
