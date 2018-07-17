require "./spec_helper"

include LA
describe LA do
  it "calls LAPACK functions directly" do
    m = [
      1.0, 0.0, 1.0,
      0.0, 4.0, 0.0,
      0.0, 0.0, 1.0,
    ]
    ipiv = Slice(Int32).new(3)
    LibLAPACKE.dgetrf(LibCBLAS::COL_MAJOR, 3, 3, m, 3, ipiv).should eq 0
    LibLAPACKE.dgetri(LibCBLAS::COL_MAJOR, 3, m, 3, ipiv).should eq 0
    m.should eq [1.0, 0.0, -1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 1.0]
  end

  it "calls functions using matrix class" do
    matrix1 = GMat[
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ]
    matrix2 = matrix1*1
    ipiv = Slice(Int32).new(3)
    LibLAPACKE.dgetrf(LibCBLAS::COL_MAJOR, 3, 3, matrix2, 3, ipiv).should eq 0
    LibLAPACKE.dgetri(LibCBLAS::COL_MAJOR, 3, matrix2, 3, ipiv).should eq 0
    (matrix1*matrix2).should eq Mat.identity(3)
  end

  it "calls functions using high level wrapper" do
    matrix1 = GMat[
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ]
    (matrix1*matrix1.inv).should eq Mat.identity(3)
  end

  it "ultra-highlevel - calls functions on a submatrices" do
    m4 = Mat.ones(4, 4) + Mat.diag([1, 2, 3, 4])
    m3 = m4[0..2, 0..2]
    (m3*m3.inv).should almost_eq Mat.identity(3)
  end

  it "raises LinAlgError on incorrect data" do
    matrix1 = GMat[
      [1, 0, 1],
      [0, 0, 0],
      [0, 0, 1],
    ]
    expect_raises(LinAlgError) do
      matrix1.inv!
    end
  end

  it "support all types" do
    matrix1 = GMat32[
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ]
    (matrix1*matrix1.inv).should eq Mat32.identity(3)

    matrix1 = GMatComplex[
      {1 + 1.i, 0, 1},
      {0, 4, 0},
      {0, 0, 1},
    ]
    (matrix1*matrix1.inv).should eq MatComplex.identity(3)
  end

  it "correctly update symmetric matrices" do
    matrix1 = GMat[
      [1, 1, 1],
      [1, 2, 3],
      [1, 3, 6],
    ]
    matrix1.detect
    matrix1.inv.should eq GMat[
      [3.0, -3.0, 1.0],
      [-3.0, 5.0, -2.0],
      [1.0, -2.0, 1.0],
    ]
  end

  it "high-level: solve linear equations" do
    a = GMat32[
      [2, 4],
      [2, 8]]
    b = GMat32[[2], [4]]
    LA.solve(a, b).should eq (a.inv * b)
  end

  it "high-level: calculate determinant" do
    a = GMat[
      [1, 2, 3],
      [4, 5, 7],
      [-1, 1, -1]]
    a.det.should eq 9

    a = GMatComplex[
      [1.i, 2, 3],
      [0, 5, 7],
      [0, 0, -1.i]]
    a.det.should eq 5
  end
  it "high-level: solve linear least square" do
    a = GMat32[
      [1, 2, 0],
      [0, 4, 3]]
    b = GMat32[[8], [18]]
    x = LA.solvels(a, b)
    x_octave = GMat32.columns([0.918032, 3.54098, 1.27869])
    x.should be_close(x_octave, 1e-3)
  end

  it "high-level: solve linear least square (complex)" do
    a = GMatComplex[
      [1, 2, 0],
      [0, 4, 3]]
    b = GMatComplex[[8], [18]]
    x, rank, s = LA.lstsq(a, b, LSMethod::SVD)
    x_octave = GMatComplex.columns([0.918032, 3.54098, 1.27869])
    x.should be_close(x_octave, 1e-5)
    rank.should eq 2
    s[0].should be_close 5.273163, 1e-5
    s[1].should be_close 1.481132, 1e-5
  end

  it "high-level: calculate singular value decomposition" do
    a = GMat[[1, 2, 3], [4, 5, 6]]
    u, s, vt = LA.svd(a)
    (u*Mat.diag(a.nrows, a.ncolumns, s)*vt).should almost_eq a
    u.detect(MatrixFlags::Orthogonal).should be_true
    vt.detect(MatrixFlags::Orthogonal).should be_true

    ac = GMatComplex[[1, 2, 3], [4, 5, 6]]
    s1 = ac.svdvals(overwrite_a: true)
    s1[0].should be_close s[0], 1e-6
    s1[1].should be_close s[1], 1e-6
  end

  # TODO - proper spec
  it "high-level: balance matrix" do
    a = GMat32[[1, 2, 0], [9, 1, 0.01], [1, 2, 10*Math::PI]]
    b, s = a.balance(separate: true)
    s.should eq GMat32[[0.5, 1, 1]]
  end

  it "high-level: LU factorization" do
    a = GMat32[
      [1, 2, 3, 4],
      [5, 6, 7, 8],
      [1, 2, 10, 0],
      [50, 6, 7, 8],
      [1, 21, 1, 0]]
    p, l, u = a.lu
    (p*l*u).should almost_eq a
    l.detect(MatrixFlags::LowerTriangular).should be_true
    l.detect(MatrixFlags::UpperTriangular).should be_false
    u.detect(MatrixFlags::UpperTriangular).should be_true
  end

  it "high-level: solve using LU" do
    a = GMat32[
      [2, 4],
      [2, 8]]
    lu = a.lu_factor
    b = GMat32[[2], [4]]
    lu.solve(b).should eq (a.inv * b)
  end

  it "high-level: cholesky decomposition" do
    a = GMatComplex[[1, -2.i], [2.i, 5]]
    c = a.cholesky(lower: true)
    (c*c.conjtranspose).should almost_eq a
  end
  it "high-level: using cholesky decomposition to solve equations" do
    a = GMatComplex[[1, -2.i], [2.i, 5]]
    chol1 = a.cholesky(lower: true, dont_clean: true)
    chol2 = a.cholesky(lower: false, dont_clean: false)
    b = GMatComplex[[2], [4]]
    x1 = chol1.cho_solve(b)
    x2 = chol1.cho_solve(b)
    x1.should eq x2
    (a*x1 - b).should eq MatComplex.zeros(2, 1)
  end

  it "high-level: hessenberg decomposition" do
    a = GMatComplex[[1, -2.i, 3], [2.i, 5, 4], [7, 0, 1.i]]
    h, q = a.hessenberg(calc_q: true)
    q.detect(MatrixFlags::Orthogonal).should be_true
    (q*h*q.conjtranspose).should almost_eq a

    a = GMat[[1, -2, 3], [2, 5, 4], [7, 0, 1]]
    h, q = a.hessenberg(calc_q: true)
    q.detect(MatrixFlags::Orthogonal).should be_true
    (q*h*q.transpose).should almost_eq a
  end

  it "high-level: schur decomposition (real argument)" do
    a = GMat[[1, -2, 3], [2, 5, 4], [7, 0, 1]]
    t, z = a.schur
    z.detect(MatrixFlags::Orthogonal).should be_true
    (z*t*z.transpose).should almost_eq a
  end

  it "high-level: schur decomposition (complex argument)" do
    a = GMatComplex[[1, -2.i, 3], [2.i, 5, 4], [7, 0, 1.i]]
    t, z = a.schur
    z.detect(MatrixFlags::Orthogonal).should be_true
    (z*t*z.conjtranspose).should almost_eq a
  end

  it "high-level: qr decomposition" do
    a = GMatComplex[[1, -2.i, 3], [2.i, 5, 4], [7, 0, 1.i]]
    q, r, pvt = a.qr
    q.detect(MatrixFlags::Orthogonal).should be_true
    (q*r).should almost_eq a
    # only r
    r1, pvt = a.qr_r
    r1.should almost_eq r
    # raw form
    q1, tau, pvt = a.qr_raw
    # pivoted TODO - spec?
    q, r, pvt = a.qr(pivoting: true)
  end

  it "high-level: rq decomposition" do
    a = GMatComplex[[1, -2.i, 3], [2.i, 5, 4], [7, 0, 1.i]]
    r, q = a.rq
    q.detect(MatrixFlags::Orthogonal).should be_true
    (r*q).should almost_eq a
    # only r
    r1 = a.rq_r
    r1.should almost_eq r
  end

  it "high-level: lq and ql decomposition" do
    a = GMatComplex[[1, -2.i, 3], [2.i, 5, 4], [7, 0, 1.i]]
    l, q = a.lq
    q.detect(MatrixFlags::Orthogonal).should be_true
    (l*q).should almost_eq a

    a = GMat32[[1, -2, 3], [2, 5, 4], [7, 0, 1]]
    q, l = a.ql
    q.detect(MatrixFlags::Orthogonal).should be_true
    (q*l).should almost_eq a
  end

  it "high-level: qz decomposition (real argument)" do
    a = GMat[
      [-2, 4, 1],
      [2, -4, 1],
      [1, 1, 1],
    ]
    b = GMat[
      [-2, 4, 0],
      [2, -3, 1],
      [2, 1, 1],
    ]
    aa, bb, vl, vr = LA.qz(a, b)
    vl.detect(MatrixFlags::Orthogonal).should be_true
    vr.detect(MatrixFlags::Orthogonal).should be_true
    (vl*aa*vr.transpose).should almost_eq a
    (vl*bb*vr.transpose).should almost_eq b
  end

  it "high-level: qz decomposition (complex argument)" do
    a = GMatComplex[
      [-2, 4.i, 1],
      [3, -4.i, 1],
      [1, 0, 1 - 1.i],
    ]
    b = GMatComplex[
      [-2, 4, 0],
      [2, -3, -2.i],
      [1, 1 - 1.i, 1 + 1.i],
    ]
    aa, bb, vl, vr = LA.qz(a, b)
    vl.detect(MatrixFlags::Orthogonal).should be_true
    vr.detect(MatrixFlags::Orthogonal).should be_true
    (vl*aa*vr.conjtranspose).should almost_eq a
    (vl*bb*vr.conjtranspose).should almost_eq b
  end

  it "high-level: calculate matrix norms" do
    a = GMat[[-4, -3, -2, -1, 0, 1, 2, 3, 4]].t
    b = a.reshape(3, 3)
    a.norm.should be_close(7.745966692414834, 1e-6)
    b.norm.should be_close(7.745966692414834, 1e-6)
    b.norm(MatrixNorm::Frobenius).should be_close(7.745966692414834, 1e-6)
    a.norm(MatrixNorm::Inf).should eq 4
    b.norm(MatrixNorm::Inf).should eq 9
    a.norm(MatrixNorm::One).should eq 20
    b.norm(MatrixNorm::One).should eq 7
  end

  it "high-level: calculate matrix rank" do
    m = GMat[
      [1, 0, 1],
      [0, 4, 0],
      [2, 3, 2],
    ]
    m.rank.should eq 2
    m.rank(method: RankMethod::QRP).should eq 2
    m.rank(method: RankMethod::SVD, overwrite_a: true).should eq 2
    m = 1e-6 * GMat[
      [1, 0, 1, 4],
      [0, 4, 0, 4],
      [1, 3, 2, 10],
    ]
    m.rank.should eq 3
  end

  it "high-level: calculate matrix rank (edge cases)" do
    m = Mat.ones(5, 1)
    m.rank.should eq 1
    m = Mat.zeros(5, 1)
    m.rank.should eq 0
    m = Mat.zeros(8, 6)
    m.rank.should eq 0
    m = Mat.ones(8, 6)
    m.rank(method: RankMethod::QRP).should eq 1
  end
end
