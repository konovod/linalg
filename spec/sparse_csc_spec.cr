require "./spec_helper"

include LA
include LA::Sparse
describe CSCMatrix do
  it "can be created empty" do
    ms = CSCMatrix(Float64).new(3, 4, 5)
    ms.should eq GMat[
      [0, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 0, 0, 0],
    ]
  end

  it "can be constructed from arrays" do
    ms = CSCMatrix(Float64).new(4, 4, raw_columns: [0, 1, 2, 3, 4], raw_rows: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0, 6.0])
    ms.should eq GMat[
      [5, 0, 0, 0],
      [0, 8, 0, 6],
      [0, 0, 3, 0],
      [0, 0, 0, 0],
    ]
    expect_raises(ArgumentError) { CSCMatrix(Float64).new(4, 3, raw_columns: [0, 1, 2, 3, 4], raw_rows: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0, 6.0]) }
    expect_raises(ArgumentError) { CSCMatrix(Float64).new(4, 4, raw_columns: [0, 1, 2, 3, 4], raw_rows: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0]) }
  end

  it "can be created from dense matrix" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    ms = CSCMatrix(Float64).new(m)
    ms.should eq m
    ms.nonzeros.should eq 6
  end

  it "can be created from COO matrix" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    m = COOMatrix(Float64).new(m)
    ms = CSCMatrix(Float64).new(m)
    ms.should eq m
    ms.nonzeros.should eq 6
  end

  it "can be created from CSR matrix" do
    m_dense = GMat[
      [-4.0, -3.0, -2.0, 0.0],
      [-1.0, 0.0, 11.0, 0.0],
      [20.0, 30.0, 4.0, 1.0],
    ]
    m_csr = CSRMatrix(Float64).new m_dense
    m_csc = CSCMatrix(Float64).new m_csr
    m_csc.should eq m_dense
    m_csr.should eq m_dense
  end

  it "can be cleared" do
    ms = CSCMatrix(Float64).new(4, 4, raw_columns: [0, 1, 2, 3, 4], raw_rows: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0, 6.0])
    ms.clear
    ms.should eq GMat[
      [0, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 0, 0, 0],
    ]
  end

  it "can be compared with matrix when empty" do
    m = GMat.eye(5)
    ms = CSCMatrix(Float64).new(*m.shape)
    ms.should_not eq m
    m.should_not eq ms
    ms.should eq GMat.new(5, 5)
    ms.clear
    ms.should_not eq m
    m.should_not eq ms
    ms.should eq GMat.new(5, 5)
  end

  it "elements can be accessed and changed" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    ms = CSCMatrix(Float64).new(m)
    ms[3, 1].should eq 2
    ms[3, 2].should eq 0
    ms[3, 1] = 0.0
    ms[3, 1].should eq 0
    expect_raises(Exception) { ms[4, 2] = 0 }
  end

  it "can be created using `diag`" do
    ms = CSCMatrix(Float64).eye(3)
    ms.should eq GMat.new [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]
    ms = CSCMatrix(Float64).diag(4, 5, [1, 2, 3])
    ms.should eq GMat.new [
      [1, 0, 0, 0, 0],
      [0, 2, 0, 0, 0],
      [0, 0, 3, 0, 0],
      [0, 0, 0, 0, 0],
    ]
    ms.nonzeros.should eq 3
    expect_raises(Exception) { CSCMatrix(Float64).diag(4, 5, [1, 2, 3, 4, 5]) }
  end

  it "support select_with_index!" do
    m1 = COOMatrix(Float64).new(10, 10)
    m1[0, 0] = 10
    m1[1, 1] = -15
    m1[5, 2] = 11
    m1[3, 6] = 12
    m1 = CSCMatrix(Float64).new(m1)

    m2 = m1.to_dense
    m2[1, 1] = 0
    m2[0, 0] = 0
    m1.should_not eq m2
    m1.select_with_index! { |v, i, j| i != j }
    m1.should eq m2
  end

  it "support select!" do
    m1 = GeneralMatrix(Float64).new(10, 10)
    m1[0, 0] = 10
    m1[1, 1] = -15
    m1[5, 2] = 11
    m1[3, 6] = 12
    m1 = CSCMatrix(Float64).new(m1)

    m2 = m1.clone
    m2[1, 1] = 0
    m1.should_not eq m2
    m1.select! { |v| v >= 0 }
    m1.should eq m2
  end

  it "support map!" do
    a = CSCMatrix(Float64).new GMat[
      [1, 2, 0, 0],
      [0, 0, 5, 1],
      [0, 1, 0, 0],
    ]
    a.map! { |v| v/2 }
    a.should eq GMat[
      [0.5, 1, 0, 0],
      [0, 0, 2.5, 0.5],
      [0, 0.5, 0, 0],
    ]
    a.map_with_index { |v, row, col| v*2 + row + col }.should eq GMat[
      [1, 2 + 1, 0, 0],
      [0, 0, 5 + 3, 1 + 4],
      [0, 1 + 3, 0, 0],
    ]
  end

  it "support #tril! and #triu!" do
    m = COOMatrix(Float64).new(10, 10)
    m[0, 0] = 10
    m[1, 1] = 15
    m[5, 2] = 11
    m[3, 6] = 12
    m2 = CSCMatrix(Float64).new(m)
    m.triu!(1)
    m2.should_not eq m
    m2.triu!(1)
    m2.should eq m
  end

  it "support #tril and #triu" do
    m = COOMatrix(Float64).new(10, 10)
    m[0, 0] = 10
    m[1, 1] = 15
    m[5, 2] = 11
    m[3, 6] = 12
    m = CSCMatrix(Float64).new(m)

    m.tril(2).should eq m.to_dense.tril(2)
    m.tril.should eq m.to_dense.tril
    m.triu(2).should eq m.to_dense.triu(2)
    m.triu.should eq m.to_dense.triu

    ma = m.triu(1)
    m.should_not eq ma
    m.triu!(1)
    m.should eq ma
  end

  it "support resize!" do
    ma = GMat[
      [1, 2, 0],
      [0, 0, 5],
      [0, 1, 0],
      [-2, -1, -3],
    ]
    mas = CSCMatrix(Float64).new ma
    ma.resize!(3, 2)
    mas.resize!(3, 2)
    mas.should eq ma
    ma.resize!(10, 10)
    mas.resize!(10, 10)
    mas.should eq ma
  end

  it "can be multiplied to other CSC" do
    m1_dense = GMat[
      [5, 0, 0],
      [0, 8, 0],
      [0, 0, 3],
      [0, 6, 0],
    ]
    m1 = CSCMatrix(Float64).new m1_dense
    m2_dense = GMat[
      [1, 0, 3, 0, 5],
      [1, 2, 0, 1, 0],
      [0, 2, 0, -6, 0],
    ]
    m2 = CSCMatrix(Float64).new m2_dense

    (m1*m2).should eq m1_dense*m2_dense
    (m1*m2).nonzeros.should eq 11
    (m1*m2).flags.should eq MatrixFlags::None
    expect_raises(ArgumentError) { m2*m1 }
  end

  it "can be transposed" do
    ma = CSCMatrix(Float32).new GMat[
      [1, 2, 0, 0],
      [0, 0, 5, 1],
      [0, 1, 0, 0],
    ]
    mat = ma.t
    mat.should be_a CSCMatrix(Float32)
    mat.should eq GMat[
      [1, 0, 0],
      [2, 0, 1],
      [0, 5, 0],
      [0, 1, 0],
    ]
    mat.transpose.should eq ma
  end
  it "can be conjtransposed" do
    ma = CSCMatrix(Complex).new GMatComplex[
      [1, 2.i, 0, 0],
      [0, 0, -5.i, 1],
      [0, 1 - 3.i, 0, 0],
    ]
    mat = ma.conjt
    mat.should be_a CSCMatrix(Complex)
    mat.should eq GMatComplex[
      [1, 0, 0],
      [-2.i, 0, 1 + 3.i],
      [0, 5.i, 0],
      [0, 1, 0],
    ]
    mat.conjtranspose.should eq ma
  end

  it "can be added to other CSC" do
    a = CSCMatrix(Float64).new GMat[
      [1, 2, 0, 0],
      [0, 0, 5, 1],
      [0, 1, 0, 0],
    ]
    b = CSCMatrix(Float64).new GMat[
      [0, -2, 0, 0],
      [0, 1, -5, 1],
      [0, 2, 0, 1],
    ]
    (a + b).should eq GMat[
      [1, 0, 0, 0],
      [0, 1, 0, 2],
      [0, 3, 0, 1],
    ]
    (a + b).nonzeros.should eq 5
    expect_raises(ArgumentError) { a.t + b }
    (-a + b + a - b).nonzeros.should eq 0
  end

  it "support norms" do
    m = CSCMatrix(Float64).new GMat[
      [-4.0, -3.0, -2.0],
      [-1.0, 0.0, 1.0],
      [2.0, 3.0, 4.0],
    ]
    m.norm.should be_close(7.745966692414834, 1e-6)
    m.norm(MatrixNorm::Frobenius).should be_close(7.745966692414834, 1e-6)
    m.norm(MatrixNorm::Inf).should eq 9
    m.norm(MatrixNorm::One).should eq 7
    m.norm(MatrixNorm::MaxAbs).should eq 4
  end

  it "can be converted to CSR using as_csr" do
    m = CSCMatrix(Float64).new GMat[
      [-4.0, -3.0, -2.0],
      [-1.0, 0.0, 1.0],
      [2.0, 3.0, 4.0],
    ]
    mt = m.as_csr
    mt.should be_a CSRMatrix(Float64)
    mt.should eq m.t
  end
end
