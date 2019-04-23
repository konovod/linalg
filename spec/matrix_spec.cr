require "./spec_helper"

include LA
describe LA::Matrix do
  it "can be created with given size" do
    m = GMat.new(10, 15)
    m.raw.size.should eq 10*15
    m.raw[0].should eq 0
  end

  it "can be created from array with given dimension" do
    m = GMat.new(3, 4, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    m.nrows.should eq 3
    m[1, 0].should eq 5
  end

  it "can be created from array with given dimension(column major)" do
    m = GMat.new(3, 4, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], true)
    m.nrows.should eq 3
    m[1, 0].should eq 2
  end

  it "can be created from array of arrays" do
    m = GMat[{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}]
    m.nrows.should eq 4
    m.ncolumns.should eq 3
    m[1, 0].should eq 4
    expect_raises(IndexError) { GMat.new({ {1, 2, 3}, {4, 5} }) }
  end

  it "can be created from a dimensions and a block" do
    m = GMat32.new(3, 3) { |i, j| i*3 + j }
    m[1, 1].should eq 4
  end

  it "can't create matrix of unsupported type" do
    # commented as it causes a compile error, as it should
    # m = Matrix(Int32).new(3, 3)
  end

  it "access to members works and is zero-based" do
    m = GMat.new(5, 4)
    m[4, 3] = 1.0
    m[4, 3].should eq 1.0
    expect_raises(IndexError) { m[5, 4] = 1.0 }
    expect_raises(IndexError) { m[0, 6] = 1.0 }
  end

  it "negative indices counts from size" do
    m = GMat.new(5, 4)
    m[3, -2] = 1.0
    m[-2, 2].should eq 1.0
    m[-2, -2].should eq 1.0
    m[3, -1].should eq 0.0
    expect_raises(IndexError) { m[3, -5] = 1.0 }
  end

  it "can multiply matrices" do
    m1 = GMat.new(3, 4) { |i, j| i + j }
    m2 = GMat.new(4, 2) { |i, j| i - j }
    m = m1*m2
    m.should eq GMat.new(3, 2, [14.0, 8.0, 20.0, 10.0, 26.0, 12.0])
    expect_raises(ArgumentError) { GMat.new(3, 4) * GMat.new(3, 4) }
  end

  it "can do sum and scalar multiply" do
    m1 = GMat.new(3, 4) { |i, j| i }
    m2 = GMat.new(3, 4) { |i, j| j }
    m = m1 + m2*2
    m.should eq GMat.new(3, 4) { |i, j| i + 2*j }
    expect_raises(ArgumentError) { GMat.new(3, 4) + GMat.new(4, 4) }
  end

  it "can be multiplied to complex number" do
    m1 = Mat32.ones(3, 4)
    (m1*(2 + 1.i)).should eq GMatComplex.new(3, 4) { 2 + 1.i }
  end

  it "scalars can be left member of multiplication, can be right member of division" do
    m1 = Mat.ones(3, 4)
    m = 3 * m1 / 2
    m.should eq GMat.new(3, 4) { |i, j| 1.5 }

    m1 = MatComplex.ones(3, 4)
    m = 2.i * m1 / 1.i
    m.should eq GMatComplex.new(3, 4) { |i, j| 2 }
  end

  it "can be added to scalar" do
    m = GMat.eye(2)
    m1 = m + 1
    m1.should eq GMat[[2, 1], [1, 2]]
    m1.flags.should eq MatrixFlags::None

    (m1 - 1.0).should eq m
    (1.i + m).should eq GMatComplex[[1 + 1.i, 1.i], [1.i, 1 + 1.i]]
  end

  it "don't reset flags if added to zero" do
    m = GMat.eye(5)
    m2 = (2.0 - 2) - m
    m2.should eq -m
    m2.flags.should eq (-m).flags
  end

  it "can checks if it is square" do
    GMat.new(3, 4).square?.should be_false
    GMatComplex.new(30, 30).square?.should be_true
  end

  it "can be initialized with zeros and ones" do
    m = Mat.identity(3).should eq GMat.new([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    m = Mat.zeros(2, 2).should eq GMat.new([[0, 0], [0, 0]])
    m = Mat.ones(1, 3).should eq GMat.new([[1, 1, 1]])
  end

  it "can be initialized with random vales" do
    rng1 = Random.new(1)
    m1 = MatComplex.rand(5, 5, rng1)
    rng2 = Random.new(1)
    m2 = MatComplex.rand(5, 5, rng2)
    m3 = MatComplex.rand(5, 5)
    m1.should eq m2
    m1.should_not eq m3
  end

  it "can be constructed using repmat" do
    a = GMat.new([[1, 2]])
    b = a.repmat(5, 3)
    b.should eq GMat.new([
      [1, 2, 1, 2, 1, 2],
      [1, 2, 1, 2, 1, 2],
      [1, 2, 1, 2, 1, 2],
      [1, 2, 1, 2, 1, 2],
      [1, 2, 1, 2, 1, 2],
    ])
    c = Mat.repmat(a, 2, 1)
    c.should eq GMat.new([
      [1, 2],
      [1, 2],
    ])
  end

  it "can be constructed with diagonal elements" do
    Mat32.diag(2, 2, 5).should eq GMat32.new([[5, 0], [0, 5]])
    Mat32.diag(3, 2) { |i| -i - 1 }.should eq GMat32.new([[-1, 0], [0, -2], [0, 0]])
    Mat32.diag(2, 3, [14, 15]).should eq GMat32.new([[14, 0, 0], [0, 15, 0]])
  end

  it "can be constructed using #column and #row" do
    Mat32.column([2, 2, 5]).should eq GMat32.new([[2], [2], [5]])
    Mat32.row([2, 2, 5]).should eq GMat32.new([[2, 2, 5]])
  end

  it "can be trasposed" do
    m = GMat.new([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    mt = m.transpose
    mt.should eq GMat.new([[1, 4, 7], [2, 5, 8], [3, 6, 9]])

    m = GMat.new([[1, 2, 3, 4], [5, 6, 7, 8]])
    mt = m.transpose
    mt.should eq GMat.new([[1, 5], [2, 6], [3, 7], [4, 8]])
  end
  it "can be trasposed inplace (square case)" do
    m = GMat.new([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    m.transpose!
    m.should eq GMat.new([[1, 4, 7], [2, 5, 8], [3, 6, 9]])
  end
  it "can be trasposed inplace (rectangular case)" do
    m = GMat.new([[1, 2, 3, 4], [5, 6, 7, 8]])
    m.transpose!
    m.should eq GMat.new([[1, 5], [2, 6], [3, 7], [4, 8]])
  end

  it "can be conjtrasposed" do
    m = GMatComplex.new([[1 + 1.i, 2 + 1.i, 3 - 1.i], [4, 5, 6 + 2.i]])
    mt = m.conjtranspose
    mt.should eq GMatComplex.new([[1 - 1.i, 4], [2 - 1.i, 5], [3 + 1.i, 6 - 2.i]])
    mt.conjtranspose.should eq m

    m = GMatComplex.new([[1, 2.i, 3], [4.i, 5, 6.i], [7, 8.i, 9 - 1.i]])
    m.conjtranspose!
    m.should almost_eq GMatComplex.new([[1, -4.i, 7], [-2.i, 5, -8.i], [3, -6.i, 9 + 1.i]])
  end

  it "has kron operation" do
    a = GMat.new([[1, 2], [3, 4]])
    b = GMat.new([[1, -1, 1]])
    Mat.kron(a, b).should eq GMat.new([
      [1, -1, 1, 2, -2, 2],
      [3, -3, 3, 4, -4, 4],
    ])
  end

  it "have tril and triu functions" do
    a = GMat32.new([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
    a.tril(-1).should eq GMat32.new([
      [0, 0, 0],
      [4, 0, 0],
      [7, 8, 0],
      [10, 11, 12],
    ])
    a.triu.should eq GMat32.new([
      [1, 2, 3],
      [0, 5, 6],
      [0, 0, 9],
      [0, 0, 0],
    ])
    a.triu(-1).should eq GMat32.new([
      [1, 2, 3],
      [4, 5, 6],
      [0, 8, 9],
      [0, 0, 12],
    ])
  end

  it "can be reshaped" do
    a = GMat.new([[1, 2, 3, 4, 5, 6]])
    a.reshape!(3, 2)
    a.should eq GMat.new([[1, 2], [3, 4], [5, 6]])
    a.reshape(2, 3).should eq GMat.new([[1, 2, 3], [4, 5, 6]])
    expect_raises(ArgumentError) do
      a.reshape(1, 4)
    end
  end

  it "can be reshaped (column major)" do
    a = GMat.new([[1, 2, 3, 4, 5, 6]])
    a.reshape!(3, 2, true)
    a.should eq GMat.new([[1, 4], [2, 5], [3, 6]])
    a.reshape(2, 3, true).should eq GMat.new([[1, 3, 5], [2, 4, 6]])
    expect_raises(ArgumentError) do
      a.reshape(1, 4, true)
    end
  end

  it "have tri function" do
    Mat.tri(3, 5, 2).should eq GMat.new(
      [[1, 1, 1, 0, 0],
       [1, 1, 1, 1, 0],
       [1, 1, 1, 1, 1]])
  end

  it "printed to string" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.to_s.should eq "\n[1.0, 2.0]\n[3.0, 4.0]\n[5.0, 6.0]\n\n"
  end

  it "inspected" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.inspect.should eq "LA::GeneralMatrix(Float64) (3x2, None):\n[1.0, 2.0]\n[3.0, 4.0]\n[5.0, 6.0]\n\n"
  end

  it "converted to array" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.to_a.should eq [1, 2, 3, 4, 5, 6]
  end
  it "converted to array of arrays" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.to_aa.should eq [[1, 2], [3, 4], [5, 6]]
  end

  it "converted to array (column major)" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.to_a(true).should eq [1, 3, 5, 2, 4, 6]
  end
  it "converted to array of arrays (column major)" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.to_aa(true).should eq [[1, 3, 5], [2, 4, 6]]
  end

  it "can be created from matrix of another type" do
    m = GMat32.new([[1, 2], [3, 4], [5, 6]])
    GMatComplex.new(m).should eq GMatComplex.new([[1, 2], [3, 4], [5, 6]])
  end

  it "can provide submatrices" do
    m = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
    m1 = m[1..3, 2..3]
    m1[0, 0].should eq m[1, 2]
    m1[2, 1] = 10
    m[3, 3].should eq 10
    expect_raises(IndexError) do
      m1[1, 3]
    end
    expect_raises(IndexError) do
      m[-6..3, 2..3]
    end
    expect_raises(IndexError) do
      m[0..3, 2..4]
    end
    m[1, 2..3].should eq GMat32.new([[7, 8]])
    m[1..2, 3].should eq GMat32.new([[8], [12]])
  end

  it "has row and column functions" do
    m = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    r1 = m.rows[2]
    r1.size.to_a.should eq [1, 4]
    r1[0, 2].should eq m[2, 2]
    r1 = m.columns[3]
    r1.size.to_a.should eq [3, 1]
    r1[1, 0].should eq m[1, 3]
  end

  it "can override submatrices" do
    m = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    m[1, 2..3] = 10
    m[1..2, 0] = 20
    m[1, 0..3] = m.rows[2]
    m[0..2, 1] = -1
    m[0..2, 3] = m.columns[2] - m.rows[1].transpose[0..2, 0]
    m.should eq GMat32.new [
      [1, -1, 3, -17],
      [20, -1, 11, 12],
      [20, -1, 11, 0],
    ]
  end

  it "can be concatenated with other matrices" do
    a = Mat.ones(2, 3)*2
    b = Mat.eye(3)
    c = a.vcat(b)
    a.vcat(b).should eq GMat.new [
      [2, 2, 2],
      [2, 2, 2],
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]
    expect_raises(ArgumentError) { a.hcat(b) }
    a = a.t
    a.cat!(b, Axis::Rows)
    a.t.should eq c
  end

  it "columns and rows can be extracted by range" do
    m = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    m.rows[1..2].should eq m[1..2, 0..3]
    m.columns[1...3].should eq m[0..2, 1..2]
  end

  it "have map and map_with_index methods" do
    m = Mat.rand(10, 15)
    m2 = m.map { |x| -x }
    (m + m2).should eq Mat.zeros(10, 15)

    m3 = m.map_with_index { |x, row, col| row < col ? -x : 0 }
    (m + m3).detect?(MatrixFlags::LowerTriangular).should be_true
  end

  it "have map! and map_with_index! methods" do
    m = Mat.tri(10, 15) * (-1.0)
    m.map! { |x| -x }.should eq Mat.tri(10, 15)
    m.should eq Mat.tri(10, 15)

    m = Mat.ones(10, 15)
    m.map_with_index! { |x, row, col| row < col ? x : 0 }
    m.detect?(MatrixFlags::UpperTriangular).should be_true
  end

  it "can evaluate trace" do
    m = MatComplex.diag([1, 2, 3, 4, 5])
    m.trace.should eq (1 + 2 + 3 + 4 + 5)
    m = Mat.ones(10, 20)
    m.trace.should eq 10
  end

  it "have Indexable diag" do
    m = GMat32.new([[-1, 2, 3, 4], [5, -6, 7, 8], [9, 10, -11, 12]])
    m.diag[1].should eq -6
    (m.diag.all? &.<(0)).should be_true
  end

  it "have Indexable diag with offset" do
    m = GMat32.new([[-1, 2, 3, 4], [5, -6, 7, 8], [9, 10, -11, 12]])
    m.diag(0).to_a.should eq [-1, -6, -11]
    m.diag(1).to_a.should eq [2, 7, 12]
    m.diag(2).to_a.should eq [3, 8]
    m.diag(3).to_a.should eq [4]
    expect_raises(ArgumentError) { m.diag(4) }
    m.diag(-1).to_a.should eq [5, 10]
    m.diag(-2).to_a.should eq [9]
    expect_raises(ArgumentError) { m.diag(-3) }
    expect_raises(ArgumentError) { m.diag(-4) }
  end

  it "don't crash when using submatrix of submatrix" do
    m4 = Mat.ones(4, 4)
    m3 = m4[0..2, 0..2]
    m2 = m3[0..1, 0..1]
    m2.should eq GMat[[1, 1], [1, 1]]
  end

  it "can be `resize!`d" do
    m = Mat.ones(2, 3)
    m_same = m.resize!(3, 2)
    m_same.should eq GMat[
      [1.0, 1.0],
      [1.0, 1.0],
      [0.0, 0.0]]
    m.resize!(2, 3)
    m.should eq GMat[
      [1.0, 1.0, 0.0],
      [1.0, 1.0, 0.0]]
  end

  it "can `add!` matrices" do
    m = GMat.eye(4)
    m.add! m
    m.should eq GMat.diag([2, 2, 2, 2])
    m2 = Mat.rand(4, 4)
    m.add! -1, m2
    (m + m2).should eq GMat.diag([2, 2, 2, 2])

    expect_raises(ArgumentError) { m.add!(GMat.rand(3, 4)) }
  end

  it "can `scale!` self by scalar" do
    m = Mat.rand(5, 8)
    m_half = m.clone
    m_half.scale! 0.5
    (m_half + m_half).should eq m
  end

  it "submatrices can be `map!`ped too" do
    m = Mat.ones(3, 2)
    x = m[0..1, 0..0]
    x.map! { |v| -v }
    m.should eq GMat[
      [-1.0, 1.0],
      [-1.0, 1.0],
      [1.0, 1.0],
    ]
  end

  it "have `each` methods" do
    m = GMat[[1, 2, 3], [4, 5, 6]]
    elements = [] of Float64
    m.each { |v| elements << v }
    elements.should eq [1, 2, 3, 4, 5, 6]

    elements1 = [] of Float64
    m.each(all: true) { |v| elements1 << v }
    elements1.should eq elements

    elements1 = [] of Float64
    m.each(all: false) { |v| elements1 << v }
    elements1.should eq elements

    m2 = Mat.zeros(*m.size)
    m.each_with_index { |v, i, j| m2[i, j] = v }
    m2.should eq m

    m.each_index(all: false) { |i, j| m[i, j] = -m[i, j] }
    m.should eq -m2
  end

  it "have `reduce` method" do
    m = GMat32[[1, 2, 5], [10, 4, 0.025]]
    m.reduce(Axis::Rows, 1) { |memo, x| memo/x }.should eq GMat32[[0.1], [1]]
    m.reduce(Axis::Columns, 0) { |memo, x| memo - x }.should eq GMat32[[-11, -6, -5.025]]
  end

  it "have `sum` and other aggregation methods" do
    m = GMat[[0.1, 0.2, 0.3], [10, 20, 30]]
    m.sum(Axis::Rows).should be_close GMat[[0.6], [60]], 1e-9
    m.sum(Axis::Columns).should be_close GMat[[10.1, 20.2, 30.3]], 1e-9
  end

  it "array of matrices can be summed without initial value" do
    [GMat.ones(3, 4), GMat.ones(3, 4)].sum.should eq 2*GMat.ones(3, 4)
    expect_raises(ArgumentError) { [GMat.ones(3, 4), GMat.ones(4, 3)].sum }
  end

  it "array of matrices can be multiplied without initial value" do
    [GMat.ones(3, 4), GMat.ones(4, 3)].product.should eq 4*GMat.ones(3, 3)
    expect_raises(ArgumentError) { [GMat.ones(3, 4), GMat.ones(3, 4)].product }
  end

  it "empty array of matrices sums to a scalar zero" do
    x = ([] of GMat32).sum
    typeof(x).should eq (Float32 | GMat32)
    x.as(Float32).should eq 0
  end

  it "empty array of matrices multiplies to a scalar one" do
    ([] of GMatComplex).product.as(Complex).should eq 1.0 + 0.i
  end

  it "complex matrix can be chopped to real" do
    GMatComplex[[1, 2, 5], [10, 4, 1e-6]].chop.should be_a (GMat | Nil)
    GMatComplex[[1, 2, 5], [10, 4, 1e-6]].chop.not_nil!.should eq GMat[[1, 2, 5], [10, 4, 1e-6]]
    GMatComplex[[1, 2, 5.i], [10, 4, 1e-6]].chop.should be_a Nil
    GMatComplex[[1, 2, 5], [10, 4, 1e-6 + 1e-6.i]].chop.should be_a Nil
    GMatComplex[[1, 2, 5], [10, 4, 1e-6 + 1e-6.i]].chop(1e-3).not_nil!.should eq GMat[[1, 2, 5], [10, 4, 1e-6]]
    GMatComplex[[1, 2, 5], [10, 4, 1e-6 + 1e-26.i]].chop.not_nil!.should eq GMat[[1, 2, 5], [10, 4, 1e-6]]
  end
end
