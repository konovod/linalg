require "./spec_helper"

include Linalg
j = Complex.new(0, 1)
describe Linalg::Matrix do
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

  it "can be created from array of arrays" do
    m = GMat.new({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12} })
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

  it "scalars can be left member of multiplication, can be right member of division" do
    m1 = Mat.ones(3, 4)
    m = 3 * m1 / 2
    m.should eq GMat.new(3, 4) { |i, j| 1.5 }

    m1 = MatComplex.ones(3, 4)
    m = j * m1 / j
    m.should eq GMatComplex.new(3, 4) { |i, j| 1 }
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

  it "can be constructed using ##column and ##row" do
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
  pending "can be trasposed inplace (rectangular case)" do
    m = GMat.new([[1, 2, 3, 4], [5, 6, 7, 8]])
    m.transpose!
    m.should eq GMat.new([[1, 5], [2, 6], [3, 7], [4, 8]])
  end

  it "can be conjtrasposed" do
    m = GMatComplex.new([[1 + j, 2 + j, 3 - j], [4, 5, 6 + 2*j]])
    mt = m.conjtranspose
    mt.should eq GMatComplex.new([[1 - j, 4], [2 - j, 5], [3 + j, 6 - 2*j]])
    mt.conjtranspose.should eq m

    m = GMatComplex.new([[1, 2*j, 3], [4*j, 5, 6*j], [7, 8*j, 9 - j]])
    m.conjtranspose!
    m.should almost_eq GMatComplex.new([[1, -4*j, 7], [-2*j, 5, -8*j], [3, -6*j, 9 + j]])
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

  it "converted to array" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.to_a.should eq [1, 2, 3, 4, 5, 6]
  end
  it "converted to array of arrays" do
    m = GMat.new([[1, 2], [3, 4], [5, 6]])
    m.to_aa.should eq [[1, 2], [3, 4], [5, 6]]
  end
  it "can be created from matrix of another type" do
    m = GMat32.new([[1, 2], [3, 4], [5, 6]])
    GMatComplex.from(m).should eq GMatComplex.new([[1, 2], [3, 4], [5, 6]])
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
      m[-1..3, 2..3]
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
    a.cat!(b, 1)
    a.t.should eq c
  end
end
