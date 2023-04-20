require "./spec_helper"

include LA
include LA::Sparse
describe COOMatrix do
  it "can be created from dense matrix" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    ms = COOMatrix(Float64).new(m)
    ms.should eq m
    ms.nonzeros.should eq 6
  end

  it "can be constructed from arrays" do
    ms = COOMatrix(Float64).new(3, 4, [0, 1, 2], [1, 1, 2], [5.0, 10.0, 15.0])
    ms.should eq GMat[
      [0, 5, 0, 0],
      [0, 10, 0, 0],
      [0, 0, 15, 0],
    ]
    expect_raises(ArgumentError) { COOMatrix(Float64).new(3, 3, [0, 1, 3], [1, 1, 2], [5.0, 10.0, 15.0]) }
    expect_raises(ArgumentError) { COOMatrix(Float64).new(3, 3, [0, 1, 2], [1, 1], [5.0, 10.0, 15.0]) }
  end

  it "can be created from another COO matrix" do
    m = GMat.eye(5)
    ms1 = COOMatrix(Float32).new(m)
    ms1[3, 1] = 1.0
    ms2 = COOMatrix(Float32).new(ms1)
    ms2.should eq ms1
    ms1[3, 1] = 2.0
    ms2.should_not eq ms1

    ms3 = COOMatrix(Float64).new(ms1)
    ms3.should eq ms1
  end

  it "can be created empty" do
    m = COOMatrix(Float64).new(2, 5)
    m.flags.diagonal?.should be_true
    m[0, 0] = 1.0
    m[1, 2] = 5.0
    m.should eq GMat[
      [1, 0, 0, 0, 0],
      [0, 0, 5, 0, 0],
    ]
    m.flags.none?.should be_true
  end

  it "elements can be added and removed" do
    m = COOMatrix(Float64).new(4, 2)
    m.nonzeros.should eq 0
    m[3, 1] = 1.0
    m.nonzeros.should eq 1
    m[3, 1] = 2.0
    m.nonzeros.should eq 1
    m.should eq GMat[
      [0, 0],
      [0, 0],
      [0, 0],
      [0, 2],
    ]
    m[3, 1] = 0.0
    m.nonzeros.should eq 0
    m.should eq GMat.new(4, 2)
  end

  it "can be cleared" do
    m = COOMatrix(Float64).new(GMat32.eye(5))
    m.clear
    m.nonzeros.should eq 0
    m.should eq GMat32.new(5, 5)
  end

  it "can be transposed" do
    m = COOMatrix(Float64).new(GMat32.eye(5))
    m[3, 2] = 20.0
    mt = m.transpose
    mt[3, 2].should eq 0
    mt[2, 3].should eq 20.0
    mt.should_not eq m
    m.transpose!
    mt.should eq m
  end

  it "can be initialized with defined number of random values" do
    rng1 = Random.new(1)
    m1 = COOMatrix(Complex).rand(50, 50, nonzero_elements: 4, rng: rng1)
    m1.nonzeros.should eq 4
    rng2 = Random.new(1)
    m2 = COOMatrix(Complex).rand(50, 50, nonzero_elements: 4, rng: rng2)
    m3 = COOMatrix(Complex).rand(50, 50, nonzero_elements: 4)
    m3.nonzeros.should eq 4
    m1.should eq m2
    m1.should_not eq m3
    expect_raises(ArgumentError) { COOMatrix(Complex).rand(5, 5, nonzero_elements: 26) }
  end

  it "can be initialized with defined percentage of random values" do
    rng1 = Random.new(1)
    m1 = COOMatrix(Complex).rand(50, 50, fill_factor: 0.1, rng: rng1)
    m1.nonzeros.should be_close (50*50*0.1), 50*50*0.05
    rng2 = Random.new(1)
    m2 = COOMatrix(Complex).rand(50, 50, fill_factor: 0.1, rng: rng2)
    m3 = COOMatrix(Complex).rand(50, 50, fill_factor: 0.1)
    m3.nonzeros.should be_close (50*50*0.1), 50*50*0.05
    m1.should eq m2
    m1.should_not eq m3
  end

  it "support arithmetic operations" do
    ma = GMat[
      [1, 2, 0],
      [0, 0, 5],
      [0, 1, 0],
      [-2, -1, -3],
    ]
    mb = GMat[
      [0, -2, 0],
      [0, 1, -5],
      [0, -1, 1],
      [10, 0, 0],
    ]
    mas = COOMatrix(Float64).new(ma)
    mbs = COOMatrix(Float64).new(mb)

    (mas + mbs).should eq(ma + mb)
    (3*mas - 6*mbs).should eq(3*ma - 6*mb)
    mas.add!(-5, mbs)
    mas.should_not eq ma
    ma.add!(-5, mb)
    mas.should eq ma
  end

  it "support arithmetic operations with dense matrices" do
    ma = GMat[
      [1, 2, 0],
      [0, 0, 5],
      [0, 1, 0],
      [-2, -1, -3],
    ]
    mb = GMat[
      [0, -2, 0],
      [0, 1, -5],
      [0, -1, 1],
      [10, 0, 0],
    ]
    mas = COOMatrix(Float64).new(ma)
    mbs = COOMatrix(Float64).new(mb)

    (ma + mbs).should eq(ma + mb)
    (mas + mb).should eq(ma + mb)
    (ma - mbs).should eq(ma - mb)
    (mas - mb).should eq(ma - mb)

    banded = BMat.new(3, 4, 2, 1) { |i, j| (i + 1)*10 + j + 1 }
    expect_raises(Exception) { mas + banded }
    expect_raises(Exception) { banded + mas }
    expect_raises(Exception) { mas - banded }
    expect_raises(Exception) { banded - mas }
  end

  it "support #tril and #triu" do
    m = COOMatrix(Float64).new(10, 10)
    m[0, 0] = 10
    m[1, 1] = 15
    m[5, 2] = 11
    m[3, 6] = 12

    m.tril(2).should eq m.to_dense.tril(2)
    m.tril.should eq m.to_dense.tril
    m.triu(2).should eq m.to_dense.triu(2)
    m.triu.should eq m.to_dense.triu

    ma = m.triu(1)
    m.should_not eq ma
    m.triu!(1)
    m.should eq ma
  end

  it "support select!" do
    m1 = COOMatrix(Float64).new(10, 10)
    m1[0, 0] = 10
    m1[1, 1] = -15
    m1[5, 2] = 11
    m1[3, 6] = 12

    m2 = m1.clone
    m2[1, 1] = 0
    m1.should_not eq m2
    m1.select! { |v| v >= 0 }
    m1.should eq m2
  end

  it "support select_with_index!" do
    m1 = COOMatrix(Float64).new(10, 10)
    m1[0, 0] = 10
    m1[1, 1] = -15
    m1[5, 2] = 11
    m1[3, 6] = 12

    m2 = m1.clone
    m2[1, 1] = 0
    m2[0, 0] = 0
    m1.should_not eq m2
    m1.select_with_index! { |v, i, j| i != j }
    m1.should eq m2
  end

  it "support resize!" do
    ma = GMat[
      [1, 2, 0],
      [0, 0, 5],
      [0, 1, 0],
      [-2, -1, -3],
    ]
    mas = COOMatrix(Float64).new ma
    ma.resize!(3, 2)
    mas.resize!(3, 2)
    mas.should eq ma
    ma.resize!(10, 10)
    mas.resize!(10, 10)
    mas.should eq ma
  end
end
