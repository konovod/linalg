require "./spec_helper"

include LA
describe LA::BandedMatrix do
  pending "is created from Mat##diag call" do
    GMat.diag([1, 2, 3]).should be_a BandedMatrix(Float64)
  end
  it "can be created by block" do
    a = BMat.new(5, 7, 2) { |i, j| i - j }
    a.should eq GMat[
      [0, -1, -2, 0, 0, 0, 0],
      [1, 0, -1, -2, 0, 0, 0],
      [2, 1, 0, -1, -2, 0, 0],
      [0, 2, 1, 0, -1, -2, 0],
      [0, 0, 2, 1, 0, -1, -2],
    ]
    a = BMat.new(7, 5, 2) { |i, j| i - j }
    a.should eq GMat[
      [0, -1, -2, 0, 0],
      [1, 0, -1, -2, 0],
      [2, 1, 0, -1, -2],
      [0, 2, 1, 0, -1],
      [0, 0, 2, 1, 0],
      [0, 0, 0, 2, 1],
      [0, 0, 0, 0, 2],
    ]
  end
  it "can be created by block with asymmetric bands" do
    a = BMat.new(5, 5, 1, 2) { |i, j| (i + 1)*10 + j + 1 }
    a.should eq GMat[
      [11, 12, 0, 0, 0],
      [21, 22, 23, 0, 0],
      [31, 32, 33, 34, 0],
      [0, 42, 43, 44, 45],
      [0, 0, 53, 54, 55],
    ]
  end
  it "can be created from provided diagonals" do
    BMat.new(5, 5, 1, 2,
      {
        [12, 23, 34, 45],
        [11, 22, 33, 44, 55],
        [21, 32, 43, 54],
        [31, 42, 53],
      }
    ).should eq GMat[
      [11, 12, 0, 0, 0],
      [21, 22, 23, 0, 0],
      [31, 32, 33, 34, 0],
      [0, 42, 43, 44, 45],
      [0, 0, 53, 54, 55],
    ]

    BMat.new(5, 7, 1, 2,
      {
        [12, 23, 34, 45, 56],
        [11, 22, 33, 44, 55],
        [21, 32, 43, 54],
        [31, 42, 53],
      }
    ).should eq GMat[
      [11, 12, 0, 0, 0, 0, 0],
      [21, 22, 23, 0, 0, 0, 0],
      [31, 32, 33, 34, 0, 0, 0],
      [0, 42, 43, 44, 45, 0, 0],
      [0, 0, 53, 54, 55, 56, 0],
    ]

    BMat.new(7, 5, 1, 2,
      {
        [12, 23, 34, 45],
        [11, 22, 33, 44, 55],
        [21, 32, 43, 54, 65],
        [31, 42, 53, 64, 75],
      }
    ).should eq GMat[
      [11, 12, 0, 0, 0],
      [21, 22, 23, 0, 0],
      [31, 32, 33, 34, 0],
      [0, 42, 43, 44, 45],
      [0, 0, 53, 54, 55],
      [0, 0, 0, 64, 65],
      [0, 0, 0, 0, 75],
    ]
    expect_raises(ArgumentError) do
      BMat.new(3, 4, 1, {[1, 2], [4, 5, 6], [7, 8]})
    end
    expect_raises(ArgumentError) do
      BMat.new(3, 4, 1, {[1, 2, 3], [4, 5, 6], [7, 8, 9]})
    end
  end

  it "can be created from general matrix" do
    a = GMat[
      [11, 12, 0, 0, 0],
      [21, 22, 23, 0, 0],
      [31, 32, 33, 34, 0],
      [0, 42, 43, 44, 45],
      [0, 0, 53, 54, 55],
      [0, 0, 0, 64, 65],
      [0, 0, 0, 0, 75],
    ]
    b = BMat.new(a)
    b.should be_a BandedMatrix(Float64)
    b.should eq a
  end

  it "can be created with given tolerance" do
    a = GMat[
      [11, 12, 0, 0, 0],
      [21, 22, 23, 0, 0],
      [31, 32, 33, 34, 0],
      [1e-6, 42, 43, 44, 45],
      [0, 0, 53, 54, 55],
      [0, 0, 0, 64, 65],
      [0, 0, 0, 0, 75],
    ]
    b = BMat.new(a)
    b.lower_band.should eq 3
    b.should eq a

    b2 = BMat.new(a, 1e-5)
    b2.lower_band.should eq 2
    b2.should_not eq a
  end

  it "can be created from band matrix of another type" do
    a = BMat32.new(5, 7, 2) { |i, j| i - j }
    a_complex = BMatComplex.new(a)
    a_complex.should eq GMatComplex[
      [0, -1, -2, 0, 0, 0, 0],
      [1, 0, -1, -2, 0, 0, 0],
      [2, 1, 0, -1, -2, 0, 0],
      [0, 2, 1, 0, -1, -2, 0],
      [0, 0, 2, 1, 0, -1, -2],
    ]
    a_float = BMat.new(a)
    a_float.should eq GMat[
      [0, -1, -2, 0, 0, 0, 0],
      [1, 0, -1, -2, 0, 0, 0],
      [2, 1, 0, -1, -2, 0, 0],
      [0, 2, 1, 0, -1, -2, 0],
      [0, 0, 2, 1, 0, -1, -2],
    ]
  end

  it "can be iterated with `each`" do
    a = BMat.new(3, 3, 1, {[1, 2], [3, 4, 5], [6, 7]})
    elements = [] of {Int32, Int32, Float64}
    a.each_with_index(all: true) do |v, i, j|
      elements << {i, j, v}
    end
    elements.should eq [
      {0, 0, 3.0},
      {0, 1, 1.0},
      {0, 2, 0.0},
      {1, 0, 6.0},
      {1, 1, 4.0},
      {1, 2, 2.0},
      {2, 0, 0.0},
      {2, 1, 7.0},
      {2, 2, 5.0},
    ]
    elements = [] of {Int32, Int32, Float64}
    a.each_with_index(all: false) do |v, i, j|
      elements << {i, j, v}
    end
    elements.should eq [
      {0, 1, 1.0},
      {1, 2, 2.0},
      {0, 0, 3.0},
      {1, 1, 4.0},
      {2, 2, 5.0},
      {1, 0, 6.0},
      {2, 1, 7.0},
    ]
  end

  it "can be mapped with map" do
    a = BMat.new(3, 4, 1, {[1, 2, 3], [4, 5, 6], [7, 8]})
    a1 = a.map { |v| -v }
    a1.should be_a BMat
    a1.should eq -a

    a.map_with_index! { |v, i, j| (i + 1)*10 + (j + 1) + v }
    a.should eq GMat[
      [11 + 4.0, 12 + 1.0, 0.0, 0.0],
      [21 + 7.0, 22 + 5.0, 23 + 2.0, 0.0],
      [0.0, 32 + 8.0, 33 + 6.0, 34 + 3.0],
    ]
  end

  it "equality check is optimized" do
    a = BMat.new(3, 4, 1, {[1, 2, 3], [4, 5, 6], [0, 0]})
    b = BMat.new(3, 4, 1, {[1, 2, 3], [4, 5, 6], [7, 8]})
    a2 = BMat.new(3, 4, 1, 0, {[1, 2, 3], [4, 5, 6]})
    c = BMat.new(3, 4, 0, 1, {[0, 0, 0], [7, 8]})

    a.should eq a2
    a.should_not eq b
    a.should_not eq c
    (a + c).should eq b
  end

  it "can be `transpose!`d inplace" do
    a = BMat.new(3, 4, 1, 2, {[1, 2, 3], [4, 5, 6], [7, 8], [9]})
    a.t!
    a.should eq GMat[
      [4.0, 7.0, 9.0],
      [1.0, 5.0, 8.0],
      [0.0, 2.0, 6.0],
      [0.0, 0.0, 3.0],
    ]
  end
end
