require "./spec_helper"

include LA
describe LA::BandedMatrix do
  it "is created from Mat##diag call" do
    Mat.diag([1, 2, 3]).should be_a BandedMatrix(Float64)
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

  it "'s bands count can be changed" do
    a = BMat.new(3, 4, 1, 2, {[1, 2, 3], [4, 5, 6], [7, 8], [9]})
    a.lower_band = 1
    a.should eq GMat[
      [4.0, 1.0, 0.0, 0.0],
      [7.0, 5.0, 2.0, 0.0],
      [0.0, 8.0, 6.0, 3.0],
    ]
    a.lower_band.should eq 1
    a.lower_band = 2
    a[2, 0].should eq 0.0
    a.lower_band.should eq 2

    a.upper_band = 0
    a.flags.should eq MatrixFlags::LowerTriangular
    a.upper_band = 1
    a.flags.should eq MatrixFlags::None

    expect_raises(ArgumentError) { a.upper_band = -1 }

    a.set_bands(0, 0)
    a.flags.should eq (MatrixFlags::LowerTriangular | MatrixFlags::UpperTriangular)
  end

  it "supports tril! and triu!" do
    a = BMat.new(GMat[
      [1, 2, 3, 0],
      [4, 5, 6, 7],
      [0, 8, 9, 1],
      [0, 0, 2, 3],
      [0, 0, 0, 4],
    ])
    b = a.t

    a.tril!(1)
    a.should eq GMat[
      [1, 2, 0, 0],
      [4, 5, 6, 0],
      [0, 8, 9, 1],
      [0, 0, 2, 3],
      [0, 0, 0, 4],
    ]
    a.flags.should eq MatrixFlags::None
    b.triu!(-1)
    b.should eq a.t
    b.flags.should eq MatrixFlags::None

    a.tril!
    a.should eq GMat[
      [1, 0, 0, 0],
      [4, 5, 0, 0],
      [0, 8, 9, 0],
      [0, 0, 2, 3],
      [0, 0, 0, 4],
    ]
    a.flags.should eq MatrixFlags::LowerTriangular
    b.triu!
    b.should eq a.t
    b.flags.should eq MatrixFlags::UpperTriangular

    a.tril!(-1)
    a.should eq GMat[
      [0, 0, 0, 0],
      [4, 0, 0, 0],
      [0, 8, 0, 0],
      [0, 0, 2, 0],
      [0, 0, 0, 4],
    ]
    a.flags.should eq MatrixFlags::LowerTriangular
    b.triu!(1)
    b.should eq a.t
    b.flags.should eq MatrixFlags::UpperTriangular
  end

  it "supports to_real and to_imag" do
    a = BMatComplex.new(GMatComplex[
      [1 + 1.i, 2.i, 3 + 3.i, 0],
      [0, 0, 6, 7],
    ])
    areal = a.to_real
    areal.should be_a BandedMatrix(Float64)
    areal.upper_band.should eq a.upper_band
    areal.lower_band.should eq a.lower_band
    areal.flags.upper_triangular?.should be_true
    areal.should eq GMat[
      [1, 0, 3, 0],
      [0, 0, 6, 7],
    ]
    aimag = a.to_imag
    aimag.should be_a BandedMatrix(Float64)
    aimag.upper_band.should eq a.upper_band
    aimag.lower_band.should eq a.lower_band
    aimag.flags.upper_triangular?.should be_true
    aimag.should eq GMat[
      [1, 2, 3, 0],
      [0, 0, 0, 0],
    ]
  end

  it " + and - with other BandedMatrix produce BandedMatrix" do
    a = BMat.new(5, 6, 2, 1) { |i, j| (i + 1)*10 + j + 1 }
    b = BMat.new(5, 6, 1, 2) { |i, j| (i + 1)*10 + j + 1 }
    r1 = a + (-b)
    r2 = a - b
    r1.should eq r2
    r1.should be_a BMat
    r1.should eq GMat[
      [0.0, 0.0, 13.0, 0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0, 24.0, 0.0, 0.0],
      [-31.0, 0.0, 0.0, 0.0, 35.0, 0.0],
      [0.0, -42.0, 0.0, 0.0, 0.0, 46.0],
      [0.0, 0.0, -53.0, 0.0, 0.0, 0.0],
    ]
  end

  it "can be transposed" do
    a = BMatComplex.new(GMatComplex[
      [1 + 1.i, 2.i, 3 + 3.i, 0],
      [0, 0, 6, 7],
    ])
    at = a.t
    at.should be_a BMatComplex
    at.should eq GMatComplex[
      [1.0 + 1.i, 0.0 + 0.i],
      [0.0 + 2.i, 0.0 + 0.i],
      [3.0 + 3.i, 6.0 + 0.i],
      [0.0 + 0.i, 7.0 + 0.i],
    ]
    aconj = a.conjt
    aconj.should be_a BMatComplex
    aconj.should eq (BMatComplex.new(at.to_real) - BMatComplex.new(at.to_imag)*1.i)
  end

  it "multiplies to normal matrix" do
    g = GMat[
      [1, 1, 0],
      [2, 2, 1],
      [0, 2, 3],
    ]
    ((g.inv)*BMat.new(g)).should almost_eq Mat.identity(3)
  end
end
