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
    b = BMat.from(a)
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
    b = BMat.from(a)
    b.lower_band.should eq 3
    b.should eq a

    b2 = BMat.from(a, 1e-5)
    b2.lower_band.should eq 2
    b2.should_not eq a
  end
end
