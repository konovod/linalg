require "./spec_helper"
include Linalg

describe Linalg::Matrix do
  it "can be constructed from matlab string" do
    Mat.from_matlab("[1,2,3; 4,5,6; 7,8,9  ] ").should eq GMat[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
  end
  it "can be printed to matlab string" do
    Mat.eye(3).to_matlab.should eq "[1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0]"
    MatComplex.eye(3).to_matlab.should eq "[1.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i; 0.0 + 0.0i, 1.0 + 0.0i, 0.0 + 0.0i; 0.0 + 0.0i, 0.0 + 0.0i, 1.0 + 0.0i]"
  end
  it "keeps precision when printing" do
    m = Mat.rand(10, 10) * 1e-19
    Mat.from_matlab(m.to_matlab).should almost_eq m
    m = Mat.rand(10, 10) * 1e+19
    Mat.from_matlab(m.to_matlab).should almost_eq m
  end

  it "can be save to and loaded from csv" do
    m = Mat.rand(10, 10).vcat((Mat.rand(10, 10)*6).expm)
    m.save_csv "test.csv"
    m1 = Mat.load_csv "test.csv"
    m1.should almost_eq m
    File.delete "test.csv"
  end
end
