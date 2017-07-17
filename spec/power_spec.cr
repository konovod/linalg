require "./spec_helper"

include LA
describe LA do
  it "can evaluate matrix integer power" do
    m = Mat.ones(3, 4)
    expect_raises(ArgumentError) do
      m**2
    end

    m = GMat[
      [1, 2, 0],
      [4, 5, 6],
      [7, 8, 9],
    ]
    (m**0).should eq Mat.eye(3)
    (m**1).should eq m
    (m**(-2)).should eq (m*m).inv
    (m**5).should eq m*m*m*m*m
  end
end
