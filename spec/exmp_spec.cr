require "./spec_helper"

include Linalg
describe Linalg do
  it "can evaluate matrix exponent" do
    m = GMat[{0, 6, 0, 0}, {0, 0, 6, 0}, {0, 0, 0, 6}, {0, 0, 0, 0}]
    pp m.expm
  end
end
