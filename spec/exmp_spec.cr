require "./spec_helper"

include Linalg
describe Linalg do
  it "can evaluate matrix exponent" do
    m = GMat[{0, 6, 0, 0}, {0, 0, 6, 0}, {0, 0, 0, 6}, {0, 0, 0, 0}]
    # m.detect MatrixFlags::Triangular
    # pp m.normAm(1)
    pp m.expm
    # pp m.expm(schur_fact: true)
    # pp (-m).expm
    # pp (m + m.t).expm
    mbad = GMat[[1, 1e9], [0, 1]]
    # mbad.detect MatrixFlags::Triangular
    pp mbad.expm
  end
end
