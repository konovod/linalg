require "./spec_helper"

include LA
include LA::Sparse
describe COOMatrix do
  m = GMat.eye(5)
  ms = COOMatrix(Float64).new(m)
  pp! m, ms
end
