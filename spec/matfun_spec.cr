require "./spec_helper"

describe LA do
  it "evaluate expm, cosm, sinm, tanm" do
    m = LA::GMatComplex[[1.0, 2.0], [-1.0, 3.0]]
    LA.expm(1.i*m).should be_close (LA::GMatComplex[
      [0.42645930 + 1.89217551.i, -2.13721484 - 0.97811252.i],
      [1.06860742 + 0.48905626.i, -1.71075555 + 0.91406299.i]]), 1e-5

    (LA.cosm(m) + 1.i*LA.sinm(m)).should almost_eq LA.expm(1.i*m)

    LA.tanm(m).should almost_eq (LA.sinm(m)*LA.cosm(m).inv)
  end

  it "evalute sinhm, coshm, tanhm" do
    a = LA::GMat[[1.0, 3.0], [1.0, 4.0]]
    LA.sinhm(a).should be_close(LA::GMat[
      [10.57300653, 39.28826594],
      [13.09608865, 49.86127247]], 1e-6)

    LA.tanhm(a).should almost_eq (LA.sinhm(a)*LA.coshm(a).inv)
  end
end
