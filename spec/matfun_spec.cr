require "./spec_helper"

describe Linalg do
  it "evaluate expm, cosm, sinm, tanm" do
    m = Linalg::GMatComplex[[1.0, 2.0], [-1.0, 3.0]]
    Linalg.expm(1.i*m).should be_close (Linalg::GMatComplex[
      [0.42645930 + 1.89217551.i, -2.13721484 - 0.97811252.i],
      [1.06860742 + 0.48905626.i, -1.71075555 + 0.91406299.i]]), 1e-5

    (Linalg.cosm(m) + 1.i*Linalg.sinm(m)).should almost_eq Linalg.expm(1.i*m)

    Linalg.tanm(m).should almost_eq (Linalg.sinm(m)*Linalg.cosm(m).inv)
  end

  it "evalute sinhm, coshm, tanhm" do
    a = Linalg::GMat[[1.0, 3.0], [1.0, 4.0]]
    Linalg.sinhm(a).should be_close(Linalg::GMat[
      [10.57300653, 39.28826594],
      [13.09608865, 49.86127247]], 1e-6)

    Linalg.tanhm(a).should almost_eq (Linalg.sinhm(a)*Linalg.coshm(a).inv)
  end
end
