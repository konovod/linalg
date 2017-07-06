require "./spec_helper"

include Linalg
describe Linalg do
  it "can evaluate matrix exponent" do
    m = GMat[{0, 6, 0, 0}, {0, 0, 6, 0}, {0, 0, 0, 6}, {0, 0, 0, 0}]
    mres = GMat[
      [1.0, 6.0, 18.0, 36.0],
      [0.0, 1.0, 6.0, 18.0],
      [0.0, 0.0, 1.0, 6.0],
      [0.0, 0.0, 0.0, 1.0],
    ]
    m.expm.should almost_eq mres
    m.expm(schur_fact: true).should almost_eq mres
    m.detect MatrixFlags::Triangular
    m.expm.should almost_eq mres
    (-m).expm.should almost_eq GMat[
      [1.0, -6.0, 18.0, -36.0],
      [0.0, 1.0, -6.0, 18.0],
      [0.0, 0.0, 1.0, -6.0],
      [0.0, 0.0, 0.0, 1.0],
    ]
    (m + m.t).expm.should almost_eq GMat[
      [2288.37752914804, 3687.89823361033, 3669.66066846686, 2258.86846809108],
      [3687.89823361033, 5958.03819761490, 5946.76670170141, 3669.66066846686],
      [3669.66066846686, 5946.76670170141, 5958.03819761491, 3687.89823361034],
      [2258.86846809108, 3669.66066846686, 3687.89823361034, 2288.37752914805],
    ]
    mbad = GMat[[1, 1e9], [0, 1]]
    mbad.expm.should almost_eq GMat[[Math::E, Math::E*1e9], [0, Math::E]]
  end

  it "apply special optimization for expm of triangular matrices" do
    mbad = GMat[[1, 1e9], [0, 1]]
    mbad.detect MatrixFlags::Triangular
    mbad.expm.should almost_eq GMat[[Math::E, Math::E*1e9], [0, Math::E]]

    m = GMat[{0, 6, 0, 0}, {0, 0, 6, 0}, {0, 0, 0, 6}, {0, 0, 0, 0}]
    (m + m.t).expm(schur_fact: true).should almost_eq GMat[
      [2288.37752914804, 3687.89823361033, 3669.66066846686, 2258.86846809108],
      [3687.89823361033, 5958.03819761490, 5946.76670170141, 3669.66066846686],
      [3669.66066846686, 5946.76670170141, 5958.03819761491, 3687.89823361034],
      [2258.86846809108, 3669.66066846686, 3687.89823361034, 2288.37752914805],
    ]
  end
end
