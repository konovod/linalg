require "spec"

module Spec
  struct AlmostEqualExpectation(T)
    def initialize(@expected_value : T)
    end

    def match(actual_value)
      actual_value.almost_eq @expected_value
    end

    def failure_message(actual_value)
      expected = @expected_value.inspect
      got = actual_value.inspect
      if expected == got
        expected += " : #{@expected_value.class}"
        got += " : #{actual_value.class}"
      end
      "Expected: #{expected}\n     got: #{got}"
    end

    def negative_failure_message(actual_value)
      "Expected: actual_value != #{@expected_value.inspect}\n     got: #{actual_value.inspect}"
    end
  end

  module Expectations
    def almost_eq(value)
      Spec::AlmostEqualExpectation.new value
    end
  end
end
