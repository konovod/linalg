module LA
  module Matrix(T)
    def assume!(flag : MatrixFlags, value : Bool = true)
      if value
        @flags |= flag
      else
        @flags &= ~flag
      end
    end

    # TODO - check for all flags
    private def check_single(flag : MatrixFlags, eps = tolerance)
      case flag
      when .symmetric?
        # TODO - 2 times less work
        return false unless square?
        each_with_index do |value, row, column|
          return false if row < column && (value - unsafe_at(column, row)).abs > eps
        end
        return true
      when .hermitian?
        {% if T == Complex %}
          return false unless square?
          each_with_index do |value, row, column|
            return false if row < column && (value.conj - unsafe_at(column, row)).abs > eps
          end
          return true
        {% else %}
          return check_single(MatrixFlags::Symmetric, eps)
        {% end %}
      when .positive_definite?
        # TODO - cleaner detection?
        return false unless square? && check_single(MatrixFlags::Hermitian, eps)
        begin
          cholesky(dont_clean: true)
          return true
        rescue LinAlgError
          return false
        end
      when .orthogonal?
        return square? && (self*self.conjt).almost_eq Matrix(T).identity(nrows), eps
      when .upper_triangular?
        each_with_index do |value, row, column|
          return false if row > column && (value).abs > eps
        end
        return true
      when .lower_triangular?
        each_with_index do |value, row, column|
          return false if row < column && (value).abs > eps
        end
        return true
      else
        return false
      end
      return false
    end

    private def detect_single(flag : MatrixFlags, eps = tolerance)
      check_single(flag, eps).tap do |ok|
        assume!(flag, ok)
      end
    end

    def detect(aflags : MatrixFlags = MatrixFlags::All, eps = tolerance)
      result = true
      {MatrixFlags::Symmetric,
       MatrixFlags::Hermitian,
       MatrixFlags::PositiveDefinite,
       MatrixFlags::Orthogonal,
       MatrixFlags::LowerTriangular,
       MatrixFlags::UpperTriangular}.each do |f|
        if aflags & f != MatrixFlags::None
          result = false unless detect_single(f, eps)
        end
      end
      if aflags.triangular?
        result = false unless detect_single(MatrixFlags::LowerTriangular, eps) || detect_single(MatrixFlags::UpperTriangular, eps)
      end
      result
    end

    def assume(aflags : MatrixFlags, value : Bool = true)
      if value
        detect(aflags)
      else
        assume! aflags, false
      end
    end

    def clear_flags
      @flags = MatrixFlags::None
    end
  end
end
