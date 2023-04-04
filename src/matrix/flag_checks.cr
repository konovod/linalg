module LA
  abstract class Matrix(T)
    # Directly set or reset matrix `flag` without check
    # See `MatrixFlags` for description of flags
    def assume!(flag : MatrixFlags, value : Bool = true)
      if value
        self.flags |= flag
      else
        self.flags &= ~flag
      end
    end

    # TODO - check for all flags
    private def check_single(flag : MatrixFlags, eps = tolerance)
      case flag
      when .symmetric?
        return false unless square?
        each_upper(diagonal: false) do |value, row, column|
          return false if (value - unsafe_fetch(column, row)).abs > eps
        end
        return true
      when .hermitian?
        {% if T == Complex %}
          return false unless square?
          each_upper(diagonal: false) do |value, row, column|
            return false if (value.conj - unsafe_fetch(column, row)).abs > eps
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
        each_lower(diagonal: false) do |value, row, column|
          return false if value.abs > eps
        end
        return true
      when .lower_triangular?
        each_upper(diagonal: false) do |value, row, column|
          return false if value.abs > eps
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

    # Detect if given `aflags` are true or flase for a matrix with tolerance `eps`
    # Update `flags` property
    # See `MatrixFlags` for description of matrix flags
    # Returns True if _all_ given flags are set
    def detect?(aflags : MatrixFlags = MatrixFlags::All, eps = tolerance)
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

    # Detect if given `aflags` are true or flase for a matrix with tolerance `eps`
    # Update `flags` property
    # See `MatrixFlags` for description of matrix flags
    # Returns `self` for method chaining
    def detect(aflags : MatrixFlags = MatrixFlags::All, eps = tolerance)
      detect? aflags, eps
      self
    end

    # Reset matrix flags to None
    # Usually is done automatically,
    # but this method could be needed if internal content was changed using `to_unsafe`
    def clear_flags
      self.flags = MatrixFlags::None
    end
  end
end
