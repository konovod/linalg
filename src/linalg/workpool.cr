module LA
  # Work arrays pool for lapack routines
  # It isn't thread safe for now because crystal isn't multithreaded
  class WorkPool
    @area = Bytes.new(1)
    @used = 0

    def get(n) : Bytes
      # reallocate(n + @used)
      @area[@used, n].tap { @used += n }
    end

    def get_f32(n) : Slice(Float32)
      get(n*sizeof(Float32)).unsafe_as(Slice(Float32))
    end

    def get_f64(n) : Slice(Float64)
      get(n*sizeof(Float64)).unsafe_as(Slice(Float64))
    end

    def get_cmplx(n) : Slice(LibCBLAS::ComplexDouble)
      get(n*sizeof(LibCBLAS::ComplexDouble)).unsafe_as(Slice(LibCBLAS::ComplexDouble))
    end

    def get_i32(n) : Slice(Int32)
      get(n*sizeof(Int32)).unsafe_as(Slice(Int32))
    end

    def release
      return if @used == 0
      # TODO - make it debug only
      raise "worksize guard failed" unless @area[@used] == 0xDE &&
                                           @area[@used + 1] == 0xAD &&
                                           @area[@used + 2] == 0xBE &&
                                           @area[@used + 3] == 0xEF
      @used = 0
    end

    def reallocate(required_size)
      required_size += 4
      n = @area.size
      if n < required_size
        while n < required_size
          n = n*2
        end
        @area = Bytes.new(n)
      end
      @area[required_size - 4] = 0xDE
      @area[required_size - 3] = 0xAD
      @area[required_size - 2] = 0xBE
      @area[required_size - 1] = 0xEF
    end
  end

  WORK_POOL = WorkPool.new
end
