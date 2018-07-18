module LA
  # Work arrays pool for lapack routines
  # It isn't thread safe for now because crystal isn't multithreaded
  class WorkPool
    @area = Bytes.new(1024)

    def get(n) : Bytes
      reallocate(n)
      @area
    end

    def get_f32(n) : Slice(Float32)
      get(n*sizeof(Float32)).unsafe_as(Slice(Float32))
    end

    def get_f64(n) : Slice(Float64)
      get(n*sizeof(Float64)).unsafe_as(Slice(Float64))
    end

    def release(ptr)
    end

    def reallocate(required_size)
      n = @area.size
      while n < required_size
        n = n*2
      end
      @area = Bytes.new(n)
    end
  end

  WORK_POOL = WorkPool.new
end
