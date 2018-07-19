require "./libLAPACKE"

module LA
  abstract class Matrix(T)
    macro lapack_util(name, worksize, *args)
      buf = alloc_real_type(worksize)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}

      {% for arg, index in args %}
        {% if !(arg.stringify =~ /^matrix\(.*\)$/) %}
          %var{index} = {{arg}}
        {% end %}
      {% end %}

      result = LibLAPACK.{{typ}}{{name}}_(
        {% for arg, index in args %}
          {% if !(arg.stringify =~ /^matrix\(.*\)$/) %}
            pointerof(%var{index}),
          {% else %}
            {{arg.stringify.gsub(/^matrix\((.*)\)$/, "(\\1)").id}},
          {% end %}
        {% end %}
        buf)
      WORK_POOL.release(buf)
      result
    end

    macro lapacke(storage, name, *args)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
      {% if T == Complex && storage.id == :or.id
           st = :un.id
         else
           st = storage
         end %}
      info = LibLAPACKE.{{typ}}{{st}}{{name}}(LibCBLAS::COL_MAJOR, {{*args}})
      raise LinAlgError.new("LAPACKE.{{typ}}{{storage}}{{name}} returned #{info}") if info != 0
    end

    macro lapack(storage, name, *args, **worksizes)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
      {% if T == Complex && storage.id == :or.id
           st = :un.id
         else
           st = storage
         end %}

       {% for arg, index in args %}
         {% if (arg.stringify =~ /^intout\(.*\)$/) %}
           {{arg.stringify.gsub(/^intout\((.*)\)$/, "\\1").id}} = 0
         {% elsif !(arg.stringify =~ /^matrix\(.*\)$/) %}
           %var{index} = {{arg}}
         {% else %}
         {% end %}
       {% end %}

       info = 0
       LibLAPACK.{{typ}}{{storage}}{{name}}_(
         {% for arg, index in args %}
           {% if (arg.stringify =~ /^intout\(.*\)$/) %}
             pointerof({{arg.stringify.gsub(/^intout\((.*)\)$/, "\\1").id}}),
           {% elsif !(arg.stringify =~ /^matrix\(.*\)$/) %}
             pointerof(%var{index}),
           {% else %}
             {{arg.stringify.gsub(/^matrix\((.*)\)$/, "(\\1)").id}},
           {% end %}
         {% end %}
         pointerof(info))

      raise LinAlgError.new("LAPACK.{{typ}}{{storage}}{{name}} returned #{info}") if info != 0
    end
  end
end
