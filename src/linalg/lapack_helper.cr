require "./libLAPACKE"

module LA
  module LapackHelper
    ARG_NORMAL = 0
    ARG_MATRIX = 1
    ARG_INTOUT = 2
  end

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

    macro lapacke(name, *args)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
      {% if T == Complex && name.stringify =~ /^or/
           name = name.stringify.gsub(/^(or)/, "un").id
         end %}
      info = LibLAPACKE.{{typ}}{{name}}(LibCBLAS::COL_MAJOR, {{*args}})
      raise LinAlgError.new("LAPACKE.{{typ}}{{name}} returned #{info}") if info != 0
    end

    macro lapack(name, *args, **worksizes)

      {%
        lapack_funcs = {
          "gebal" => {3 => LapackHelper::ARG_MATRIX, 5 => LapackHelper::ARG_INTOUT, 6 => LapackHelper::ARG_INTOUT, 7 => LapackHelper::ARG_MATRIX},
          "gesv"  => {3 => LapackHelper::ARG_MATRIX, 5 => LapackHelper::ARG_MATRIX, 6 => LapackHelper::ARG_MATRIX},
          "getrf" => {3 => LapackHelper::ARG_MATRIX, 5 => LapackHelper::ARG_MATRIX},
          "getrs" => {4 => LapackHelper::ARG_MATRIX, 6 => LapackHelper::ARG_MATRIX, 7 => LapackHelper::ARG_MATRIX},
          "posv"  => {4 => LapackHelper::ARG_MATRIX, 6 => LapackHelper::ARG_MATRIX},
          "potrf" => {3 => LapackHelper::ARG_MATRIX},
          "potri" => {3 => LapackHelper::ARG_MATRIX},
          "potrs" => {4 => LapackHelper::ARG_MATRIX, 6 => LapackHelper::ARG_MATRIX},
          "trtri" => {4 => LapackHelper::ARG_MATRIX},
          "trtrs" => {6 => LapackHelper::ARG_MATRIX, 8 => LapackHelper::ARG_MATRIX},
        }
      %}


      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
      {% func_data = lapack_funcs[name.stringify] %}

      {% if T == Complex && name.stringify =~ /^or/
           name = name.stringify.gsub(/^(or)/, "un").id
         end %}

      {% for arg, index in args %}
        {% argtype = func_data[index + 1] %}
        {% if argtype == LapackHelper::ARG_MATRIX %}
        {% elsif argtype == LapackHelper::ARG_INTOUT %}
          {{arg}} = 0
        {% else %}
        %var{index} = {{arg}}
        {% end %}
      {% end %}

       info = 0
       LibLAPACK.{{typ}}{{name}}_(
         {% for arg, index in args %}
         {% argtype = func_data[index + 1] %}
         {% if argtype == LapackHelper::ARG_MATRIX %}
           {{arg}},
         {% elsif argtype == LapackHelper::ARG_INTOUT %}
           pointerof({{arg}}),
         {% else %}
          pointerof(%var{index}),
         {% end %}
         {% end %}
         pointerof(info))

      raise LinAlgError.new("LAPACK.{{typ}}{{name}} returned #{info}") if info != 0
    end
  end
end
