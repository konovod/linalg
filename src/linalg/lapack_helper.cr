require "./libLAPACKE"

module LA
  module LapackHelper
    ARG_NORMAL = 0
    ARG_MATRIX = 1
    ARG_INTOUT = 2

    WORK_NONE           = 0
    WORK_DETECT         = 1
    WORK_DETECT_SPECIAL = 2
    WORK_EMPTY          = 3
    WORK_PARAM1         = 4
    WORK_PARAM2         = 5
  end

  abstract class Matrix(T)
    include LapackHelper

    private macro alloc_real_type(size)
      {% if T == Float32 %} WORK_POOL.get_f32({{size}}) {% else %} WORK_POOL.get_f64({{size}}) {% end %}
    end

    private macro alloc_type(size)
      {% if T == Complex %} WORK_POOL.get_cmplx({{size}}) {% elsif T == Float32 %} WORK_POOL.get_f32({{size}}) {% else %} WORK_POOL.get_f64({{size}}) {% end %}
    end

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
      {% if T == Complex
           name = name.stringify.gsub(/^(or)/, "un").id
         end %}
      info = LibLAPACKE.{{typ}}{{name}}(LibCBLAS::COL_MAJOR, {{*args}})
      raise LinAlgError.new("LAPACKE.{{typ}}{{name}} returned #{info}") if info != 0
    end

    macro lapack(name, *args, worksize = nil)

      {%
        lapack_args = {
          "gebal" => {3 => ARG_MATRIX, 5 => ARG_INTOUT, 6 => ARG_INTOUT, 7 => ARG_MATRIX},
          # "gees"  => {3 => ARG_MATRIX, 5 => ARG_MATRIX, 7 => ARG_INTOUT, 8 => ARG_MATRIX, 9 => ARG_MATRIX}, SPECIAL CASE
          # "geev"  => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 7 => ARG_MATRIX, 9 => ARG_MATRIX}, SPECIAL CASE
          "gehrd" => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "gels"  => {5 => ARG_MATRIX, 7 => ARG_MATRIX},
          "gelsd" => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 8 => ARG_MATRIX, 10 => ARG_INTOUT},
          "gelsy" => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 8 => ARG_MATRIX, 10 => ARG_INTOUT},
          "geqp3" => {3 => ARG_MATRIX, 5 => ARG_MATRIX, 6 => ARG_MATRIX},
          "geqrf" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "gerqf" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "gelqf" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "geqlf" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "gesdd" => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 7 => ARG_MATRIX, 9 => ARG_MATRIX},
          "gesv"  => {3 => ARG_MATRIX, 5 => ARG_MATRIX, 6 => ARG_MATRIX},
          "getrf" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "getri" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "getrs" => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 7 => ARG_MATRIX},
          # "gges" => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 8 => ARG_MATRIX, 10 => ARG_INTOUT, 11 => ARG_MATRIX, 12 => ARG_MATRIX, 13 => ARG_MATRIX, 15 => ARG_MATRIX}, SPECIAL CASE
          # "ggev" => {},SPECIAL CASE
          "heevr" => {5 => ARG_MATRIX, 12 => ARG_INTOUT, 13 => ARG_MATRIX, 14 => ARG_MATRIX, 16 => ARG_MATRIX},
          "hegvd" => {5 => ARG_MATRIX, 7 => ARG_MATRIX, 9 => ARG_MATRIX},
          "hesv"  => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 7 => ARG_MATRIX},
          "hetrf" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "hetri" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "orghr" => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "orgqr" => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "orgrq" => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "orglq" => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "orgql" => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "posv"  => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "potrf" => {3 => ARG_MATRIX},
          "potri" => {3 => ARG_MATRIX},
          "potrs" => {4 => ARG_MATRIX, 6 => ARG_MATRIX},
          "syevr" => {5 => ARG_MATRIX, 12 => ARG_INTOUT, 13 => ARG_MATRIX, 14 => ARG_MATRIX, 16 => ARG_MATRIX},
          "sygvd" => {5 => ARG_MATRIX, 7 => ARG_MATRIX, 9 => ARG_MATRIX},
          "sysv"  => {4 => ARG_MATRIX, 6 => ARG_MATRIX, 7 => ARG_MATRIX},
          "sytrf" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "sytri" => {3 => ARG_MATRIX, 5 => ARG_MATRIX},
          "trtri" => {4 => ARG_MATRIX},
          "trtrs" => {6 => ARG_MATRIX, 8 => ARG_MATRIX},
        }

        lapack_worksize = {
          "hetri" => {"cwork" => WORK_PARAM1},
          # "lantr" => {"rwork" => WORK_PARAM1},
          # "lanhe" => {"rwork" => WORK_PARAM1},
          # "lange" => {"rwork" => WORK_PARAM1},
          # "lansy" => {"rwork" => WORK_PARAM1},
        }
      %}


      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
      {% func_args = lapack_args[name.stringify] %}
      {% func_worksize = lapack_worksize[name.stringify] %}

      {% if func_worksize %}
        {% if func_worksize["cwork"] %}
          {% if func_worksize["cwork"] == WORK_PARAM1 %}
            %csize = {{worksize[0]}}
          {% elsif func_worksize["cwork"] == WORK_PARAM2 %}
            %csize = {{worksize[1]}}
          {% end %}
          %cbuf = alloc_type(%csize)
        {% end %}

        {% if func_worksize["rwork"] %}
          {% if func_worksize["rwork"] == WORK_PARAM1 %}
            %rsize = {{worksize[0]}}
          {% elsif func_worksize["rwork"] == WORK_PARAM2 %}
            %rsize = {{worksize[1]}}
          {% end %}
          %rbuf = alloc_real_type(%rsize)
        {% end %}

        {% if func_worksize["iwork"] %}
          {% if func_worksize["iwork"] == WORK_PARAM1 %}
            %isize = {{worksize[0]}}
          {% elsif func_worksize["iwork"] == WORK_PARAM2 %}
            %isize = {{worksize[1]}}
          {% end %}
          %ibuf = WORK_POOL.get_i32(%isize)
        {% end %}
      {% end %}


      {% if T == Complex
           name = name.stringify.gsub(/^(or)/, "un").id
         end %}

      {% for arg, index in args %}
        {% argtype = func_args[index + 1] %}
        {% if argtype == ARG_MATRIX %}
        {% elsif argtype == ARG_INTOUT %}
          {{arg}} = 0
        {% else %}
        %var{index} = {{arg}}
        {% end %}
      {% end %}

       info = 0
       LibLAPACK.{{typ}}{{name}}_(
         {% for arg, index in args %}
         {% argtype = func_args[index + 1] %}
         {% if argtype == ARG_MATRIX %}
           {{arg}},
         {% elsif argtype == ARG_INTOUT %}
           pointerof({{arg}}),
         {% else %}
          pointerof(%var{index}),
         {% end %}
         {% end %}

         {% if func_worksize %}
           {% if func_worksize["cwork"] %}
              %cbuf,
           {% end %}
         {% end %}

         pointerof(info))

         {% if func_worksize %}
           {% if func_worksize["cwork"] %}
              WORK_POOL.release(%cbuf)
           {% end %}
         {% end %}
      raise LinAlgError.new("LAPACK.{{typ}}{{name}} returned #{info}") if info != 0
    end
  end
end
