using DataFrames
using Logging

function data_symmetric_lowerTriangle_to_upperTriangle(in_data::DataFrame, pars::Dict{Any, Any}; transform::Function=x->x)
  @assert length(pars) == 1
  @assert haskey(pars, "N")
  @assert typeof(pars["N"]) == Int

  Logging.configure(level=INFO)

  N = convert(Int, pars["N"])
  @assert N > 0 prerequisite("N > 0", ("N", N))
  Nᶜ = data_size_symmetric_matrix(N)
  @assert Nᶜ == size(in_data, 2) || Nᶜ == size(in_data, 2) - 1 prerequisite("Nᶜ == size(in_data, 2){ - 1}?", ("Nᶜ", Nᶜ), ("size(in_data, 2)", size(in_data, 2)))
  argcolumn = Nᶜ == size(in_data, 2) - 1

  out_data = createDataFrameWithColumns(in_data.colindex.names)
  @assert (argcolumn ? size(out_data, 2) == Nᶜ + 1 : size(out_data, 2) == Nᶜ) prerequisite("Out data ize is wrong", ("Nᶜ", Nᶜ), ("out-size", size(out_data, 2)), ("arg-column", argcolumn))

  Nᵈᵃᵗᵃ = size(in_data, 1)
  for n = 1:Nᵈᵃᵗᵃ
    logb = IOBuffer()
    if argcolumn
      #print(logb, "arg⁺ →  ")
    end
    out_row = Vector{Float64}(argcolumn ? Nᶜ + 1 : Nᶜ)
    for i = 1:N
      Nᵐʳᵒʷ = i - 1
      for j = 1:Nᵐʳᵒʷ
        cˡ = data_column_lower_matrix(i, j)
        cᵘ = data_column_upper_matrix(j, i, N)
        if argcolumn
          out_row[1] = in_data[n, 1]
        end
        out_row[argcolumn ? cᵘ + 1 : cᵘ] = transform(in_data[n, argcolumn ? cˡ + 1 : cˡ])
        #print(logb, "($(cᵘ) ←  $(cˡ); $i, $j) ")
      end
    end
    push!(out_data, out_row)
    #info("transform trace ►$(takebuf_string(logb))█")
  end
  return out_data
end
