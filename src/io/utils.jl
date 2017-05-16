using Logging

function mirror_antisymmetric_data(in_file::AbstractString, out_file::AbstractString, N::Int; header=true)
  file_data_transform(in_file, out_file,
    data_symmetric_lowerTriangle_to_upperTriangle,
    Dict{Any, Any}("N"=>N);
    element_transform=x->-x, header=header)
end

function merge_data_piece(in_file_main::AbstractString, in_file_piece::AbstractString, out_file::AbstractString;
  argframe::Tuple{Float64, Float64}=(-1.0, -1.0),
  piece_ycolumns::Vector{Int}=Vector{Int}())
  main_data = load_data(in_file_main; header=true)
  piece_data = load_data(in_file_piece; header=true)
  out_data = merge_function_table_piece_to_main_one(main_data, piece_data;
              argframe=argframe,
              piece_ycolumns=piece_ycolumns)
  save_data(out_data, out_file; header=true)
end

function join_data_files_sequencially(files::Vector{String}, out_file::AbstractString; header=true)
  Logging.configure(level=INFO)
  info("Join the data files")
  data_items = load_multiple_data(files; header=header)
  merged_data =
    merge_data_frames_sorted_sequencially(
      data_items; points_to_smooth = 5,
      piece_merging = (df1::DataFrame, df2::DataFrame, x1::Float64, x2::Float64) ->
      begin
        Logging.configure(level=INFO)
        dtype = typeof(df1[1,1])
        result = createDataFrameWithColumns(df1.colindex.names; tp=dtype)

        dff1 = convert_functionDataFrame_to_Functions(df1; extrapolation=:nearest)
        @assert !isempty(dff1)
        dff2 = convert_functionDataFrame_to_Functions(df2; extrapolation=:nearest)
        @assert !isempty(dff2)
        funcs = Vector{Function}()
        for i = 1:length(dff1)
          σ = sign(dff2[i](x2) - dff1[i](x1))
          x₀ = (x2 + x1) / 2; α = σ * (x2 - x1) / e
          push!(funcs, sigmoid_of(dff1[i], dff2[i], x₀, α))
        end
        Δx = (x2 - x1) / 20
        info("Smoothing an area [$x1, $x2] with step Δx=$Δx.")
        for x = x1:Δx:x2
          row = Vector{Float64}(length(funcs) + 1)
          row[1] = x
          for i=2:length(row)
            row[i] = funcs[i-1](x)
          end
          push!(result, row)
        end
        return result
      end
    )
  info("Done joining")
  info("Savig the result to '$out_file'")
  save_data(merged_data, out_file; header=header)
  info("Done")
  return merged_data
end

function extend_data_in_interval(xᵃ::Float64, Δx::Float64, xᵇ::Float64, data::DataFrame; extrapolation::Symbol=:nearest)
  @assert xᵃ < xᵇ prerequisite("xᵃ < xᵇ", ("xᵃ", xᵃ), ("xᵇ", xᵇ))
  @assert Δx < xᵇ - xᵃ prerequisite("Δx < xᵇ - xᵃ", ("Δx", Δx), ("xᵇ - xᵃ", xᵇ - xᵃ))
  @assert extrapolation ≠ :error prerequisite("extrapolation ≠ :error", ("extrapolation", extrapolation))
  xmᵃ = data[1,1]; xmᵇ = data[end,1]
  if xmᵃ <= xᵃ < xᵇ <= xmᵇ
    info("The interval [$xᵃ, $xᵇ] is within the inverval [$xmᵃ, $xmᵇ] - no need to extend data.")
    return data
  end
  info("Extending the interval [$xmᵃ, $xmᵇ] with [$xᵃ, $xᵇ]...")

  info("Making spline...")
  funcs = convert_functionDataFrame_to_Functions(data; extrapolation=extrapolation)
  @assert length(funcs) == size(data, 2)
  info("Making spline... done")

  result_data = createDataFrameWithColumns(data.colindex.names; tp=Float64)
  if xᵃ < xmᵃ
    X = xmᵃ - xᵃ > Δx ? collect(xᵃ:Δx:xmᵃ) : [xᵃ]
    for x in X
      row = Vector{Float64}(length(funcs)+1)
      row[1] = x
      for i = 2:length(row)
        row[i] = funcs[i-1](x)
      end
      push!(result_data, row)
    end
  end
  append!(result_data, data)
  if xmᵇ < xᵇ
    X = xᵇ - xmᵇ > Δx ? collect(xmᵇ+Δx:Δx:xᵇ) : [xᵇ]
    if xᵇ ∉ X
      push!(X, xᵇ)
    end
    for x in X
      row = Vector{Float64}(length(funcs)+1)
      row[1] = x
      for i = 2:length(row)
        row[i] = funcs[i-1](x)
      end
      push!(result_data, row)
    end
  end
end
