using Logging

function resample_data_file(file_data::AbstractString, file_out::AbstractString; points_in_each_selection::Int=20, maximum_function_angle::Float64=π/8, keep_tail::Int=10, header=true)
  @assert isfile(file_data) prerequisite("Not a file or the file with name '$file_data' does not exist.")
  data = load_data(file_data; header=header)
  newdata = resample_data_frame(data, Set([:derivative_limit]),
      Dict(
        :points_in_each_selection=>points_in_each_selection,
        :maximum_function_angle=>maximum_function_angle,
        :keep_tail=>keep_tail
      )
    )
  save_data(newdata, file_out; header=header)
end

function round_data_file(file_data::AbstractString, file_out::AbstractString, roundings::Set{Tuple{Float64, Float64, Float64}}; header=true)
  @assert isfile(file_data) prerequisite("Not a file or the file with name '$file_data' does not exist.")
  @assert !isempty(roundings) prerequisite("Roundings can't be empty.")
  data = load_data(file_data; header=header)
  round_data!(data, roundings)
  save_data(data, file_out; header=header)
end

function round_data_file_columns(file_data::AbstractString, file_out::AbstractString,
    roundings::Set{Tuple{Float64, Float64, Float64}}; header=true)
  @assert isfile(file_data) prerequisite("Not a file or the file with name '$file_data' does not exist.")
  @assert !isempty(roundings) prerequisite("Roundings can't be empty.")
  data = load_data(file_data; header=header)
  round_data!(data, roundings=roundings)
  save_data(data, file_out; header=header)
end

function round_data_file_columns(file_data::AbstractString, file_out::AbstractString,
    roundings::Dict{Int, Set{Tuple{Float64, Float64, Float64}}}; header=true)
  @assert isfile(file_data) prerequisite("Not a file or the file with name '$file_data' does not exist.")
  @assert !isempty(roundings) prerequisite("Roundings can't be empty.")
  data = load_data(file_data; header=header)
  round_data!(data; roundings=roundings)
  save_data(data, file_out; header=header)
end

function create_table_with_constants(grid::Vector{Float64}, constants::Vector{Float64}; header=true, column_prefix::AbstractString="Y")
  @assert !isempty(grid) prerequisite("The provided grid array is empty.")
  @assert !isempty(constants) prerequisite("The provided constant array is empty.")

  column_names = Vector{Symbol}(length(constants)+1)
  column_names[1] = "X"
  for i = 1:length(constants)
    column_names[i+1] = "$(column_prefix)_$i"
  end
  data = createDataFrameWithColumns(column_names; tp=Float64)

  L = length(constants)
  for x in collect(grid[:,1])
    row = Vector{Float64}(L+1)
    row[1] = x
    for i = 1:L
      row[i+1] = constants[i]
    end
    push!(data, row)
  end

  return data
end

function create_table_with_constants!(file_grid::AbstractString, constants::Vector{Float64}, file_out::AbstractString; header=true, column_prefix::AbstractString="Y")
  @assert isfile(file_grid) prerequisite("The grid file does not exist", ("file_grid", file_grid))

  df_grid = load_data(file_grid; header=header)
  @assert !isempty(df_grid) prerequisite("The grid file must have some data")

  grid = collect(Float64, df_grid[:, 1])
  @assert !isempty(grid) prerequisite("The grid file must have some data")

  data = create_table_with_constants(grid, constants; header=header, column_prefix=column_prefix)

  save_data(data, file_out; header=true)
end

function create_table_with_constants!(
    file_grid::AbstractString, file_matrix_constants::AbstractString,
    file_out::AbstractString, matrix_size::Tuple{Int, Int}; header=true, column_prefix::AbstractString="Y")
  df_grid = load_data(file_grid; header=header)
  @assert !isempty(df_grid) prerequisite("The grid file must have some data")

  grid = collect(Float64, df_grid[:, 1])
  @assert !isempty(grid) prerequisite("The grid file must have some data")

  mat_constants = load_matrix(file_matrix_constants)
  @assert !isnull(mat_constants)
  constants = zeros(Float64, length(get(mat_constants)))
  N = size(get(mat_constants), 1); L = size(get(mat_constants), 2)
  @assert N == matrix_size[1] && L == matrix_size[2] prerequisite("Unexpected matrix size.", ("N", N), ("L", L), ("Nex", matrix_size[1]), ("Lex", matrix_size[2]))
  @assert N == L prerequisite("Only square matrices are supported: N=L", ("N", N), ("L", L))
  for i = 1:N, j=1:L
    l = mvec(i, j, N)
    constants[l] = get(mat_constants)[i, j]
  end

  data = create_table_with_constants(grid, constants; header=header, column_prefix=column_prefix)

  save_data(data, file_out; header=true)
end

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

function merge_data_file_columns_to_main(
  main_file::AbstractString, extra_file::AbstractString, out_file::AbstractString,
  mapping::Dict{Int, Int};
  extra_boundaries::Nullable{Tuple{Float64, Float64}}=Nullable{Tuple{Float64, Float64}}(),
  main_extrapolation::Symbol=:zero, extra_extrapolation::Symbol=:zero, extra_transform::Function=y->y,
  smoothing_minpoints::Tuple{Int, Int}=(5, 5), smoothing_sharpness::Float64=10.0,
  header=true)
  Logging.configure(level=INFO)
  info("Merging columns from '$extra_file' to '$main_file'")
  info("Column mapping rules: $mapping")
  main = load_data(main_file; header=header)
  extra = load_data(extra_file; header=header)
  out = merge_data_frame_to_main_by_columns(main, extra;
    mapping_columns=Nullable{Dict{Int, Int}}(mapping), extra_boundaries=extra_boundaries,
    main_extrapolation=main_extrapolation, extra_extrapolation=extra_extrapolation, extra_transform=extra_transform,
    smoothing_minpoints=smoothing_minpoints, smoothing_sharpness=smoothing_sharpness)
  save_data(out, out_file; header=header)
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
