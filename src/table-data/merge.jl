using DataFrames
using Logging

"""
Sort a given set of data frames and create a new data frame out of them.
"""
function merge_data_frames_sorted_sequencially(
  data::Vector{DataFrame}; points_to_smooth::Int = 5,
  piece_merging::Function=(df1::DataFrame, df2::DataFrame, x1::Float64, x2::Float64) ->
  begin createDataFrameWithColumns(df1.colindex.names; tp=typeof(df1[1,1])) end)
  Logging.configure(level=INFO)
  if isempty(data)
    info("The vector of data items is empty - nothing to do.")
    return DataFrame()
  elseif length(data) == 1
    info("The vector of data items has a single data frame simply - nothing to do.")
    return data[1]
  end
  data_map = Dict{Tuple{Float64, Float64}, DataFrame}()
  info("Merging the data pieces of size '$(length(data))'.")
  for data_piece in data
    # the first column is considered as an argument one
    Rᵃ = data_piece[1,1]
    Rᵇ = data_piece[end,1]
    @assert Rᵃ < Rᵇ prerequisite("Rᵃ < Rᵇ", ("Rᵃ", Rᵃ), ("Rᵇ", Rᵇ))
    info("A data piece is defined in the argument interval [$Rᵃ, $Rᵇ].")
    data_map[(Rᵃ, Rᵇ)] = data_piece
  end
  info("Data with the folowing intervals are going to be ordered and merged: $(sort(collect(keys(data_map))))")
  intervals = sort(collect(keys(data_map)))
  # ---
  for i = 1:length(intervals)-1
    interval_1 = intervals[i]; interval_2 = intervals[i+1]
    Rᵃ¹ = interval_1[1]; Rᵇ¹ = interval_1[2]
    Rᵃ² = interval_2[1]; Rᵇ² = interval_2[2]
    info("Checking data in intervals [$Rᵃ¹, $Rᵇ¹] and [$Rᵃ², $Rᵇ²]...")
    @assert Rᵃ¹ < Rᵇ² prerequisite("Rᵃ¹ < Rᵇ² for $i, $(i+1)", ("Rᵃ¹", Rᵃ¹), ("Rᵇ²", Rᵇ²))
    @assert Rᵇ¹ <= Rᵃ² prerequisite("Rᵇ¹ ≤ Rᵃ² for $i, $(i+1)", ("Rᵇ¹", Rᵇ¹), ("Rᵃ²", Rᵃ²))
    cols_1 = size(data_map[interval_1], 2); cols_2 = size(data_map[interval_2], 2)
    @assert cols_1 == cols_1 prerequisite("cols¹ == cols² for $i, $(i+1); Data are not mergeable.", ("cols¹", cols_1), ("cols²", cols_2))
    col_names_1 = data_map[interval_1].colindex.names; col_names_2 = data_map[interval_2].colindex.names
    @assert col_names_1 == col_names_2 prerequisite("col_names_1 == col_names_2 for $i, $(i+1)", ("col_names_1", col_names_1), ("col_names_2", col_names_2))
  end
  info("Seems that the data are mergeable. Let's try it...")
  # ---
  # ---
  data_type = typeof(data_map[intervals[1]][1,1]) # TODO: check all column data types before merge
  data_frame = createDataFrameWithColumns(data_map[intervals[1]].colindex.names; tp=data_type)
  info("Going to do a merge of $(length(data)) pieces...")
  for i = 1:length(intervals)-1
    interval_1 = intervals[i]; interval_2 = intervals[i+1]
    Rᵃ¹ = interval_1[1]; Rᵇ¹ = interval_1[2]
    Rᵃ² = interval_2[1]; Rᵇ² = interval_2[2]
    info("Merging data intervals [$Rᵃ¹, $Rᵇ¹] and [$Rᵃ², $Rᵇ²]...")
    data_1 = data_map[interval_1][1:end-points_to_smooth-1,:]
    info("First data block: $(size(data_1))")
    data_2 = data_map[interval_2][points_to_smooth+1:end,:]
    info("Second data block: $(size(data_2))")
    x1 = data_map[interval_1][end-points_to_smooth,1]
    x2 = data_map[interval_2][points_to_smooth,1]
    info("A smoothing area is [$x1, $x2].")
    data_12 = piece_merging(data_map[interval_1], data_map[interval_2], x1, x2)
    append!(data_frame, data_1)
    append!(data_frame, data_12)
    append!(data_frame, data_2)
    info("Done merging [$Rᵃ¹, $Rᵇ¹] and [$Rᵃ², $Rᵇ²].")
  end
  info("Done merging of $(length(data)) pieces.")
  return data_frame
end

"""
We assume that the first column in data frame object contains data points.
"""
function merge_function_table_piece_to_main_one(main::DataFrame, piece::DataFrame;
    argframe::Tuple{Float64, Float64}=(-1.0, -1.0),
    piece_ycolumns::Vector{Int}=Vector{Int}(1:size(piece, 2)))
  Nₘ = size(main, 1); L = size(main, 2) - 1
  Nₚ = size(piece, 1); Lₚ = size(piece, 2) - 1
  Xᵖ = convert(Vector{Float64}, piece[:,1])
  Xᵐ = convert(Vector{Float64}, main[:,1])
  xₐ₁ = argframe[1]; xₐ₂ = argframe[2]
  x₁ₚ = first(Xᵖ); x₂ₚ = last(Xᵖ)
  x₁ₘ = first(Xᵐ); x₂ₘ = last(Xᵐ)
  @assert Lₚ == L prerequisite("Lₚ == L", ("Lₚ", Lₚ), ("L", L))
  @assert x₁ₚ <= xₐ₁ && xₐ₂ <= x₂ₚ prerequisite("x₁ₚ <= xₐ₁ && xₐ₂ <= x₂ₚ", ("x₁ₚ", x₁ₚ), ("xₐ₁", xₐ₁), ("x₂ₚ", x₂ₚ), ("xₐ₂", xₐ₂))
  @assert x₁ₘ <= x₁ₚ && x₂ₚ <= x₂ₘ prerequisite("x₁ₘ <= x₁ₚ && x₂ₚ <= x₂ₘ", ("x₁ₘ", x₁ₘ), ("x₁ₚ", x₁ₚ), ("x₂ₚ", x₂ₚ), ("x₂ₘ", x₂ₘ))

  @assert length(piece_ycolumns) <= L prerequisite("length(piece_ycolumns) <= L", ("length(piece_ycolumns)", length(piece_ycolumns)), ("L", L))
  for col in piece_ycolumns
    @assert 1 <= col <= L prerequisite("1 ≤ col ≤ L in piece_ycolumns", ("col", col), ("L", L))
  end
  pycols = isempty(piece_ycolumns) ? Vector{Int}(1:1:L) : Vector{Int}(piece_ycolumns)

  ixₚˡᵉᶠᵗ = findlast(x -> x < xₐ₁, Xᵖ) + 1; ixₚʳⁱᵍʰᵗ = findfirst(x -> x > xₐ₂, Xᵖ) - 1
  xₚˡᵉᶠᵗ = Xᵖ[ixₚˡᵉᶠᵗ]; xₚʳⁱᵍʰᵗ = Xᵖ[ixₚʳⁱᵍʰᵗ]
  ixₘˡᵉᶠᵗ = findlast(x -> x < xₚˡᵉᶠᵗ, Xᵐ); ixₘʳⁱᵍʰᵗ = findfirst(x -> x > xₚʳⁱᵍʰᵗ, Xᵐ)
  xₘˡᵉᶠᵗ = Xᵐ[ixₘˡᵉᶠᵗ]; xₘʳⁱᵍʰᵗ = Xᵐ[ixₘʳⁱᵍʰᵗ]
  @assert xₘˡᵉᶠᵗ < xₚˡᵉᶠᵗ prerequisite("xₘˡᵉᶠᵗ < xₚˡᵉᶠᵗ", ("xₘˡᵉᶠᵗ", xₘˡᵉᶠᵗ), ("xₚˡᵉᶠᵗ", xₚˡᵉᶠᵗ))
  @assert xₚʳⁱᵍʰᵗ < xₘʳⁱᵍʰᵗ prerequisite("xₚʳⁱᵍʰᵗ < xₘʳⁱᵍʰᵗ", ("xₚʳⁱᵍʰᵗ", xₚʳⁱᵍʰᵗ), ("xₘʳⁱᵍʰᵗ", xₘʳⁱᵍʰᵗ))

  Fᵐ = convert_functionDataFrame_to_Functions(main)
  Fᵖ = convert_functionDataFrame_to_Functions(piece; skip=Vector{Int}(filter(ix -> ix ∉ pycols, 1:1:L)))

  newdata = createDataFrameWithColumns(main.colindex.names)
  for k = 1:ixₘˡᵉᶠᵗ
    row = data_row_to_vector(main, k)
    push!(newdata, row)
  end

  # --- boundary area
  Δxₘˡᵉᶠᵗ = (xₚˡᵉᶠᵗ - xₘˡᵉᶠᵗ) / 10.0
  for x = xₘˡᵉᶠᵗ+Δxₘˡᵉᶠᵗ:Δxₘˡᵉᶠᵗ:xₚˡᵉᶠᵗ-Δxₘˡᵉᶠᵗ
    row = Vector{Float64}(L+1); row[1] = x
    for l = 1:L
      row[l+1] = l ∈ pycols ? (Fᵐ[l](x) + Fᵖ[l](x))/2.0 : Fᵐ[l](x)
    end
    push!(newdata, row)
  end
  # ---

  for kₚ = ixₚˡᵉᶠᵗ:ixₚʳⁱᵍʰᵗ
    x = Xᵖ[kₚ]
    row = Vector{Float64}(L+1); row[1] = x
    for lₚ = 1:Lₚ
      row[lₚ+1] = lₚ ∈ pycols ? Fᵖ[lₚ](x) : Fᵐ[lₚ](x)
    end
    push!(newdata, row)
  end

  # --- boundary area
  Δxₘʳⁱᵍʰᵗ = (xₘʳⁱᵍʰᵗ - xₚʳⁱᵍʰᵗ) / 10.0
  for x = xₚʳⁱᵍʰᵗ+Δxₘʳⁱᵍʰᵗ:Δxₘʳⁱᵍʰᵗ:xₘʳⁱᵍʰᵗ-Δxₘʳⁱᵍʰᵗ
    row = Vector{Float64}(L+1); row[1] = x
    for l = 1:L
      row[l+1] = l ∈ pycols ? (Fᵐ[l](x) + Fᵖ[l](x))/2.0 : Fᵐ[l](x)
    end
    push!(newdata, row)
  end
  # ---

  for k = ixₘʳⁱᵍʰᵗ:Nₘ
    row = data_row_to_vector(main, k)
    push!(newdata, row)
  end

  return newdata
end
