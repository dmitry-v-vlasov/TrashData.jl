using DataFrames

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
