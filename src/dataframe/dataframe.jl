using DataFrames
using Logging

function createDataFrameWithColumns(column_names::Vector{Symbol}; tp::DataType=Any)
  df = DataFrame()
  for col_name in column_names
    df[col_name] = tp[]
  end
  return df
end

function load_data(file::AbstractString; header=true)
  readtable(
    file,
    header = header, separator = ' ',
    allowcomments = true, commentmark = '#',
    skipblanks = true, encoding = :utf8, normalizenames = false
  )
end

function save_data(data::DataFrame, file::AbstractString; header=true)
  writetable(file, data; separator=' ', quotemark=' ', header=header, nastring="EMPTY")
end

function load_multiple_data(files::Vector{String}; header=true)
  Logging.configure(level=INFO)
  info("Multiple data files loading: $files")
  data = Vector{DataFrame}()
  for file in files
    info("Loading data file '$file'...")
    data_piece = readtable(
      file,
      header = header, separator = ' ',
      allowcomments = true, commentmark = '#',
      skipblanks = true, encoding = :utf8, normalizenames = false
    )
    push!(data, data_piece)
    info("Done")
  end
  return data
end

function data_row_to_vector(data::DataFrame, irow::Int)
  L = size(data, 2)
  v = Vector{Float64}(L)
  for l = 1:L
    v[l] = data[irow, l]
  end
  return v
end

function data_column_to_vector(data::DataFrame, icol::Int)
  L = size(data, 1)
  v = Vector{Float64}(L)
  for l = 1:L
    v[l] = data[l, icol]
  end
  return v
end

function round_data!(data::DataFrame;
  roundings::Set{Tuple{Float64, Float64, Float64}}=Set{Tuple{Float64, Float64, Float64}}(),
  shift::Int=1)
  if isempty(data)
    return
  end
  @assert !isempty(roundings) prerequisite("Roundings can't be empty.")
  @assert 0 <= shift <= size(data, 2)-1 prerequisite("0 ≤ shift ≤ size(data, 2)-1", ("shift", shift), ("size(data, 2)-1", size(data, 2)-1))

  dict_roundings=Dict{Int, Tuple{Float64, Float64, Float64}}()
  N = size(data, 2)
  for i = shift:size(data)
    dict_roundings[i] = roundings
  end
  round_data!(data; roundings=dict_roundings, shift=shift)
end

function round_data!(data::DataFrame;
  roundings::Dict{Int, Set{Tuple{Float64, Float64, Float64}}}=Dict{Int, Set{Tuple{Float64, Float64, Float64}}}(),
  shift::Int=1)
  if isempty(data)
    return
  end
  @assert !isempty(roundings) prerequisite("Roundings can't be empty.")
  @assert maximum(keys(roundings)) <= size(data, 2) prerequisite("maximum(keys(roundings)) ≤ size(data, 2)", ("maximum(keys(roundings))", maximum(keys(roundings)), ("size(data, 2)", size(data, 2))))
  @assert 0 <= shift <= size(data, 2)-1 prerequisite("0 ≤ shift ≤ size(data, 2)-1", ("shift", shift), ("size(data, 2)-1", size(data, 2)-1))
  @assert shift <= maximum(keys(roundings)) prerequisite("shift ≤ maximum(keys(roundings))", ("shift", shift), ("maximum(keys(roundings))", maximum(keys(roundings))))

  N = size(data, 1); M = size(data, 2)
  for i=1:N, j=1+shift:M
    v = data[i, j]
    for rounding in roundings
      ι, ρˢ = rounding.first, rounding.second
      if ι == j
        for ρ in ρˢ
          υ, ϵᵃᵇˢ, τ = ρ[1], ρ[2], ρ[3]
          if abs(v - υ) <= ϵᵃᵇˢ
            data[i, j] = τ
          end
        end
      end
    end
  end
end
