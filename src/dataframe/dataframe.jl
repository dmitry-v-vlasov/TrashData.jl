using DataFrames

function createDataFrameWithColumns(column_names::Vector{Symbol})
  df = DataFrame()
  for col_name in column_names
    df[col_name] = Any[]
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

function data_row_to_vector(data::DataFrame, irow::Int)
  L = size(data, 2)
  v = Vector{Float64}(L)
  for l = 1:L
    v[l] = data[irow, l]
  end
  return v
end
