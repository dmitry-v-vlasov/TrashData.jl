function assert_square_matrix_indices(i::Int, j::Int, N::Int)
  @assert i > 0 && j > 0 prerequisite("i > 0 & j > 0", ("i", i), ("j", j))
  @assert N > 0 prerequisite("N > 0", ("N", N))
  @assert i <= N && j <= N prerequisite("i ≤ N & j ≤ N", ("i", i), ("j", j), ("N", N))
end

function prerequisite(text::AbstractString, pars::Tuple{Any, Any}...)
  @assert length(pars) > 0 "Prerequisite parameter list mus t not be empty."
  "Prerequisite '$text', found: {$(join(map(p->"$(p[1])=$(p[2])", pars), ", "))}"
end

function data_column_upper_matrix(i::Int, j::Int, N::Int)
  assert_square_matrix_indices(i, j, N)
  @assert i < j prerequisite("i < j", ("i", i), ("j", j))
  return convert(Int, (2*(N - 1) - i) * (i - 1) / 2 + (j - 1))
end

function data_column_lower_matrix(i::Int, j::Int)
  assert_square_matrix_indices(i, j, typemax(Int))
  @assert i > j prerequisite("i > j", ("i", i), ("j", j))
  return convert(Int, convert(Int, (1 + (i - 2)) * (i - 2) / 2 + j))
end

function data_size_symmetric_matrix(N::Int)
  @assert N > 0 prerequisite("N > 0", ("N", N))
  convert(Int, N*(N - 1) / 2)
end
