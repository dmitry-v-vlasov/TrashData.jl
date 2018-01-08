const TYPES_INT = Set([Int8, Int16, Int32, Int64, Int128])
const TYPES_UINT = Set([UInt8, UInt16, UInt32, UInt64, UInt128])
const ℤ = TYPES_INT ∪ TYPES_UINT
const ℝ = Set([Float16, Float32, Float64])

function assert_square_matrix_indices(i::Int, j::Int, N::Int)
  @assert i > 0 && j > 0 prerequisite("i > 0 & j > 0", ("i", i), ("j", j))
  @assert N > 0 prerequisite("N > 0", ("N", N))
  @assert i <= N && j <= N prerequisite("i ≤ N & j ≤ N", ("i", i), ("j", j), ("N", N))
end

function assert_elements_in_set(elements::Set, set::Set)
  if isempty(elements) return end
  @assert length(elements) <= length(set) prerequisite("|{elements}| ≤ |{set}|", ("elements", elements), ("set", set))
  for element in elements
    @assert element ∈ set prerequisite("element ∈ {set}", ("element", element), ("set", set))
  end
end

function prerequisite(text::AbstractString, pars::Tuple{Any, Any}...)
  details =  if length(pars) > 0
      "found: {$(join(map(p->"$(p[1])=$(p[2])", pars), ", "))}"
    else
      "please fix it."
    end
  "Prerequisite '$text', $details"
end

function exists(f::Function, collection::Set)
  return findfirst(f, collection) > 0
end

function exists(f::Function, collection::Vector)
  return findfirst(f, collection) > 0
end

function exists(f::Function, collection::Array)
  return findfirst(f, collection) > 0
end

function frexp10(value::Int)
  return frexp10(Float64(value))
end

function frexp10(value::Float64)
  v = value
  p = 0.0
  while v < 1.0
    v *= 10
    p -= 1.0
  end
  while v > 10.0
    v /= 10
    p += 1.0
  end
  return (v, p)
end

function weighted_mean_max(values::Vector{Float64})
  return weighted_mean(values; weightf=it->it[1])
end

function weighted_mean(values::Vector{Float64}; weightf::Function=it::Tuple{Int, Float64}->1.0)
  e_values = enumerate(values)
  return mapreduce(it->weightf(it)*it[2], +, e_values)/mapreduce(it->weightf(it), +, e_values)
end

function mean_base10(values::Float64...; meanf::Function = mean)
  vals = collect(values)
  return mean_base10(vals; meanf = meanf)
end

function mean_base10(values::Vector{Float64}; meanf::Function = mean)
  values_parts = frexp10.(values)
  s_values = map(v->v[1], values_parts); e_values = map(v->v[2], values_parts)
  sv = meanf(s_values); ev = meanf(e_values)
  return sv * 10^ev
end

function mean_base2(values::Float64...; meanf::Function = mean)
  vals = collect(values)
  return average_base2(vals; meanf = meanf)
end

function mean_base2(values::Vector{Float64}; meanf::Function = mean)
  s_values = map(significand, values); e_values = collect(Float64, map(exponent, values))
  sv = meanf(s_values); ev = meanf(e_values)
  return sv * 2^ev
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

"""
Convert vector element number to a pair of N×N matrix indices.
"""
function mpos(l::Int, N::Int)
  n, r = divrem(l, N)
  return (r == 0) ? n : n + 1, (r == 0) ? N : r
end

"""
Convert a pair of N×N matrix indices to a one-dimentional array index.
"""
function mvec(i::Int, j::Int, N::Int)
  return N*(i - 1) + j
end
