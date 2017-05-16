import Dierckx
using DataFrames

const EXTRAPOLATION_MAP = Dict{Symbol, AbstractString}(
  :nearest=>"nearest",
  :zero=>"zero",
  :extrapolate=>"extrapolate",
  :error=>"error"
)

"""
We assume that the first column in data frame object contains data points.
"""
function convert_functionDataFrame_to_Functions(data::DataFrame; skip::Vector{Int}=Vector{Int}(), extrapolation::Symbol=:extrapolate)
  Logging.configure(level=INFO)
  @assert haskey(EXTRAPOLATION_MAP, extrapolation) "Unknown extrapolation type '$extrapolation'."
  N = size(data, 2) - 1
  functions = Vector{Function}(N)
  X = convert(Vector{Float64}, data[:, 1])
  info("Making a vector of spline functions with size $N.")
  for l = 1:N
    if (!isempty(skip) && l ∉ skip) || isempty(skip)
      Y = convert(Vector{Float64}, data[l+1])
      spline = Dierckx.Spline1D(X, Y;  w=ones(length(X)), k=3, bc=EXTRAPOLATION_MAP[extrapolation], s=0.0)
      func = x -> begin Dierckx.evaluate(spline, x) end
      functions[l] = func
    else
      functions[l] = x -> begin error("Column $l is skipped intentionally.") end
    end
  end
  return functions
end
