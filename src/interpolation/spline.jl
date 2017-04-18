import Dierckx
using DataFrames

"""
We assume that the first column in data frame object contains data points.
"""
function convert_functionDataFrame_to_Functions(data::DataFrame)
  return convert_functionDataFrame_to_Functions(data; skip=Vector{Int}())
end

function convert_functionDataFrame_to_Functions(data::DataFrame; skip::Vector{Int}=Vector{Int}())
  N = size(data, 2) - 1
  functions = Vector{Function}(N)
  X = convert(Vector{Float64}, data[:, 1])
  for l = 1:N
    if (!isempty(skip) && l ∉ skip) || isempty(skip)
      Y = convert(Vector{Float64}, data[l+1])
      spline = Dierckx.Spline1D(X, Y;  w=ones(length(X)), k=3, bc="extrapolate", s=0.0)
      func = x -> begin Dierckx.evaluate(spline, x) end
      functions[l] = func
    else
      functions[l] = x -> begin error("Column $l is skipped intentionally.") end
    end
  end
  return functions
end
