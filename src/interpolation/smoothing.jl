function smoothing_within_interval(X::Vector{Float64}, Y::Function, interval::Tuple{Float64, Float64})
  @assert !isempty(X)
  @assert !isempty(interval)
  @assert length(interval) == 2
  @assert interval[1] ≠ interval[2]

  σˣ = sign(X[end] - X[1])
  @assert issorted(X; rev = σˣ < 0)

  σⁱ = sign(interval[2] - interval[1])

  _X = if σˣ > 0 X else reverse(X) end
  _Y = if σˣ > 0 Y else x->Y(x - X[end]) end
  _interval = if σⁱ > 0 interval else reverse(interval) end
end
