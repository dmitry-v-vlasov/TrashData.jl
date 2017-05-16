function sigmoid_of(f1::Function, f2::Function, x₀, α)
  return x -> (1 - sigmoid(x, x₀, α))*f1(x) + sigmoid(x, x₀, α)*f2(x)
end

const sigmoid = (x, x₀, α) -> 1/(1+exp(-(x - x₀)/α))
