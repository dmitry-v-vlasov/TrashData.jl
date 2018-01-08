function sigmoid_of(f1::Function, f2::Function, x₀, α)
  return x -> begin
    Logging.configure(level=INFO)
    sf = (1 - sigmoid(x, x₀, α))*f1(x) + sigmoid(x, x₀, α)*f2(x)
    if sf === NaN || sf === NaN64 || sf === NaN32 || sf === NaN16
      error("SIGMOID NaN: x=$x, x₀=$x₀, α=$α, f1(x)=$(f1(x)), f2(x)=$(f2(x)), σ(x)=$(sigmoid(x, x₀, α))")
    end
    return sf
  end
end

function sigmoid_of_name(f1::Function, f2::Function, x₀, α, name::AbstractString)
  return x -> begin
    Logging.configure(level=INFO)
    sf = (1 - sigmoid(x, x₀, α))*f1(x) + sigmoid(x, x₀, α)*f2(x)
    if sf === NaN || sf === NaN64 || sf === NaN32 || sf === NaN16
      error("SIGMOID NaN: x=$x, x₀=$x₀, α=$α, f1(x)=$(f1(x)), f2(x)=$(f2(x)), σ(x)=$(sigmoid(x, x₀, α)), name=$name")
    end
    return sf
  end
end

const sigmoid = (x, x₀, α) -> 1/(1+exp(-(x - x₀)/α))
