using DataFrames
using DataStructures
using Logging

using ProgressMeter


const RESAMPLE_STRATEGIES = Set([:derivative_limit, :step_fixed, :step_min_max, :step_min, :step_max])
"""
Make resampling of a given fuction table data using a stragegy.
Supported strategies:
  - :derivative_limit, parameters - :angle_max, :last_points_number
  - :step_fixed, parameters - :step_value
  - :step_min_max, parameters - :step_min_value, :step_max_value
  - :step_min, parameters - :step_min_value
  - :step_max, parameters - :step_max_value
"""
function resample_data_frame(data::DataFrame, strategies::Set{Symbol}, parameters::Dict{Symbol, Any})
  Logging.configure(level=INFO)
  info("Data resampling with strategies $strategies...")
  @assert !isempty(data) prerequisite("The input data must not be empty.")
  @assert size(data)[2] >= 2 prerequisite("The input data must have two columns at least.")
  @assert !isempty(strategies) prerequisite("The strategy set must not be empty.")
  @assert issubset(strategies, RESAMPLE_STRATEGIES)
  @assert (:derivative_limit ∈ strategies) && length(strategies) == 1 prerequisite("Only :derivative_limit is supported at the moment.")

  @assert (:points_in_each_selection ∈ keys(parameters)) prerequisite(":points_in_each_selection ∈ {parameters}", ("parameters", parameters))
  @assert (typeof(parameters[:points_in_each_selection]) ∈ ℤ && parameters[:points_in_each_selection] > 0) prerequisite(":points_in_each_selection ∈ ℕ", (":points_in_each_selection", parameters[:points_in_each_selection]))
  points_in_each_selection = convert(UInt, parameters[:points_in_each_selection])
  @assert points_in_each_selection > 10 prerequisite(":points_in_each_selection > 10", (":points_in_each_selection", points_in_each_selection))

  # @assert (:preffered_step ∈ keys(parameters)) prerequisite(":preffered_step ∈ {parameters}", ("parameters", parameters))
  # @assert (typeof(parameters[:preffered_step]) ∈ ℝ) prerequisite(":preffered_step ∈ ℝ", )
  # preffered_step = convert(Float64, parameters[:preffered_step])
  # @assert preffered_step > 0 prerequisite(":preffered_step > 0", (":preffered_step", preffered_step))
  # @assert preffered_step > 1e-9 prerequisite(":preffered_step > 10⁻⁹", (":preffered_step", preffered_step))

  @assert (:maximum_function_angle ∈ keys(parameters)) prerequisite(":maximum_function_angle ∈ {parameters}", ("parameters", parameters))
  @assert (typeof(parameters[:maximum_function_angle]) ∈ ℝ) prerequisite("parameters[:maximum_function_angle] ∈ ℝ")
  maximum_function_angle = convert(Float64, parameters[:maximum_function_angle])
  @assert 0 < maximum_function_angle < π/2 prerequisite("0 < :maximum_function_angle < π/2 ($(π/2))", (":maximum_function_angle", maximum_function_angle))

  @assert (:keep_tail ∈ keys(parameters)) prerequisite(":keep_tail ∈ {parameters}", ("parameters", parameters))
  @assert (typeof(parameters[:keep_tail]) ∈ ℤ && parameters[:keep_tail] >= 0) prerequisite("parameters[:keep_tail] ∈ ℤ && > 0", ("keep_tail", parameters[:keep_tail]))
  keep_tail = convert(UInt, parameters[:keep_tail])

  info("Resampling data of length $(size(data)[1])...")
  result = resample_data_frame_derivative_limit(data, points_in_each_selection, maximum_function_angle, keep_tail)
  info("Done!")
  return result
end

function resample_data_frame_derivative_limit(data::DataFrame, Nᵖ::UInt, αᵐᵃˣ::Float64, tail::UInt)
  Logging.configure(level=INFO)
  info("Starting resample_data_frame_derivative_limit procedure...")
  X = Vector(data[:, 1])
  @assert issorted(X)
  F = data[:, 2:end]
  ΔX = filter(Δx -> 1e-9 < Δx < 1e2, X[2:end] - X[1:end-1])
  Δxᵐⁱⁿ = minimum(ΔX); Δxᵐᵃˣ = maximum(ΔX)
  Δxᵐᵉᵃⁿ = mean_base10(Δxᵐⁱⁿ, Δxᵐᵃˣ; meanf=weighted_mean_max)
  info("The grid steps: Δxᵐⁱⁿ=$Δxᵐⁱⁿ, Δxᵐᵃˣ=$Δxᵐᵃˣ, Δxᵐᵉᵃⁿ=$Δxᵐᵉᵃⁿ")
  αᵐᵃˣ¹ = if π/4 <= αᵐᵃˣ  π/2
    αᵐᵃˣ / 2
  elseif π/6 <= αᵐᵃˣ < π/4
    αᵐᵃˣ / 3
  elseif π/8 <= αᵐᵃˣ < π/6
    αᵐᵃˣ / 4
  elseif π/24 <= αᵐᵃˣ < π/8
    αᵐᵃˣ / 5
  else
    αᵐᵃˣ / 6
  end
  αᵐᵃˣ² = if π/4 <= αᵐᵃˣ  π/2
    αᵐᵃˣ
  elseif π/6 <= αᵐᵃˣ < π/4
    golden * αᵐᵃˣ
  elseif π/8 <= αᵐᵃˣ < π/6
    2 * αᵐᵃˣ
  elseif π/24 <= αᵐᵃˣ < π/8
    3 * αᵐᵃˣ
  else
    4 * αᵐᵃˣ
  end
  info("Angles: αᵐᵃˣ=$αᵐᵃˣ ($(rad2deg(αᵐᵃˣ))°), αᵐᵃˣ¹=$αᵐᵃˣ¹ ($(rad2deg(αᵐᵃˣ¹))°), π - αᵐᵃˣ¹ = $(π - αᵐᵃˣ¹) ($(rad2deg(π - αᵐᵃˣ¹))°)")

  info("Making spline interpolation for better resampling...")
  Fs = convert_functionDataFrame_to_Functions(data)
  info("Done...")
  Nᶜ = size(data)[2] - 1
  newdata = createDataFrameWithColumns(data.colindex.names; tp=Float64)
  cᵖ = 1
  i = 1
  ixˡᵃˢᵗ = 1
  #p = Progress(length(X), 1)
  while i <= length(X) - tail
    #ProgressMeter.update!(p, i)
    if cᵖ >= Nᵖ
      i¹ = i - cᵖ + 1; iᵉⁿᵈ = i
      # --- differences in the p-th interval of length Nᵖ ---
      ΔXᵖ = X[i¹+1:iᵉⁿᵈ] - X[i¹:iᵉⁿᵈ-1]
      ΔFᵖ = zeros(Float64, length(ΔXᵖ), Nᶜ)
      for j ∈ 1:Nᶜ
        ΔFᵖ[:, j] = collect(Float64, F[i¹+1:iᵉⁿᵈ, j] - F[i¹:iᵉⁿᵈ-1, j])
      end
      # -----------------------------------------------------

      # --- derivatives in the p-th interval ----------------
      tanαᵖ = zeros(Float64, size(ΔFᵖ)[1], Nᶜ)
      for j ∈ 1:Nᶜ
        tanαᵖ[:, j] = ΔFᵖ[:, j] ./ ΔXᵖ
      end
      αᵖ = atan(tanαᵖ)
      # -----------------------------------------------------

      # --- angles and their deltas in the p-th interval ----
      Δαᵖ = αᵖ[2:end, :] - αᵖ[1:end-1, :]
      ΣΔαᵖ = ones(Float64, 1, size(Δαᵖ)[1]) * Δαᵖ
      ΣΔαᵖ⁺ = ones(Float64, 1, size(Δαᵖ)[1]) * abs(Δαᵖ)
      ΔsignΣΔαᵖ = ones(Float64, 1, size(Δαᵖ)[1]) * sign(Δαᵖ)
      # -----------------------------------------------------

      if exists(ΣΔαᵖⱼ -> ΣΔαᵖⱼ >= αᵐᵃˣ, ΣΔαᵖ)
        for iᵖ = (i-Nᵖ+1):i
          push!(newdata, data_row_to_vector(data, convert(Int, iᵖ)))
        end
        ixˡᵃˢᵗ = i
      elseif exists(ΣΔαᵖ⁺ⱼ -> π/2 - αᵐᵃˣ <= ΣΔαᵖ⁺ⱼ <= π/2 + αᵐᵃˣ || π - αᵐᵃˣ¹ <= ΣΔαᵖ⁺ⱼ <= π + αᵐᵃˣ¹, ΣΔαᵖ⁺) && exists(ΔsignΣΔαᵖⱼ -> abs(ΔsignΣΔαᵖⱼ) <= Nᵖ/2+1, ΔsignΣΔαᵖ)
        info("Possible peak in [$(X[i-Nᵖ+1]), $(X[i])]; ΔsignΣΔαᵖ=$ΔsignΣΔαᵖ; max(Δαᵖ⁺)=$(maximum(ΣΔαᵖ⁺)) ($(rad2deg(maximum(ΣΔαᵖ⁺)))°)")
        for iᵖ = (i-Nᵖ+1):i
          push!(newdata, data_row_to_vector(data, convert(Int, iᵖ)))
        end
        ixˡᵃˢᵗ = i
      elseif exists(αᵖᵢⱼ -> abs(αᵖᵢⱼ) >= αᵐᵃˣ², αᵖ)
        xˡᵃˢᵗ = data[i-Nᵖ+1, 1]
        xᶜᵘʳʳᵉⁿᵗ = data[i, 1]
        xᵟ = (xˡᵃˢᵗ + xᶜᵘʳʳᵉⁿᵗ) / 2
        newrow = Vector{Float64}(Nᶜ + 1)
        newrow[1] = xᵟ
        for j ∈ 2:(Nᶜ+1)
          newrow[j] = Fs[j - 1](xᵟ)
        end
        push!(newdata, newrow)
        push!(newdata, data_row_to_vector(data, i))
      else
        push!(newdata, data_row_to_vector(data, i))
      end
      cᵖ = 0
    end
    cᵖ += 1;
    i += 1
  end
  for i ∈ (length(X) - tail + 1):length(X)
    push!(newdata, data_row_to_vector(data, convert(Int, i)))
  end
  info("The result has sizes $(size(newdata))")

  return newdata
end
