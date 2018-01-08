using DataFrames
using DataStructures
using Logging

Base.length(it::DataStructures.SDMIncludeLast) = it.last - it.first + 1

"""
Merge columns from a dataframe to another one which is considered main.
The firts data column in each data frame is for arguments.
So do no use in mapping object the first column numver (it is argument column).
"""
function merge_data_frame_to_main_by_columns(main::DataFrame, extra::DataFrame;
    mapping_columns::Nullable{Dict{Int, Int}}=Nullable{Dict{Int, Int}}(),
    extra_boundaries::Nullable{Tuple{Float64, Float64}}=Nullable{Tuple{Float64, Float64}}(),
    main_extrapolation::Symbol=:zero, extra_extrapolation::Symbol=:zero, extra_transform::Function=y->y,
    smoothing_minpoints::Tuple{Int, Int}=(5, 5), smoothing_sharpness::Float64=-1.0)
  Logging.configure(level=INFO)
  info("Merging the selected extra table data to main ones. The mapping rules are defined as:\n$mapping_columns")
  info("Main table has size $(size(main)); Extra table has size $(size(extra))")
  info("Checking the input argument...")
  mapping = Dict{Int, Int}()
  if (isnull(mapping_columns) || isempty(get(mapping_columns))) && size(extra, 2) == size(main, 2)
    for i = 1:size(main, 2)-1
      mapping[i] = i
    end
    info("Empty mapping is defined and the main column quantity equals to the extra one. So we create one to one mapping rule by defaut.")
  elseif isnull(mapping_columns) && size(extra, 2) ≠ size(main, 2)
    error("No mapping is not specified for unequal column quantities in 'main'($(size(main, 2)-1)) and 'columns'($(size(extra, 2)-1))")
  elseif !isnull(mapping_columns)
    mapping = get(mapping_columns)
    @assert all(colindex -> 1 <= colindex <= size(extra, 2)-1, keys(mapping))
      prerequisite("A columns index in mapping must be not greater than the column dataframe column size", ("indices", keys(mapping)));
    @assert all(mainindex -> 1 <= mainindex <= size(main, 2)-1, values(mapping))
      prerequisite("A main index in mapping must be not greater than the main dataframe column size", ("indices", values(mapping)));
  end
  info("Done...")
  # ----

  main_mapping = map(reverse, mapping)
  info("The maping rules are reversed: $main_mapping")
  extra_b = isnull(extra_boundaries) ? (extra[1,1], extra[end,1]) : get(extra_boundaries)
  info("Going to make merging for extra data in boundaries [$(extra_b[1]), $(extra_b[2])]")

  info("Making spline functions out of main and extra data tables...")
  FC = convert_functionDataFrame_to_Functions(extra; extrapolation=extra_extrapolation)
  Fᵉ = map(fcf->begin x->extra_transform(fcf(x)) end, FC)
  Fᵐ = convert_functionDataFrame_to_Functions(main; extrapolation=main_extrapolation)
  info("Testing the data functions:")
  info("Testing Fᵐ:")
  for i=1:length(Fᵐ)
    fᵐ = Fᵐ[i]
    xarg = (main[1,1] + main[end,1]) / 2
    info("Fᵐ[$i]($xarg) = $(fᵐ(xarg))")
  end
  info("Testing Fᵉ:")
  for i=1:length(Fᵉ)
    fᵉ = Fᵉ[i]
    xarg = (extra_b[1] + extra_b[2]) / 2
    info("Fᵉ[$i]($xarg) = $(fᵉ(xarg))")
  end
  info("Done...")
  info("Making a function array mixture according to the inverted mapping rules $(main_mapping)...")
  Lᵐ = length(Fᵐ)
  Fᵐᵉ = copy(Fᵐ)
  for i = 1:length(Fᵐ)
    if haskey(main_mapping, i)
      iᶜ = main_mapping[i]
      Fᵐᵉ[i] = Fᵉ[iᶜ]
      info("Mapping rule $i → $iᶜ; the extra spline function is substituted.")
    end
  end
  info("Testing Fᵐᵉ:")
  for i=1:length(Fᵐᵉ)
    fᵐᵉ = Fᵐᵉ[i]
    xarg = (extra_b[1] + extra_b[end]) / 2
    info("Fᵐᵉ[$i]($xarg) = $(fᵐᵉ(xarg))")
  end
  info("Done...")

  function_dict = Dict(:main=>Fᵐ, :extra=>Fᵉ)
  info("The most complicated part... Making a merged argument grid in order to reach maximum precision and smoothing...")
  intervals = calculate_intervals_with_smoothing(
    [
      :main=>collect(main[:,1]),
      :extra=>(isnull(extra_boundaries) ? collect(extra[:,1]) : collect(extra[extra_b[1] .<= extra[:,1] .<= extra_b[2], 1]))
    ];
    smoothing_minpoints=smoothing_minpoints
  )
  info("The following intervals are evaluated for a good new grid datatable:")
  for interval in intervals
    info("Interval: type - $(interval.first); table name - $(interval.second[1]); range - [$(interval.second[2][1]), $(interval.second[2][end])]")
  end

  info("Making new merged data table...")
  newmain = createDataFrameWithColumns(main.colindex.names; tp=Float64)
  for interval in intervals
    info("Working with interval: type - $(interval.first); table name - $(interval.second[1]); range - [$(interval.second[2][1]), $(interval.second[2][end])] ...")
    interval_type = interval.first
    related_intervals = interval.second[1]
    X = interval.second[2]

    @assert interval_type == :plain || interval_type == :smoothing
    if interval_type == :plain
      @assert length(related_intervals) >= 1
    elseif interval_type == :smoothing
      if length(related_intervals) ≠ 2
        warn("length(related_intervals) ≠ 2, related_intervals=$related_intervals")
      end
    end

    info("Selecting functions array...")
    funcs = begin
      if interval_type == :plain
        if length(related_intervals) == 1 && related_intervals[1] == :main
          info("Main data table functions are selected")
          Fᵐ
        elseif (length(related_intervals) == 1 && related_intervals[1] == :extra) ||
          (length(related_intervals) == 2 && :extra ∈ related_intervals)
          info("Merged data table functions are selected")
          Fᵐᵉ
        end
      elseif interval_type == :smoothing
        Fˢ = copy(Fᵐ)
        info("Smoothed data table functions are going to be evaluated...")
        for i = 1:length(Fˢ)
          if haskey(main_mapping, i) && length(related_intervals) == 2
            iᶜ = main_mapping[i]
            info("Making a smoothing function for the mapping rule $i → $iᶜ")
            xˡ = X[1]; xʳ = X[end]
            info("The smoothing interval is [$xˡ, $xʳ]")
            m_e = related_intervals[1] == :main && related_intervals[2] == :extra
            info("The sequence of smoothing intervals: $(related_intervals[1]), $(related_intervals[2]); $(m_e ? "main → extra" : "extra → main")")
            fˡ, fʳ = m_e ? (Fᵐ[i](xʳ), Fᵉ[iᶜ](xˡ)) : (Fᵉ[iᶜ](xʳ), Fᵐ[i](xˡ))

            if fˡ === NaN || fʳ === NaN
              warn("fˡ=$fˡ, fʳ=$fʳ, m_e=$m_e; xˡ=$xˡ, xʳ=$xʳ")
            end

            σf = sign(fʳ - fˡ)
            αˢʰᵃʳᵖ = if smoothing_sharpness < 0 10 else smoothing_sharpness end
            x₀ = (X[1] + X[end])/2; α = σf * abs(X[end]-X[1]) / αˢʰᵃʳᵖ
            info("Sigmoid parameters: σf=$σf, x₀=$x₀, α=$α; (σf=$σf)")
            smoothing = if σf > 0
                if m_e
                  info("σf > 0, Fᵐ[$i], Fᵉ[$iᶜ]")
                  sigmoid_of_name(Fᵐ[i], Fᵉ[iᶜ], x₀, α, "Fᵐ[$i], Fᵉ[$iᶜ]")
                else
                  info("σf > 0, Fᵉ[$iᶜ], Fᵐ[$i]")
                  sigmoid_of_name(Fᵉ[iᶜ], Fᵐ[i], x₀, α, "Fᵉ[$iᶜ], Fᵐ[$i]")
                end
              else
                if m_e
                  info("σf < 0, Fᵉ[$iᶜ], Fᵐ[$i]")
                  sigmoid_of_name(Fᵉ[iᶜ], Fᵐ[i], x₀, α, "Fᵉ[$iᶜ], Fᵐ[$i]")
                else
                  info("σf < 0, Fᵐ[$iᶜ], Fᵉ[$i]")
                  sigmoid_of_name(Fᵐ[i], Fᵉ[iᶜ], x₀, α, "Fᵐ[$i], Fᵉ[$iᶜ]")
                end
              end
            info("Smooting function for the rule $i → $iᶜ is done; σ($xˡ)=$(smoothing(xˡ)) ↔∿↔ σ($xʳ)=$(smoothing(xʳ))")
            Fˢ[i] = smoothing
          else
            if length(related_intervals) ≠ 2
              warn("No smoothing for the related intervals $related_intervals for function $i; using main function data...")
            end
          end
        end
        info("Smoothing functions done: $Fˢ")
        Fˢ
      else
        error("The following interval types are supported: $(:plain), $(:smoothing)")
      end
    end

    info("Filling the new data table in the interval [$(X[1]), $(X[end])]")
    info("Funcs array: $funcs")
    for x in X
      row = zeros(Float64, Lᵐ+1)
      row[1] = x
      for i = 1:Lᵐ
        row[i+1] = funcs[i](x)
      end
      push!(newmain, row)
    end
    info("Done...")

  end
  info("Done!")
  return newmain
end


# -----------------------------------------------------------------------------

type Edge
  value::Float64
  intervals::Dict{Symbol, Symbol}
  middles::Set{Symbol}
end

function calculate_intervals_with_smoothing(VX::Vector{Pair{Symbol, Vector{Float64}}}; smoothing_minpoints::Tuple{Int, Int}=(5, 5))
  Logging.configure(level=INFO)
  info("Calculation of intervals with smoothing...")
  info("The given interval names: $(map(it->it.first, VX))")
  info("Interval properties:")
  for ititem in VX
    info("Interval: name - $(ititem.first); range - [$(ititem.second[1]), $(ititem.second[end])]")
  end
  # ------
  # =====
  info("Checking the input data...")
  @assert !isempty(VX) && length(VX) > 1
  VX_denotions = map(vx -> vx.first, VX)
  @assert length(VX_denotions) == length(unique(VX_denotions)) prerequisite("Vectors symbols must be unique", ("symbols", VX_denotions))
  for (i, VXⁱ) in enumerate(VX)
    @assert !isempty(VXⁱ) && length(VXⁱ.second) > 1 prerequisite(
      "The $i-th interval $(VXⁱ.first) [$(VXⁱ.second[1]), $(VXⁱ.second[end]))] must have more than one element", ("length", length(VXⁱ.second)))
  end
  σVX = map(vx -> sign(vx.second[end] - vx.second[1]), VX)
  @assert abs(reduce((sign1, sign2) -> sign1 + sign2, σVX)) == length(VX) prerequisite(
    "All VX intervals must have the same non-zero sign.", ("σVX", σVX))
  for (i, VXⁱ) in enumerate(VX)
    @assert issorted(VXⁱ.second; rev = σVX[i] == -1) prerequisite(
      "The $i-th interval $(VXⁱ.first) [$(VXⁱ.second[1]), $(VXⁱ.second[end]))] must be sorted.", ("σVX[$i]", σVX[i]))
  end
  info("Done...")
  # =====
  σ = 1
  MX = Dict(VX)
  grid = unique(sort(collect(Base.flatten(map(vx -> vx.second, VX))); rev = σ == -1))
  info("Merged grid has length $(length(grid))")

  intervals = Dict(map(it->it.first=>(it.second[1], it.second[end]), VX))
  info("Intervals dictionary: $intervals")
  edges_raw = sort(collect(
      Base.flatten(
        map(
          it -> [(it.first, :begin, it.second[1]), (it.first, :end, it.second[end])],
          VX
        )
      )
    );
    lt = (e1, e2) -> e1[3] < e2[3], rev = σ == -1
  )
  edge_values = unique(sort(map(edge->edge[3], edges_raw); rev = σ == -1))
  info("Edge points: $edges_raw")
  info("Edge values: $edge_values")

  edges = Vector{Edge}()
  for edge_value in edge_values
    info("Making Edge object for the edge value $edge_value")
    iedges_by_intervals = collect(find(e->abs(e[3] - edge_value) < 1e-8, edges_raw))
    edges_by_intervals = getindex(edges_raw, iedges_by_intervals)
    info("Edges by intervals: indices - $iedges_by_intervals, edges - $edges_by_intervals")
    intervals_for_edge = Dict(map(e->e[1]=>e[2], edges_by_intervals))
    middles = Set{Symbol}(keys(filter((ik, iv)->iv[1] < edge_value < iv[2], intervals)))
    edge = Edge(edge_value, intervals_for_edge, middles)
    push!(edges, edge)
  end
  info("Edge objects: $edges")

  intervals_out = Vector{Pair{Symbol, Tuple{Vector{Symbol}, Vector{Float64}}}}()
  edge_previous = edges[1]; grid_previous = 1
  for ie = 1:length(edges)
    edge = edges[ie]
    info("EDGE: $edge")
    first_edge = all(edge_type->edge_type==:begin, values(edge.intervals)) && isempty(edge.middles)
    last_edge = all(edge_type->edge_type==:end, values(edge.intervals)) && isempty(edge.middles)
    if first_edge
      edge_previous = edge
      grid_previous = 1
    elseif last_edge
      Xˡ = grid[grid_previous:end]
      Xˡ_names = collect(
        keys(
          filter(
            (k, v) -> begin
              matches = v[1] <= Xˡ[1] < v[2] && v[1] < Xˡ[end] <= v[2]
              if matches
                info("The interval $k → $v matches the left grid interval [$(Xˡ[1]), $(Xˡ[end])].")
                true
              else
                ΔX = abs(Xˡ[end] - Xˡ[1]); Δiv = abs(v[2] - v[1])
                ϵ = (1 - min(ΔX, Δiv)/max(ΔX, Δiv))*100
                if ϵ < 5 && Xˡ[end] <= v[2]
                  warn("The interval $k → $v matches the left grid interval [$(Xˡ[1]), $(Xˡ[end])] partially; ϵ=$ϵ, Xˡ[end]=$(Xˡ[end]), v[2]=$(v[2]).")
                  true
                else
                  info("The interval $k → $v _does_not_ match the left grid interval [$(Xˡ[1]), $(Xˡ[end])]; ϵ=$ϵ, Xˡ[end]=$(Xˡ[end]), v[2]=$(v[2]).")
                  false
                end
              end
            end,
            intervals
          )
        )
      )

      info("LAST INTERVAL: [$(Xˡ[1]), $(Xˡ[end])], names: $(Xˡ_names)")

      push!(intervals_out, :plain=>(Xˡ_names, Xˡ))
    else
      info("Handling a middle edge: $(edge)")
      left_names = map(it->it.first, collect(filter((ik, iv)->iv==:end, edge.intervals)))
      right_names = map(it->it.first, collect(filter((ik, iv)->iv==:begin, edge.intervals)))
      info("On the left: $left_names")
      info("On the right: $right_names")
      middle_intervals = isempty(edge.middles) ? Dict{Symbol, Vector{Float64}}() : filter((ik, iv)->ik ∈ edge.middles, MX)
      info("Middles: $(keys(middle_intervals))")

      Δxˡ = isempty(left_names) ? 0 : maximum(
        map(
          it->begin
            Δxⁱᵗ = abs(it.second[end]-it.second[end-1])
            info("Collecting the left step interval for the interval [$(it.second[1]), $(it.second[end])] with $(length(it.second)) points, Δxⁱᵗ=$Δxⁱᵗ")
            Δxⁱᵗ
          end
        , collect(filter((ik, iv)->ik ∈ left_names, MX))))
      Δxʳ = isempty(right_names) ? 0 : maximum(
        map(
          it->begin
            Δxⁱᵗ = abs(it.second[1]-it.second[2])
            info("Collecting the right step interval for the interval [$(it.second[1]), $(it.second[end])] with $(length(it.second)) points, Δxⁱᵗ=$Δxⁱᵗ")
            Δxⁱᵗ
          end
        , collect(filter((ik, iv)->ik ∈ right_names, MX))))
      info("Δxˡ=$Δxˡ, Δxʳ=$Δxʳ")

      Δxᵐˡ = isempty(middle_intervals) ? 0 : maximum(map(it->begin
                            ixl = findlast(x->x < edge.value, it)
                            return abs(it[ixl+1] - it[ixl])
                          end,
                          values(middle_intervals)))
      Δxᵐʳ = isempty(middle_intervals) ? 0 : maximum(map(it->begin
                            ixl = findfirst(x->x > edge.value, it)
                            return abs(it[ixl] - it[ixl-1])
                          end,
                          values(middle_intervals)))
      info("Δxᵐˡ=$Δxᵐˡ, Δxᵐʳ=$Δxᵐʳ")
      Δxˢ = maximum([Δxˡ, Δxʳ, Δxᵐˡ, Δxᵐʳ])
      info("Δxˢ=$(Δxˢ)")

      ix = findfirst(x->abs(x-edge.value) < 1e-8, grid); @assert ix > 0
      xˡ = edge.value-Δxˢ
      info("Looking for xˡ=$xˡ (Δxˢ=$Δxˢ) in the grid [$(grid[1]), $(grid[end])]")
      ixˡ = findlast(x-> x < xˡ, grid)
      if ixˡ <= 0
        warn("ixˡ ≤ 0 → ixˡ=$ixˡ for edge $(edge.value) over the grid [$(grid[1]), $(grid[end])]; Assuming ixˡ=1...")
        ixˡ = 1
      end
      xʳ = edge.value+Δxˢ
      info("Looking for xʳ=$xʳ (Δxˢ=$Δxˢ) in the grid [$(grid[1]), $(grid[end])]")
      ixʳ = findfirst(x->x > xʳ, grid); @assert ixʳ > 0
      if ixʳ <= 0
        warn("ixʳ ≤ 0 → ixʳ=$ixʳ for edge $(edge.value) over the grid [$(grid[1]), $(grid[end])]; Assuming ixʳ=1...")
        ixʳ = 1
      end
      info("ix=$ix, ixˡ=$ixˡ, ixʳ=$ixʳ")

      pᵐⁱⁿ⁻ˡ = smoothing_minpoints[1]; pᵐⁱⁿ⁻ʳ = smoothing_minpoints[2]
      pˡ = abs(ix - ixˡ); pˡ = pˡ < pᵐⁱⁿ⁻ˡ ? pᵐⁱⁿ⁻ˡ : pˡ
      pʳ = abs(ixʳ - ix); pʳ = pʳ < pᵐⁱⁿ⁻ʳ ? pᵐⁱⁿ⁻ʳ : pʳ
      info("pˡ=$pˡ, pʳ=$pʳ")

      ixˡ = if ixˡ == 1 && pˡ > 0 ixˡ else ix - pˡ end
      ixʳ = if ixʳ == 1 && pʳ > 0 ixʳ elseif ixˡ == 1 && pˡ > 0 1 else ix + pʳ end
      info("ix=$ix, ixˡ=$ixˡ, ixʳ=$ixʳ; xˡ=$(grid[ixˡ]), xʳ=$(grid[ixʳ]), x=$(grid[ix])")

      Xˡ = if ixˡ == 1 && pˡ > 0 (Float64)[] else grid[grid_previous:ixˡ-1] end
      if isempty(Xˡ)
        warn("Left interval for the point x=$(grid[ix]) is empty - skipping it and a smoothing too...")
      else
        info("Looking for names for the interval [$(Xˡ[1]), $(Xˡ[end])] in the intervals dict: $intervals")
        Xˡ_names = collect(
          keys(
            filter(
              (k, v) -> begin
                matches = v[1] <= Xˡ[1] < v[2] && v[1] < Xˡ[end] <= v[2]
                if matches
                  info("The interval $k → $v matches the left grid interval [$(Xˡ[1]), $(Xˡ[end])].")
                  true
                else
                  ΔX = abs(Xˡ[end] - Xˡ[1]); Δiv = abs(v[2] - v[1])
                  ϵ = (1 - min(ΔX, Δiv)/max(ΔX, Δiv))*100
                  if ϵ < 5 && Xˡ[end] <= v[2]
                    warn("The interval $k → $v matches the left grid interval [$(Xˡ[1]), $(Xˡ[end])] partially; ϵ=$ϵ, Xˡ[end]=$(Xˡ[end]), v[2]=$(v[2]).")
                    true
                  else
                    info("The interval $k → $v _does_not_ match the left grid interval [$(Xˡ[1]), $(Xˡ[end])]; ϵ=$ϵ, Xˡ[end]=$(Xˡ[end]), v[2]=$(v[2]).")
                    false
                  end
                end
              end,
              intervals
            )
          )
        )
        info("LEFT INTERVAL: [$(Xˡ[1]), $(Xˡ[end])], names: $(Xˡ_names)")
        push!(intervals_out, :plain=>(unique(Xˡ_names), Xˡ))


        Xˢ = collect(linspace(grid[ixˡ], grid[ixʳ], 3*abs(ixʳ-ixˡ+1)))
        info("Smooting interval for the point $(grid[ix]) has $(length(Xˢ)) point(s).")
        Xˢ_names = Vector{Symbol}()
        append!(Xˢ_names, Xˡ_names)
        Xʳ_min = reduce((it1, it2)->begin
                  i1 = it1.second
                  i2 = it2.second
                  return abs(i1[end] - i1[1]) < abs(i2[end]-i2[1]) ? it1 : it2
                end, collect(filter((k, v)->v[1] < Xˢ[end] < v[2], intervals)))
        info("Right smoothing minimal interval: $(Xʳ_min)")
        Xʳ_name = Xʳ_min.first
        info("Smoothing names on the right: $(Xʳ_name)")
        append!(Xˢ_names, [Xʳ_name])
        info("SMOOTHING INTERVAL: [$(Xˢ[1]), $(Xˢ[end])], names: $(Xʳ_name)")
        push!(intervals_out, :smoothing=>(unique(Xˢ_names), Xˢ))
      end

      edge_previous = edge
      grid_previous = if ixˡ == 1 && pˡ > 0 ixʳ else ixʳ + 1 end
    end
  end
  info("Done ⌣")
  return intervals_out
end

# -----------------------------------------------------------------------------

function get_affected_interval_pair(interval::Tuple{Float64, Float64}, intervals::Dict{Symbol, Tuple{Float64, Float64}})
  Logging.configure(level=INFO)
  info("Evaluating closest interval pair for the interval values $interval")
  left = interval[1]; right = interval[2]
  left_diffs = map(it->(it.first, abs(it.second[2] - left)), intervals)
  right_diffs = map(it->(it.first, abs(right - it.second[1])), intervals)

  min_left = Nullable{Pair{Symbol, Tuple{Float64, Float64}}}()
  min_right = Nullable{Pair{Symbol, Tuple{Float64, Float64}}}()

  left_result = reduce((d1, d2) -> d1[2] < d2[2] ? d1 : d2, left_diffs)
  right_result = reduce((d1, d2) -> d1[2] < d2[2] ? d1 : d2, right_diffs)

  min_left = Nullable{Pair{Symbol, Tuple{Float64, Float64}}}(left_result.first=>intervals[left_result.first])
  min_right = Nullable{Pair{Symbol, Tuple{Float64, Float64}}}(right_result.first=>intervals[right_result.first])

  info("Done...")
  return min_left, min_right
end

function get_affected_interval_minimal(interval::Tuple{Float64, Float64}, intervals::Dict{Symbol, Tuple{Float64, Float64}})
  Logging.configure(level=INFO)
  info("Evaluating minimal interval for the interval values $interval")
  @assert interval[1] ≠ interval[2]
  entries = map(it->it[2], enumerate(intervals))
  min_interval = Nullable{Pair{Symbol, Tuple{Float64, Float64}}}()
  for entry in entries[2:end]
    name = entry.first
    b = entry.second
    if (b[1] <= interval[1] <= b[2] && b[1] <= interval[2] <= b[2]) || (b[2] <= interval[1] <= b[1] && b[2] <= interval[2] <= b[1])
      if isnull(min_interval)
        min_interval = Nullable{Pair{Symbol, Tuple{Float64, Float64}}}(entry)
      else
        mb = get(min_interval).second
        if abs(b[2] - b[1]) < abs(mb[2] - mb[1])
          min_interval = Nullable{Pair{Symbol, Tuple{Float64, Float64}}}(entry)
        end
      end
    end
  end
  info("Done...")
  return min_interval
end

function get_affected_intervals(interval::Tuple{Float64, Float64}, intervals::Dict{Symbol, Tuple{Float64, Float64}})
  @assert interval[1] ≠ interval[2]
  interval_names = Vector{Symbol}()
  entries = map(it->it[2], enumerate(intervals))
  for entry in entries
    name = entry.first
    b = entry.second
    if (b[1] <= interval[1] <= b[2] && b[1] <= interval[2] <= b[2]) || (b[2] <= interval[1] <= b[1] && b[2] <= interval[2] <= b[1])
      push!(interval_names, name)
    end
  end
  return interval_names
end

"""
Sort a given set of data frames and create a new data frame out of them.
"""
function merge_data_frames_sorted_sequencially(
  data::Vector{DataFrame}; points_to_smooth::Int = 5,
  piece_merging::Function=(df1::DataFrame, df2::DataFrame, x1::Float64, x2::Float64) ->
  begin createDataFrameWithColumns(df1.colindex.names; tp=typeof(df1[1,1])) end)
  Logging.configure(level=INFO)
  if isempty(data)
    info("The vector of data items is empty - nothing to do.")
    return DataFrame()
  elseif length(data) == 1
    info("The vector of data items has a single data frame simply - nothing to do.")
    return data[1]
  end
  data_map = Dict{Tuple{Float64, Float64}, DataFrame}()
  info("Merging the data pieces of size '$(length(data))'.")
  for data_piece in data
    # the first column is considered as an argument one
    Rᵃ = data_piece[1,1]
    Rᵇ = data_piece[end,1]
    @assert Rᵃ < Rᵇ prerequisite("Rᵃ < Rᵇ", ("Rᵃ", Rᵃ), ("Rᵇ", Rᵇ))
    info("A data piece is defined in the argument interval [$Rᵃ, $Rᵇ].")
    data_map[(Rᵃ, Rᵇ)] = data_piece
  end
  info("Data with the folowing intervals are going to be ordered and merged: $(sort(collect(keys(data_map))))")
  intervals = sort(collect(keys(data_map)))
  # ---
  for i = 1:length(intervals)-1
    interval_1 = intervals[i]; interval_2 = intervals[i+1]
    Rᵃ¹ = interval_1[1]; Rᵇ¹ = interval_1[2]
    Rᵃ² = interval_2[1]; Rᵇ² = interval_2[2]
    info("Checking data in intervals [$Rᵃ¹, $Rᵇ¹] and [$Rᵃ², $Rᵇ²]...")
    @assert Rᵃ¹ < Rᵇ² prerequisite("Rᵃ¹ < Rᵇ² for $i, $(i+1)", ("Rᵃ¹", Rᵃ¹), ("Rᵇ²", Rᵇ²))
    @assert Rᵇ¹ <= Rᵃ² prerequisite("Rᵇ¹ ≤ Rᵃ² for $i, $(i+1)", ("Rᵇ¹", Rᵇ¹), ("Rᵃ²", Rᵃ²))
    cols_1 = size(data_map[interval_1], 2); cols_2 = size(data_map[interval_2], 2)
    @assert cols_1 == cols_1 prerequisite("cols¹ == cols² for $i, $(i+1); Data are not mergeable.", ("cols¹", cols_1), ("cols²", cols_2))
    col_names_1 = data_map[interval_1].colindex.names; col_names_2 = data_map[interval_2].colindex.names
    @assert col_names_1 == col_names_2 prerequisite("col_names_1 == col_names_2 for $i, $(i+1)", ("col_names_1", col_names_1), ("col_names_2", col_names_2))
  end
  info("Seems that the data are mergeable. Let's try it...")
  # ---
  # ---
  data_type = typeof(data_map[intervals[1]][1,1]) # TODO: check all column data types before merge
  data_frame = createDataFrameWithColumns(data_map[intervals[1]].colindex.names; tp=data_type)
  info("Going to do a merge of $(length(data)) pieces...")
  for i = 1:length(intervals)-1
    interval_1 = intervals[i]; interval_2 = intervals[i+1]
    Rᵃ¹ = interval_1[1]; Rᵇ¹ = interval_1[2]
    Rᵃ² = interval_2[1]; Rᵇ² = interval_2[2]
    info("Merging data intervals [$Rᵃ¹, $Rᵇ¹] and [$Rᵃ², $Rᵇ²]...")
    data_1 = data_map[interval_1][1:end-points_to_smooth-1,:]
    info("First data block: $(size(data_1))")
    data_2 = data_map[interval_2][points_to_smooth+1:end,:]
    info("Second data block: $(size(data_2))")
    x1 = data_map[interval_1][end-points_to_smooth,1]
    x2 = data_map[interval_2][points_to_smooth,1]
    info("A smoothing area is [$x1, $x2].")
    data_12 = piece_merging(data_map[interval_1], data_map[interval_2], x1, x2)
    append!(data_frame, data_1)
    append!(data_frame, data_12)
    append!(data_frame, data_2)
    info("Done merging [$Rᵃ¹, $Rᵇ¹] and [$Rᵃ², $Rᵇ²].")
  end
  info("Done merging of $(length(data)) pieces.")
  return data_frame
end

"""
We assume that the first column in data frame object contains data points.
"""
function merge_function_table_piece_to_main_one(main::DataFrame, piece::DataFrame;
    argframe::Tuple{Float64, Float64}=(-1.0, -1.0),
    piece_ycolumns::Vector{Int}=Vector{Int}(1:size(piece, 2)))
  Nₘ = size(main, 1); L = size(main, 2) - 1
  Nₚ = size(piece, 1); Lₚ = size(piece, 2) - 1
  Xᵖ = convert(Vector{Float64}, piece[:,1])
  Xᵐ = convert(Vector{Float64}, main[:,1])
  xₐ₁ = argframe[1]; xₐ₂ = argframe[2]
  x₁ₚ = first(Xᵖ); x₂ₚ = last(Xᵖ)
  x₁ₘ = first(Xᵐ); x₂ₘ = last(Xᵐ)
  @assert Lₚ == L prerequisite("Lₚ == L", ("Lₚ", Lₚ), ("L", L))
  @assert x₁ₚ <= xₐ₁ && xₐ₂ <= x₂ₚ prerequisite("x₁ₚ <= xₐ₁ && xₐ₂ <= x₂ₚ", ("x₁ₚ", x₁ₚ), ("xₐ₁", xₐ₁), ("x₂ₚ", x₂ₚ), ("xₐ₂", xₐ₂))
  @assert x₁ₘ <= x₁ₚ && x₂ₚ <= x₂ₘ prerequisite("x₁ₘ <= x₁ₚ && x₂ₚ <= x₂ₘ", ("x₁ₘ", x₁ₘ), ("x₁ₚ", x₁ₚ), ("x₂ₚ", x₂ₚ), ("x₂ₘ", x₂ₘ))

  @assert length(piece_ycolumns) <= L prerequisite("length(piece_ycolumns) <= L", ("length(piece_ycolumns)", length(piece_ycolumns)), ("L", L))
  for col in piece_ycolumns
    @assert 1 <= col <= L prerequisite("1 ≤ col ≤ L in piece_ycolumns", ("col", col), ("L", L))
  end
  pycols = isempty(piece_ycolumns) ? Vector{Int}(1:1:L) : Vector{Int}(piece_ycolumns)

  ixₚˡᵉᶠᵗ = findlast(x -> x < xₐ₁, Xᵖ) + 1; ixₚʳⁱᵍʰᵗ = findfirst(x -> x > xₐ₂, Xᵖ) - 1
  xₚˡᵉᶠᵗ = Xᵖ[ixₚˡᵉᶠᵗ]; xₚʳⁱᵍʰᵗ = Xᵖ[ixₚʳⁱᵍʰᵗ]
  ixₘˡᵉᶠᵗ = findlast(x -> x < xₚˡᵉᶠᵗ, Xᵐ); ixₘʳⁱᵍʰᵗ = findfirst(x -> x > xₚʳⁱᵍʰᵗ, Xᵐ)
  xₘˡᵉᶠᵗ = Xᵐ[ixₘˡᵉᶠᵗ]; xₘʳⁱᵍʰᵗ = Xᵐ[ixₘʳⁱᵍʰᵗ]
  @assert xₘˡᵉᶠᵗ < xₚˡᵉᶠᵗ prerequisite("xₘˡᵉᶠᵗ < xₚˡᵉᶠᵗ", ("xₘˡᵉᶠᵗ", xₘˡᵉᶠᵗ), ("xₚˡᵉᶠᵗ", xₚˡᵉᶠᵗ))
  @assert xₚʳⁱᵍʰᵗ < xₘʳⁱᵍʰᵗ prerequisite("xₚʳⁱᵍʰᵗ < xₘʳⁱᵍʰᵗ", ("xₚʳⁱᵍʰᵗ", xₚʳⁱᵍʰᵗ), ("xₘʳⁱᵍʰᵗ", xₘʳⁱᵍʰᵗ))

  Fᵐ = convert_functionDataFrame_to_Functions(main)
  Fᵖ = convert_functionDataFrame_to_Functions(piece; skip=Vector{Int}(filter(ix -> ix ∉ pycols, 1:1:L)))

  newdata = createDataFrameWithColumns(main.colindex.names)
  for k = 1:ixₘˡᵉᶠᵗ
    row = data_row_to_vector(main, k)
    push!(newdata, row)
  end

  # --- boundary area
  Δxₘˡᵉᶠᵗ = (xₚˡᵉᶠᵗ - xₘˡᵉᶠᵗ) / 10.0
  for x = xₘˡᵉᶠᵗ+Δxₘˡᵉᶠᵗ:Δxₘˡᵉᶠᵗ:xₚˡᵉᶠᵗ-Δxₘˡᵉᶠᵗ
    row = Vector{Float64}(L+1); row[1] = x
    for l = 1:L
      row[l+1] = l ∈ pycols ? (Fᵐ[l](x) + Fᵖ[l](x))/2.0 : Fᵐ[l](x)
    end
    push!(newdata, row)
  end
  # ---

  for kₚ = ixₚˡᵉᶠᵗ:ixₚʳⁱᵍʰᵗ
    x = Xᵖ[kₚ]
    row = Vector{Float64}(L+1); row[1] = x
    for lₚ = 1:Lₚ
      row[lₚ+1] = lₚ ∈ pycols ? Fᵖ[lₚ](x) : Fᵐ[lₚ](x)
    end
    push!(newdata, row)
  end

  # --- boundary area
  Δxₘʳⁱᵍʰᵗ = (xₘʳⁱᵍʰᵗ - xₚʳⁱᵍʰᵗ) / 10.0
  for x = xₚʳⁱᵍʰᵗ+Δxₘʳⁱᵍʰᵗ:Δxₘʳⁱᵍʰᵗ:xₘʳⁱᵍʰᵗ-Δxₘʳⁱᵍʰᵗ
    row = Vector{Float64}(L+1); row[1] = x
    for l = 1:L
      row[l+1] = l ∈ pycols ? (Fᵐ[l](x) + Fᵖ[l](x))/2.0 : Fᵐ[l](x)
    end
    push!(newdata, row)
  end
  # ---

  for k = ixₘʳⁱᵍʰᵗ:Nₘ
    row = data_row_to_vector(main, k)
    push!(newdata, row)
  end

  return newdata
end
