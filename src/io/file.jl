using DataFrames

function file_data_transform(in_file::AbstractString, out_file::AbstractString, transform::Function, pars::Dict{Any, Any}; element_transform::Function=x->x, header=true)
  in_data = load_data(in_file; header=header)
  out_data = transform(in_data, pars; transform=element_transform)
  save_data(out_data, out_file; header=header)
end
