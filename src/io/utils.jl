
function transform_antisymmetric_data(in_file::AbstractString, out_file::AbstractString, N::Int; header=true)
  file_data_transform(in_file, out_file,
    data_symmetric_lowerTriangle_to_upperTriangle,
    Dict{Any, Any}("N"=>N);
    element_transform=x->-x, header=header)
end
