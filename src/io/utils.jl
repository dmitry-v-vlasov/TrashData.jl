
function mirror_antisymmetric_data(in_file::AbstractString, out_file::AbstractString, N::Int; header=true)
  file_data_transform(in_file, out_file,
    data_symmetric_lowerTriangle_to_upperTriangle,
    Dict{Any, Any}("N"=>N);
    element_transform=x->-x, header=header)
end

function merge_data_piece(in_file_main::AbstractString, in_file_piece::AbstractString, out_file::AbstractString;
  argframe::Tuple{Float64, Float64}=(-1.0, -1.0),
  piece_ycolumns::Vector{Int}=Vector{Int}())
  main_data = load_data(in_file_main; header=true)
  piece_data = load_data(in_file_piece; header=true)
  out_data = merge_function_table_piece_to_main_one(main_data, piece_data;
              argframe=argframe,
              piece_ycolumns=piece_ycolumns)
  save_data(out_data, out_file; header=true)
end
