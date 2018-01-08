using DataFrames

function file_data_transform(in_file::AbstractString, out_file::AbstractString, transform::Function, pars::Dict{Any, Any}; element_transform::Function=x->x, header=true)
  in_data = load_data(in_file; header=header)
  out_data = transform(in_data, pars; transform=element_transform)
  save_data(out_data, out_file; header=header)
end

function load_matrix(file_name::AbstractString, N::Int, L::Int)
  if !isfile(file_name)
    return Nullable{Array{Float64, 2}}()
  end
  M = if N > 0 && L >0
      readdlm(file_name, ' ';
        header=false, skipstart=0, skipblanks=true,
        use_mmap=false, quotes=false, dims=(N, L),
        comments=true, comment_char='#')
    else
      readdlm(file_name, ' ';
        header=false, skipstart=0, skipblanks=true,
        use_mmap=false, quotes=false,
        comments=true, comment_char='#')
    end
  return Nullable{Array{Float64, 2}}(M)
end

function load_matrix(file_name::AbstractString, N::Int)
  return load_matrix(file_name, N, N)
end

function load_matrix(file_name::AbstractString)
  return load_matrix(file_name, -1, -1)
end
