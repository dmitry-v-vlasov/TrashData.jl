module TrashData

include("util/util.jl")
include("matrix-data/triangle.jl")
include("io/file.jl")
include("io/utils.jl")
include("dataframe/dataframe.jl")


export data_size_symmetric_matrix, data_column_upper_matrix, data_column_lower_matrix
export data_symmetric_lowerTriangle_to_upperTriangle
export file_data_transform
export load_data, save_data

export createDataFrameWithColumns

export transform_antisymmetric_data

end # module
