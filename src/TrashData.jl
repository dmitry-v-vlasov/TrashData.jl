module TrashData

include("util/util.jl")
include("matrix-data/triangle.jl")
include("dataframe/dataframe.jl")
include("interpolation/functions.jl")
include("interpolation/spline.jl")
include("interpolation/smoothing.jl")
include("table-data/merge.jl")
include("table-data/resample.jl")
include("io/file.jl")
include("io/utils.jl")

export ResampleStrategy

export data_size_symmetric_matrix, data_column_upper_matrix, data_column_lower_matrix
export data_symmetric_lowerTriangle_to_upperTriangle
export file_data_transform
export load_data, save_data, load_multiple_data
export round_data!, round_data_file

export createDataFrameWithColumns

export mirror_antisymmetric_data

export sigmoid_of, sigmoid

# ---
export convert_functionDataFrame_to_Functions
export merge_data_piece, merge_data_frames_sorted, join_data_files_sequencially
export merge_data_file_columns_to_main

export resample_data_frame, resample_data_file

export create_table_with_constants, create_table_with_constants!

export mpos, mvec

end # module
