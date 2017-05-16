# TrashData.jl

```Julia
using TrashData
mirror_antisymmetric_data("data.ddr_aq", "data.ddr_aq_t", 11; header=false)
```

```Julia
using TrashData
merge_data_piece("ddr-CaH.dsv", "cahp1.ddr_aq_t", "merged.dsv"; argframe=(29.6336, 29.6345), piece_ycolumns=[7,8,16,17,24,25,31,32,37,38,42,43,46,47,50,51,52,53,54])
```

For large data you've better use
```sh
split -C 150M data.ddr data.ddr_
```

```Julia
using TrashData
join_data_files_sequencially(["output-89-9.77927-transformation-matrix.dsv", "output-89-29.63403-transformation-matrix.dsv"], "result-t.dsv"; header=true);
```
