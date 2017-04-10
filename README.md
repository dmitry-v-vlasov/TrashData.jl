# TrashData.jl

```Julia
using TrashData
mirror_antisymmetric_data("data.ddr_aq", "data.ddr_aq_t", 11; header=false)
```

For large data you've better use
```sh
split -C 150M data.ddr data.ddr_
```
