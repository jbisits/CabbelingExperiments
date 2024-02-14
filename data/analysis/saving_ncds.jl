# Script to practise saving to nc dataset

using NCDatasets

## Create a blank data set
NCDataset(joinpath(@__DIR__, "test_saving.nc"), "c") do ds

    defDim(ds, "time", 100)
    defDim(ds, "cumulative_volume", 1000)
    defVar(ds, "sortedC", Float64, ("cumulative_volume", "time"))

end

## See what is saved, all values are just the `fill value`
ds = NCDataset(joinpath(@__DIR__, "test_saving.nc"))
ds[:sortedC][:, 1]
close(ds)

## Append data to the saved variable
NCDataset(joinpath(@__DIR__, "test_saving.nc"), "a") do ds
    ds[:sortedC][:, 1] = 1.0:1000.0
end

## Check what is now saved, it is the values entered above
ds = NCDataset(joinpath(@__DIR__, "test_saving.nc"))
ds[:sortedC][:, 1]
close(ds)

## Fill the rest of the columns
NCDataset(joinpath(@__DIR__, "test_saving.nc"), "a") do ds

    for t ∈ 2:100
        ds[:sortedC][:, t] = 1.0:1000.0
    end

end

## check data, has worked!
ds = NCDataset(joinpath(@__DIR__, "test_saving.nc"))
ds[:sortedC][1, :]
close(ds)
