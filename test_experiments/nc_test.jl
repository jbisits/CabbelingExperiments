using NCDatasets
# This creates a new NetCDF file /tmp/test.nc.
# The mode "c" stands for creating a new file (clobber)
ds = NCDataset(joinpath(@__DIR__, "outputs/test.nc"),"c")

# Define the dimension "lon" and "lat" with the size 100 and 110 resp.
defDim(ds,"lon",100)
defDim(ds,"lat",110)

# Define a global attribute
ds.attrib["title"] = "this is a test file"

# Define the variables temperature
v = defVar(ds,"temperature",Float32,("lon","lat"))

# Generate some example data
data = [Float32(i+j) for i = 1:100, j = 1:110]

# write a single column
v[:,1] = data[:,1]

# write a the complete data set
v[:,:] = data

# write attributes
v.attrib["units"] = "degree Celsius"
v.attrib["comments"] = "this is a string attribute with Unicode Ω ∈ ∑ ∫ f(x) dx"

close(ds)
