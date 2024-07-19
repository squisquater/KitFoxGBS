This assumes you have your conductance surface already generated. Code for this is in the LandGenR section but repeated here.

Putting these files in the following directory:
```
/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/kf_omniscape_20240717
```
## Load and plot the raster files
Activate the conda environment
```
micromamba activate LandGen
```
Activate R
```
R
```
Install & load required packages
```
library(LandGenCourse)
library(terra)
library(gdistance)
library(sf)
library(spdep)
```
Load and plot the raster files
```
## Habitat Suitability Model
kfsuit <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/ESRP-kfsuit-continuous-modified-gdal-5-reprojected-1000x1000-nodata.tif")
```
Remove NA raster cells
```
# Use ifel to set 0 values to NA
kfsuit <- ifel(kfsuit == 0, NA, kfsuit)

png("kfsuit-20240719.png", width = 800, height = 600)
terra::plot(kfsuit)
dev.off()
```
```
## Roads Layer
roads <- rast("/group/ctbrowngrp2/sophiepq/KitFoxGBS/LandGenR/InterstateHwy5_rasterized.tif")

### mask the roads layer to only include overlapping raster cells.
roads_clipped <- mask(roads, kfsuit)

png("roads-20240719.png", width = 800, height = 600)
terra::plot(roads_clipped)
dev.off()

## Reclassify roads
# Define the reclassification matrix
# The first column is the current value, and the second column is the new value
reclass_matrix <- matrix(c(0, 0, 1, 50), ncol = 2, byrow = TRUE)

# Apply the reclassification
roads_reclassified <- classify(roads_clipped, rcl = reclass_matrix)

# Plot it
png("roads-20240719-reclass.png", width = 800, height = 600)
terra::plot(roads_reclassified)
dev.off()
```

Create a cost surface
```
kfsuit.cost <- (100-kfsuit[[1]])*0.5
roads.cost <-  roads_clipped[[1]]*50

cost <- (kfsuit.cost + roads.cost)

png("KFfull.CostSurface-20240719.png", width = 800, height = 600)
terra::plot(cost[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```

Create a conductance surface
```

Create a single landscape conductance raster and plot it
```
conductance1 <- (100-cost[[1]])

png("KFfull.ConductanceSurface-20240719.png", width = 800, height = 600)
terra::plot(conductance1[[1]])
terra::points(points, pch=21, col="black", bg="white", cex=2)
dev.off()
```
Save the raster layer for use in Omniscape
```
writeRaster(conductance1, "KFfull_ConductanceSurface_20240719.tif", overwrite=TRUE)
```



## Run Omniscape

Planning to run this on the HPC 

Create a conda environment called Omniscape
```
micromamba create --name Omniscape
```
Activate the environment
```
micromamba activate Omniscape
```
Install the program Julia
```
micromamba install conda-forge::julia
```
Activate julia
```
julia
```
You should see a command prompt returned like this
```
              _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.10.2 (2024-03-01)
 _/ |\__'_|_|_|\__'_|  |  https://github.com/conda-forge/julia-feedstock
|__/                   |

julia>
```
Install the omniscape package using Julia (you should only need to do this the first time.
```
using Pkg; Pkg.add(["Omniscape", "Rasters", "Plots", "DataFrames"])
```
Load the packages
```
using Omniscape, Rasters, Plots
```

I am getting some warnings while running this that I may need to address later. tbd.
```
Warning: CHOLMOD version incompatibility
│ 
│ Julia was compiled with CHOLMOD version 4.0.4. It is
│ currently linked with version 3.0.14.
│ This might cause Julia to terminate when working with
│ sparse matrix factorizations, e.g. solving systems of
│ equations with \.
│ 
│ It is recommended that you use Julia with the same major
│ version of CHOLMOD as the one used during the build, or
│ download the generic binaries from www.julialang.org,
│ which ship with the correct versions of all dependencies.
└ @ SparseArrays.CHOLMOD ~/micromamba/envs/Omniscape/share/julia/stdlib/v1.10/SparseArrays/src/solvers/cholmod.jl:206
  2 dependencies successfully precompiled in 11 seconds. 174 already precompiled.
  2 dependencies had output during precompilation:
```
### NEW 20240717
Load in the .tif file and plot it
```
plot(Raster("/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/kf_omniscape_20240719/KFfull_ConductanceSurface_20240719.tif"), title = "KitFoxFullModel", xlabel = "Easting", ylabel = "Northing", size = (700,640))

savefig("KitFoxFullModel.png")
```
Load the array using Omniscape's internal read_raster() function 
```
land_cover_KF, wkt, transform = Omniscape.read_raster("/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/kf_omniscape_20240719/KFfull_ConductanceSurface_20240719.tif", Float64)
```
Create a resistance reclassification table. Not sure whether this is necessary since I already have my data classified continuously but I want to convert the 0 entries to missing which will prohibit current from moving through these areas.

First I'll create a dataframe of unique values in the raster layer
```
using Plots
using ArchGDAL
using DataFrames

# Load the raster file
raster_path = "/group/ctbrowngrp2/sophiepq/KitFoxGBS/Omniscape/kf_omniscape_20240719/KFfull_ConductanceSurface_20240719.tif"
dataset = ArchGDAL.read(raster_path)
raster_band = ArchGDAL.getband(dataset, 1)

# Convert the raster band to an array
raster_array = ArchGDAL.read(raster_band)

# Find unique values and convert to decimal
unique_values = unique(raster_array)
unique_values_decimal = Int.(unique_values)
println("Unique values in the raster layer (decimal):")
println(unique_values_decimal)

# Create a DataFrame for unique values
unique_values_df = DataFrame(OriginalValue = unique_values_decimal)

# Define the reclassification for the value 0
reclass_mapping = DataFrame(OriginalValue = [0], ReclassifiedValue = ["missing"])

# Add unchanged mappings for other values
unchanged_values = setdiff(unique_values_decimal, [0])
unchanged_mappings = DataFrame(OriginalValue = unchanged_values, ReclassifiedValue = unchanged_values)

# Combine the mappings
reclass_table = vcat(reclass_mapping, unchanged_mappings)

# Ensure ReclassifiedValue column has consistent types (Union{Missing, Float64})
reclass_table.ReclassifiedValue = Union{Missing, Float64}[x == "missing" ? missing : float(x) for x in reclass_table.ReclassifiedValue]

# Print the final reclassification table
println(reclass_table)
```

Note: I used the following code in the QGIC python console to print out all the unique raster values
```
from osgeo import gdal
import numpy as np

dataset = gdal.Open('path/to/your/raster/file.tif')
band = dataset.GetRasterBand(1)
array = band.ReadAsArray()
unique_values = np.unique(array)
print(unique_values)

#[  0.   5.  10.  15.  20.  30.  40.  50.  60.  80.  85.  90.  95.]
```

```
reclass_table = [
    0.0 missing;
    1.0 1.0;
    2.0 2.0;
    3.0 3.0;
    4.0 4.0;
    5.0 5.0;
    6.0 6.0;
    7.0 7.0;
    8.0 8.0;
    9.0 9.0;
    10.0 10.0;
    11.0 11.0;
    12.0 12.0;
    13.0 13.0;
    14.0 14.0;
    15.0 15.0;
    16.0 16.0;
    17.0 17.0;
    18.0 18.0;
    19.0 19.0;
    20.0 20.0;
    21.0 21.0;
    22.0 22.0;
    23.0 23.0;
    24.0 24.0;
    25.0 25.0;
    26.0 26.0;
    27.0 27.0;
    28.0 28.0;
    29.0 29.0;
    30.0 30.0;
    31.0 31.0;
    32.0 32.0;
    33.0 33.0;
    34.0 34.0;
    35.0 35.0;
    36.0 36.0;
    37.0 37.0;
    38.0 38.0;
    39.0 39.0;
    40.0 40.0;
    41.0 41.0;
    42.0 42.0;
    43.0 43.0;
    44.0 44.0;
    45.0 45.0;
    46.0 46.0;
    47.0 47.0;
    48.0 48.0;
    49.0 49.0;
    50.0 50.0;
    51.0 51.0;
    52.0 52.0;
    53.0 53.0;
    54.0 54.0;
    55.0 55.0;
    56.0 56.0;
    57.0 57.0;
    58.0 58.0;
    59.0 59.0;
    60.0 60.0;
    61.0 61.0;
    62.0 62.0;
    63.0 63.0;
    64.0 64.0;
    65.0 65.0;
    66.0 66.0;
    67.0 67.0;
    68.0 68.0;
    69.0 69.0;
    70.0 70.0;
    71.0 71.0;
    72.0 72.0;
    73.0 73.0;
    74.0 74.0;
    75.0 75.0;
    76.0 76.0;
    77.0 77.0;
    78.0 78.0;
    79.0 79.0;
    80.0 80.0;
    81.0 81.0;
    82.0 82.0;
    83.0 83.0;
    84.0 84.0;
    85.0 85.0;
    86.0 86.0;
    87.0 87.0;
    88.0 88.0;
    89.0 89.0;
    90.0 90.0;
    91.0 91.0;
    92.0 92.0;
    93.0 93.0;
    94.0 94.0;
    95.0 95.0;
    96.0 96.0;
    97.0 97.0;
    98.0 98.0;
    99.0 99.0;
    100.0 100.0
]
```

Next, we define the configuration options for this model run.
Note that for animals, the radius should correspond to a typical dispersal distance or the distance an animal is expected to move within its habitat to access resources.
The radius is in grid cells. Our raster layer is 
```
config = Dict{String, String}(
    "radius" => "260",
    "block_size" => "195",
    "project_name" => "kf_omniscape_output",
    "source_from_resistance" => "true",
    "r_cutoff" => "50", # Only high quality pixels should be sources
    "reclassify_resistance" => "true",
    "calc_normalized_current" => "true",
    "calc_flow_potential" => "true",
    "write_reclassified_resistance"=> "true"
)

```
Finally, compute connectivity using run_omniscape()
```
currmap, flow_pot, norm_current = run_omniscape(config,
                                                land_cover_KF,
                                                reclass_table = reclass_table,
                                                wkt = wkt,
                                                geotransform = transform,
                                                write_outputs = true)
```
Now plot the current
```
current = Raster("kf_omniscape_output_2/cum_currmap.tif")
plot(current,
     title = "Cumulative Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (600, 550))

savefig("KFcurrent.png")
```

Next map flow potential. This shows what connectivity looks like under "null" conditions (resistance equals 1 for the whole landscape).
```
fp = Raster("kf_omniscape_output/flow_potential.tif")
plot(fp,
     title = "Flow Potential", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))

savefig("KFflow_potential.png")
```

Finally, map normalized current flow, which is calculated as cumulative current divided by flow potential.
```
normalized_current = Raster("kf_omniscape_output/normalized_cum_currmap.tif")

plot(normalized_current,
     title = "Normalized Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))

savefig("KFnormalized_current.png")


### I'll work through an example to start
Next, download the landcover data we'll use in this example, and plot it:
```
url_base = "https://raw.githubusercontent.com/Circuitscape/datasets/main/"
# Download the NLCD tile used to create the resistance surface and load it
download(string(url_base, "data/nlcd_2016_frederick_md.tif"),
         "nlcd_2016_frederick_md.tif")

# Plot the landcover data
values = [11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95]
palette = ["#476BA0", "#DDC9C9", "#D89382", "#ED0000", "#AA0000",
           "#b2b2b2", "#68AA63", "#1C6330", "#B5C98E", "#CCBA7C",
           "#E2E2C1", "#DBD83D", "#AA7028", "#BAD8EA", "#70A3BA"]

plot(Raster("nlcd_2016_frederick_md.tif"),
     title = "Land Cover Type", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(palette, (values .- 12) ./ 84, categorical = true),
     size = (700, 640))

# You'll get the following error because we're dn't have a GUI to plot this
# GKS: cannot open display - headless operation mode active.
# Just go ahead and save it using the following:

savefig("land_cover_type.png")
```
Now, load the array using Omniscape's internal read_raster() function or a function from a GIS Julia package of your choice. read_raster() returns a tuple with the data array, a wkt string containing geographic projection info, and an array containing geotransform values. We'll use the wkt and geotransform later.
```
land_cover, wkt, transform = Omniscape.read_raster("nlcd_2016_frederick_md.tif", Float64)
```
The next step is to create a resistance reclassification table that defines a resistance value for each land cover value. Land cover values go in the left column, and resistance values go in the right column. In this case, we are modeling forest connectivity, so forest classes receive the lowest resistance score of one. Other "natural" land cover types are assigned moderate values, and human-developed land cover types are assigned higher values. Medium- to high-intensity development are given a value of missing, which denotes infinite resistance (absolute barriers to movement).

```
reclass_table = [
    0 missing; # Water or Forested 
    5 950; 
    10 900; 
    15 850; 
    20 800; 
    30 700; 
    40 600;
    50 500;
    60 400;
    80 200;
    85 150;
    90 100;
    95 50;
]
```
Next, we define the configuration options for this model run.
```
config = Dict{String, String}(
    "radius" => "100",
    "block_size" => "21",
    "project_name" => "KF_omniscape_output",
    "source_from_resistance" => "true",
    "r_cutoff" => "50", # Only HighQuality Habitat should be a source
    "reclassify_resistance" => "true",
    "calc_normalized_current" => "true",
    "calc_flow_potential" => "true"
)
```
Finally, compute connectivity using run_omniscape()
```
currmap, flow_pot, norm_current = run_omniscape(config,
                                                land_cover,
                                                reclass_table = reclass_table,
                                                wkt = wkt,
                                                geotransform = transform,
                                                write_outputs = true)
```
Now plot the current
```
current = Raster("md_nlcd_omniscape_output/cum_currmap.tif")
plot(current,
     title = "Cumulative Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (600, 550))

savefig("current.png")
```
Next map flow potential. This shows what connectivity looks like under "null" conditions (resistance equals 1 for the whole landscape).
```
fp = Raster("md_nlcd_omniscape_output/flow_potential.tif")
plot(fp,
     title = "Flow Potential", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))

savefig("flow_potential.png")
```
*Note: Flow potential shows what connectivity would look like in the absence of barriers to movement. The blocking that you can see is an artifact of setting a large block_size to make the example run faster. Set a smaller block_size to reduce/remove this issue.*

Finally, map normalized current flow, which is calculated as cumulative current divided by flow potential.
```
normalized_current = Raster("md_nlcd_omniscape_output/normalized_cum_currmap.tif")

plot(normalized_current,
     title = "Normalized Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))

savefig("normalized_current.png")
```

## Now let's work through some kit fox data!

Load in the .tif file and plot it
```
plot(Raster("KitFox-ESARPmodel-Raster30x30-NEWmod01.tif"), title = "KitFoxSuitModel", xlabel = "Easting", ylabel = "Northing", size = (700,640))

savefig("KitFoxSuitModel.png")
```
Load the array using Omniscape's internal read_raster() function 
```
land_cover_KF, wkt, transform = Omniscape.read_raster("KitFox-ESARPmodel-Raster30x30-NEWmod01.tif", Float64)
```
Create a resistance reclassification table. Not sure whether this is necessary since I already have my data classified continuously but I want to convert the 0 entries to missing which will prohibit current from movign through these areas. This includes only water and forested areas as per the habitat suitability model from Cypher et al. 2013.

Note: I used the following code in the QGIC python console to print out all the unique raster values
```
from osgeo import gdal
import numpy as np

dataset = gdal.Open('path/to/your/raster/file.tif')
band = dataset.GetRasterBand(1)
array = band.ReadAsArray()
unique_values = np.unique(array)
print(unique_values)

#[  0.   5.  10.  15.  20.  30.  40.  50.  60.  80.  85.  90.  95.]
```

```
reclass_table = [
    1. missing; # Water or Forested
    5 950; 
    10 900; 
    15 850; 
    20 800; 
    30 700; 
    40 600;
    50 500;
    60 400;
    80 200;
    85 150;
    90 100;
    95 50;
]
```
Next, we define the configuration options for this model run.
Note that for animals, the radius should correspond to a typical dispersal distance or the distance an animal is expected to move within its habitat to access resources.
The radius is in grid cells. Our raster layer is 
```
config = Dict{String, String}(
    "radius" => "260",
    "block_size" => "195",
    "project_name" => "kf_omniscape_output",
    "source_from_resistance" => "true",
    "r_cutoff" => "50", # Only high quality pixels should be sources
    "reclassify_resistance" => "true",
    "calc_normalized_current" => "true",
    "calc_flow_potential" => "true",
    "write_reclassified_resistance"=> "true"
)

```
Finally, compute connectivity using run_omniscape()
```
currmap, flow_pot, norm_current = run_omniscape(config,
                                                land_cover_KF,
                                                reclass_table = reclass_table,
                                                wkt = wkt,
                                                geotransform = transform,
                                                write_outputs = true)
```
Now plot the current
```
current = Raster("kf_omniscape_output_2/cum_currmap.tif")
plot(current,
     title = "Cumulative Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (600, 550))

savefig("KFcurrent.png")
```

Next map flow potential. This shows what connectivity looks like under "null" conditions (resistance equals 1 for the whole landscape).
```
fp = Raster("kf_omniscape_output/flow_potential.tif")
plot(fp,
     title = "Flow Potential", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))

savefig("KFflow_potential.png")
```

Finally, map normalized current flow, which is calculated as cumulative current divided by flow potential.
```
normalized_current = Raster("kf_omniscape_output/normalized_cum_currmap.tif")

plot(normalized_current,
     title = "Normalized Current Flow", xlabel = "Easting", ylabel = "Northing",
     seriescolor = cgrad(:inferno, [0, 0.005, 0.03, 0.06, 0.09, 0.14]),
     size = (700, 640))

savefig("KFnormalized_current.png")
```
