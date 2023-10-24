# [F]ast [E]stimated [E]ffective [M]igration [S]urface

## Trying a different approach since there were some conflicts in the existing feems install instructions.

Create a file called feems.txt
```
nano feems.txt
```
Copy and paste the following
```
# This file may be used to create an environment using:
# $ conda create --name <env> --file <this file>
# platform: linux-64
@EXPLICIT
https://repo.anaconda.com/pkgs/main/linux-64/_libgcc_mutex-0.1-main.conda
https://repo.anaconda.com/pkgs/main/linux-64/blas-1.0-openblas.conda
https://conda.anaconda.org/conda-forge/linux-64/ca-certificates-2022.9.24-ha878542_0.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/ld_impl_linux-64-2.38-h1181459_1.conda
https://repo.anaconda.com/pkgs/main/linux-64/libgfortran4-7.5.0-ha8ba4b0_17.conda
https://repo.anaconda.com/pkgs/main/linux-64/libstdcxx-ng-11.2.0-h1234567_1.conda
https://repo.anaconda.com/pkgs/main/linux-64/libgfortran-ng-7.5.0-ha8ba4b0_17.conda
https://repo.anaconda.com/pkgs/main/linux-64/libgomp-11.2.0-h1234567_1.conda
https://repo.anaconda.com/pkgs/main/linux-64/_openmp_mutex-5.1-1_gnu.conda
https://repo.anaconda.com/pkgs/main/linux-64/libgcc-ng-11.2.0-h1234567_1.conda
https://repo.anaconda.com/pkgs/main/linux-64/c-ares-1.18.1-h7f8727e_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/expat-2.4.9-h6a678d5_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/geos-3.8.0-he6710b0_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/giflib-5.2.1-h7b6447c_0.conda
https://conda.anaconda.org/conda-forge/linux-64/gmp-6.2.1-h58526e2_0.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/icu-58.2-he6710b0_3.conda
https://repo.anaconda.com/pkgs/main/linux-64/jpeg-9e-h7f8727e_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/lerc-3.0-h295c915_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libdeflate-1.8-h7f8727e_5.conda
https://repo.anaconda.com/pkgs/main/linux-64/libev-4.33-h7f8727e_1.conda
https://repo.anaconda.com/pkgs/main/linux-64/libffi-3.4.2-h6a678d5_6.conda
https://repo.anaconda.com/pkgs/main/linux-64/libopenblas-0.3.18-hf726d26_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libwebp-base-1.2.4-h5eee18b_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libxcb-1.15-h7f8727e_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/lz4-c-1.9.3-h295c915_1.conda
https://conda.anaconda.org/conda-forge/linux-64/metis-5.1.0-h58526e2_1006.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/ncurses-6.3-h5eee18b_3.conda
https://repo.anaconda.com/pkgs/main/linux-64/nspr-4.33-h295c915_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/openssl-1.1.1s-h7f8727e_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/pcre-8.45-h295c915_0.conda
https://conda.anaconda.org/conda-forge/linux-64/tbb-2021.5.0-h924138e_1.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/xz-5.2.8-h5eee18b_0.conda
https://conda.anaconda.org/conda-forge/linux-64/yaml-0.2.5-h7f98852_2.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/zlib-1.2.13-h5eee18b_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/glib-2.69.1-he621ea3_2.conda
https://conda.anaconda.org/conda-forge/linux-64/libblas-3.9.0-13_linux64_openblas.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/libedit-3.1.20210910-h7f8727e_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libevent-2.1.12-h8f2d780_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libllvm10-10.0.1-hbcb73fb_5.conda
https://repo.anaconda.com/pkgs/main/linux-64/libnghttp2-1.46.0-hce63b2e_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libpng-1.6.37-hbc83047_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libssh2-1.10.0-h8f2d780_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libxml2-2.9.14-h74e7548_0.conda
https://conda.anaconda.org/conda-forge/linux-64/mpfr-4.1.0-h9202a9a_1.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/readline-8.2-h5eee18b_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/tk-8.6.12-h1ccaba5_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/zstd-1.5.2-ha4553b6_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/dbus-1.13.18-hb2f20db_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/freetype-2.12.1-h4a9f257_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/gstreamer-1.14.0-h28cd5cc_2.conda
https://repo.anaconda.com/pkgs/main/linux-64/krb5-1.19.2-hac12032_0.conda
https://conda.anaconda.org/conda-forge/linux-64/libcblas-3.9.0-13_linux64_openblas.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/libclang-10.0.1-default_hb85057a_2.conda
https://conda.anaconda.org/conda-forge/linux-64/liblapack-3.9.0-13_linux64_openblas.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/libtiff-4.4.0-hecacb30_2.conda
https://repo.anaconda.com/pkgs/main/linux-64/libxkbcommon-1.0.1-hfa300c1_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libxslt-1.1.35-h4e12654_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/sqlite-3.40.0-h5082296_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/fontconfig-2.14.1-hef1e5e3_0.conda
https://conda.anaconda.org/conda-forge/linux-64/gsl-2.6-he838d99_2.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/gst-plugins-base-1.14.0-h8213a91_2.conda
https://repo.anaconda.com/pkgs/main/linux-64/libcurl-7.86.0-h91b91d3_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/libpq-12.9-h16c4e8d_3.conda
https://repo.anaconda.com/pkgs/main/linux-64/libwebp-1.2.4-h11a3e52_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/nss-3.74-h0370c37_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/python-3.8.15-h7a1cb2a_2.conda
https://conda.anaconda.org/conda-forge/linux-64/suitesparse-5.10.1-h9e50725_1.tar.bz2
https://conda.anaconda.org/conda-forge/noarch/attrs-22.1.0-pyh71513ae_1.tar.bz2
https://conda.anaconda.org/conda-forge/noarch/certifi-2022.9.24-pyhd8ed1ab_0.tar.bz2
https://repo.anaconda.com/pkgs/main/noarch/cycler-0.11.0-pyhd3eb1b0_0.conda
https://repo.anaconda.com/pkgs/main/noarch/decorator-5.1.1-pyhd3eb1b0_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/et_xmlfile-1.1.0-py38h06a4308_0.conda
https://repo.anaconda.com/pkgs/main/noarch/jdcal-1.4.1-pyhd3eb1b0_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/joblib-1.1.1-py38h06a4308_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/kiwisolver-1.4.2-py38h295c915_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/numpy-base-1.22.3-py38hb8be1f0_0.conda
https://conda.anaconda.org/conda-forge/noarch/pkgutil-resolve-name-1.3.10-pyhd8ed1ab_0.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/ply-3.11-py38_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/proj-7.0.1-h59a7b90_1.conda
https://repo.anaconda.com/pkgs/main/linux-64/pyparsing-3.0.9-py38h06a4308_0.conda
https://repo.anaconda.com/pkgs/main/noarch/pyshp-2.1.3-pyhd3eb1b0_0.conda
https://conda.anaconda.org/conda-forge/linux-64/python_abi-3.8-2_cp38.tar.bz2
https://conda.anaconda.org/conda-forge/noarch/pytz-2022.6-pyhd8ed1ab_0.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/qt-main-5.15.2-h327a75a_7.conda
https://repo.anaconda.com/pkgs/main/noarch/six-1.16.0-pyhd3eb1b0_1.conda
https://repo.anaconda.com/pkgs/main/noarch/threadpoolctl-2.2.0-pyh0d69192_0.conda
https://repo.anaconda.com/pkgs/main/noarch/toml-0.10.2-pyhd3eb1b0_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/tornado-6.2-py38h5eee18b_0.conda
https://conda.anaconda.org/conda-forge/noarch/typing_extensions-4.4.0-pyha770c72_0.tar.bz2
https://repo.anaconda.com/pkgs/main/noarch/wheel-0.37.1-pyhd3eb1b0_0.conda
https://repo.anaconda.com/pkgs/main/noarch/xlrd-2.0.1-pyhd3eb1b0_0.conda
https://conda.anaconda.org/conda-forge/noarch/zipp-3.11.0-pyhd8ed1ab_0.conda
https://conda.anaconda.org/conda-forge/noarch/importlib-metadata-5.1.0-pyha770c72_0.conda
https://conda.anaconda.org/conda-forge/noarch/importlib_resources-5.10.1-pyhd8ed1ab_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/numpy-1.22.3-py38h7a5d4dd_0.conda
https://repo.anaconda.com/pkgs/main/noarch/openpyxl-3.0.7-pyhd3eb1b0_0.conda
https://repo.anaconda.com/pkgs/main/noarch/packaging-21.3-pyhd3eb1b0_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/pyproj-2.6.1.post1-py38h61f852b_1.conda
https://conda.anaconda.org/conda-forge/linux-64/pyrsistent-0.18.1-py38h0a891b7_1.tar.bz2
https://repo.anaconda.com/pkgs/main/noarch/python-dateutil-2.8.2-pyhd3eb1b0_0.conda
https://conda.anaconda.org/conda-forge/linux-64/pyyaml-5.4.1-py38h497a2fe_1.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/qt-webengine-5.15.9-hd2b0992_4.conda
https://repo.anaconda.com/pkgs/main/linux-64/setuptools-65.5.0-py38h06a4308_0.conda
https://conda.anaconda.org/conda-forge/noarch/svgwrite-1.4.3-pyhd8ed1ab_0.tar.bz2
https://conda.anaconda.org/conda-forge/noarch/jsonschema-4.17.3-pyhd8ed1ab_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/matplotlib-base-3.2.2-py38hef1b27d_0.conda
https://repo.anaconda.com/pkgs/main/noarch/networkx-2.4-py_1.conda
https://conda.anaconda.org/conda-forge/linux-64/pandas-1.4.2-py38h47df419_1.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/pip-22.2.2-py38h06a4308_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/qtwebkit-5.212-h4eab89a_4.conda
https://repo.anaconda.com/pkgs/main/linux-64/scipy-1.5.0-py38habc2bb6_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/shapely-1.8.4-py38h81ba7c5_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/sip-6.6.2-py38h6a678d5_0.conda
https://repo.anaconda.com/pkgs/main/linux-64/cartopy-0.18.0-py38hc576cba_1.conda
https://conda.anaconda.org/conda-forge/noarch/patsy-0.5.3-pyhd8ed1ab_0.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/pyqt5-sip-12.11.0-py38h6a678d5_1.conda
https://repo.anaconda.com/pkgs/main/linux-64/scikit-learn-0.23.1-py38h7ea95a0_0.conda
https://conda.anaconda.org/conda-forge/linux-64/scikit-sparse-0.4.4-py38h956cf04_1004.tar.bz2
https://conda.anaconda.org/conda-forge/linux-64/tskit-0.4.1-py38h6c62de6_0.tar.bz2
https://conda.anaconda.org/conda-forge/linux-64/msprime-1.0.0-py38hbf02952_0.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/pyqt-5.15.7-py38h6a678d5_1.conda
https://conda.anaconda.org/conda-forge/linux-64/statsmodels-0.12.2-py38h6c62de6_0.tar.bz2
https://repo.anaconda.com/pkgs/main/linux-64/matplotlib-3.2.2-0.conda
```
Create the conda environment using this file and then activate it
```
micromamba create --name feems --file feems.txt
```
I get lots of yellow warnings here but let's see how we do.

Download the feems software from github in your bin
```
git clone https://github.com/NovembreLab/feems
```

cd into the feems directoty and install the software
```
cd feems
pip install .
```

## Example Tutorial
See the following [notebook](https://github.com/NovembreLab/feems/blob/main/docsrc/notebooks/getting-started.ipynb) for an example workflow

## You can run through this tutorial interactively on farm
```
srun -p bmh --time=10:00:00 --nodes=1 --cpus-per-task 4 --mem 50GB --pty /bin/bash
```
Activate the feems environment
```
micromamba activate feems
```
## You'll need to open up the python interpreter
```
python

## should see something like this
Python 3.8.15 (default, Nov 24 2022, 15:19:38) 
[GCC 11.2.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
```
## Import the required packages and feems
```
# base
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink

# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"
```
## Read in data
```
data_path = pkg_resources.resource_filename("feems", "data/")

#I'll use the same data path for the kit foxes
```
## Read the plink formatted genotype data and impute any missing SNPs with the mean at each SNP:
```
## Wolf Data ##
(bim, fam, G) = read_plink("{}/wolvesadmix".format(data_path))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

## Kit Fox Data ##
(bim, fam, G) = read_plink("{}/AllKF".format(data_path))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))
```
For preparing the graph inputs to run feems you have two options:

1. Prepare your own input files 
2. Use the feems function prepare_graph_inputs which intersects a discrete global grid (DGG) with the sample range

The tutorial runs through the latter option. More specifically *"We read the sample coordinates, coordinates of the outer polygon that defines the habitat of the sample and a discrete global grid file which has laid down a triangular grid that is uniformly spaced on earth. We then intersect this global grid with the outer file to define the graph that we use to optimize:"*

Note that you can use this [website](http://www.birdtheme.org/useful/v3tool.html) to help generate the outer coordinates file 
```
## Wolf Data ##
# setup graph
coord = np.loadtxt("{}/wolvesadmix.coord".format(data_path))  # sample coordinates
outer = np.loadtxt("{}/wolvesadmix.outer".format(data_path))  # outer coordinates


##Kit Fox Data ##
# setup graph
coord = np.loadtxt("{}/AllKF.coord".format(data_path))  # sample coordinates
outer = np.loadtxt("{}/AllKF.outer".format(data_path))  # outer coordinates

grid_path = "{}/grid_100.shp".format(data_path)  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=True, 
                                             buffer=0,
                                             outer=outer)
```

## Setup the spatialgraph object
```
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)
```

```
## For Wolf Data ##
projection = ccrs.EquidistantConic(central_longitude=-108.842926, central_latitude=66.037547)

## For KitFox Data ##
import cartopy.crs as ccrs
projection = ccrs.AlbersEqualArea(central_longitude=-120, central_latitude=40,
                                  standard_parallels=(34, 42))
```

```
import matplotlib.pyplot as plt

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=10, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)

# Save the figure to an output file
fig.savefig("output_spatialgraphobject.png")

```

## Fit feems
Next we fit a the feems model where we allow a weight to be estimated for every edge, which is encoded in a large adjacency matrix , while encouraging nearby edges to be smooth. We initialize at the fit from the null model and fix the estimate of the residual variance for the more complex optimization:
```
sp_graph.fit(lamb = 20.0)
```
Now we can visualize the weighted graph:
```
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar()

# Save the figure to an output file
fig.savefig("output_fitfeems.png")
```
Lets now try a different regularization setting that isn't as smooth:
```
sp_graph.fit(lamb = 2.0)
```
```
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar()

# Save the figure to an output file
fig.savefig("output_fitfeems_lam2.png")
```

## Choose a lamdda-value using cross-validation 
```
LoadCVFromDisk = False

# define grid
# Publication grid for this dataset
#lamb_grid = np.geomspace(1e-6, 1e2, 20)[::-1]
# Exploratory grid: 
lamb_grid = np.geomspace(1e-4, 1e1, 10)[::-1]

# run cross-validation
if not LoadCVFromDisk:
    cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)
    pickle.dump(cv_err,open("cv_err.pkl","wb"))
```
Plot the CV error
```
LoadCVFromDisk = TRUE
if LoadCVFromDisk: 
    cv_err = pickle.load(open("cv_err.pkl","rb"))
    
# average over folds
mean_cv_err = np.mean(cv_err, axis=0)

# argmin of cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])

fig, ax = plt.subplots(dpi=300)
ax.plot(np.log10(lamb_grid), mean_cv_err, ".");
ax.set_xlabel("log10(lambda)");
ax.set_ylabel("L2 CV Error");
ax.axvline(np.log10(lamb_cv), color = "orange")
lamb_cv

# Save the figure to an output file
fig.savefig("output_CVerr.png")

```
