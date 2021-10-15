# NOAA Technical report on sea level: code and data [DO NOT DISTRIBUTE]
This repository contains almost all code and data needed to compute the NOAA technical report sea-level projections. It's currently still work-in-progress, but it should be completed in December. 

## Contents
### Code
The code can be found in the `Code` folder. The folder contains the following files:

- `run_NOAA_tr.jl`: This is the main routine that (eventually) calls all the functions to compute the projections.
- `ComputeRegionalObs.jl`: Reads and processes the regional observations
- `GlobalProjections.jl`: Reads and processes the GMSL observations and projections, and computes the trajectory
- `RegionalProjections.jl`: Reads and processes the regional projections and computes the trajectories
- `LocalProjections.jl`: Reads and processes the local projections at the tide-gauge locations and computes the trajectories
- `Simple_plots.jl`: Makes some simple plots and saves the time series for plotting by GMT. Currently in a state of UTTER chaos
- `convert_NCA5_file.jl`: Convert the NCA5 projections file from CSV to NetCDF
- `Hector.jl`: Wrapper for [Hector](http://segal.ubi.pt/hector/), used to compute trends and uncertainties.
- `Masks.jl`: Defines the region masks for each region.

### Data
The Data directory contains both the input data sets and the output. 

- `CDS_monthly_1993_2019.nc`: Gridded altimetry. Obtained from Copernicus CDS
- `ClimIdx`: Map with climate indices (NAO, PDO, MEI) used to remove internal variability
- `GIA_Caron_stats_05.nc`: GIA (Caron et al. 2018)
- `GMSL_TPJAOS_5.0_199209_202106.txt`: Altimetry GMSL (GSFC, Beckley et al.)
- `GMSL_ensembles.nc`: GMSL reconstruction from tide-gauges (Frederikse et al. 2020)
- `NCA5_RSL_projections_vlm.csv`: NCA5 sea-level projections
- `NCA5_RSL_projections_vlm_gmsl.nc`: NCA5 sea-level projections (GMSL)
- `NCA5_RSL_projections_vlm_grid.nc`: NCA5 sea-level projections (Gridded)
- `NCA5_RSL_projections_vlm_tg.nc`: NCA5 sea-level projections (TG locations)
- `NOAA_TR_global_projections.nc`: GMSL projections, trajectory, and observations for the NOAA TR
- `NOAA_TR_regional_projections.nc`: Regional projections and trajectories for the NOAA TR
- `NOAA_TR_regional_obs.nc`: Regional observations for the NOAA TR
- `US_tg_monthly.xls`: NOAA TG Data
- `US_tg_monthly.xlsx`: NOAA TG Data
- `enso_correction.mat`: GMSL correction for ENSO/PDO
- `grd_1992_2020.nc`: Seafloor deformation due to contemporary GRD effects
- `region_mask.nc`: Mask with all regions

### GMT
This directory contains the [GMT](https://www.generic-mapping-tools.org/) scripts to make the final plots

### Dependencies
A lot. All scripts are in Julia. TBD

# Todo
`[ ]` Add dependencies

`[ ]` Reformat convert_NCA5_file

`[ ]` Clean up GMT routines

`[ ]` Reformat into proper package

`[ ]` Altimetry trend maps