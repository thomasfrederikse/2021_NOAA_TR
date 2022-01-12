# Interagency Technical report on sea level: code and data for Section 2
This repository contains almost all code and data needed to compute the Interagency technical report.

## Contents
### Code
The `Code` folder contains all the computer code used to read and analyze the observations and the projections, and to generate the trajectories. 

To run this code, you need [Julia](https://julialang.org/). The code requires the Julia packages CSV, Interpolations, JSON,  LoopVectorization,  MAT,  NCDatasets, NetCDF, Plots, XLSX, LinearAlgebra, and Statistics.
They can be installed by pressing ] at the Julia REPL and typing
`add CSV Interpolations JSON  LoopVectorization  MAT  NCDatasets NetCDF Plots XLSX LinearAlgebra Statistics`. This program also requires [Hector](http://segal.ubi.pt/hector/). Hector needs to be installed or compiled. In the file `Hector.jl` update the path to the Hector executable on lines 30 and 104. 

Run `Run_TR.jl` in the REPL or run `julia Run_TR.jl` from the command line to run the projections. The projections are then written to the `.\Data` directory. 

The folder contains the following files:
- `Run_TR.jl`: This is the main routine that (eventually) calls all the functions to compute the projections.
- `ConvertNCA5ToGrid.jl`: Converts the original NCA5 projections to a set of netCDF files that's used throughout this code
- `ProcessObservations.jl`: Reads and processes the tide-gauge and altimetry observations
- `GlobalProjections.jl`: Reads and processes the GMSL observations and projections, and computes the trajectory
- `RegionalProjections.jl`: Reads and processes the regional projections and computes the trajectories
- `LocalProjections.jl`: Reads and processes the local projections at the tide-gauge locations and computes the trajectories
- `GriddedProjections.jl`: Reads the gridded NCA5 projections and add a GMSL baseline correction for the 2005 vs 2000 baseline
- `SaveFigureData.jl`: Reads the results and writes text files for GMT
- `Hector.jl`: Wrapper for [Hector](http://segal.ubi.pt/hector/), used to compute trends and uncertainties.
- `Masks.jl`: Defines the region masks for each region.

### Data
The Data directory contains the input data sets. The pre-computed output is also there. 

#### Input data and references:
The following input data has been included for the analysis. Please appropriately cite the source data if used.  
- `CDS_monthly_1993_2020.nc`: Monthly-mean sea level (1993-2020) from gridded altimetry. Obtained from [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-sea-level-global). This dataset contains modified Copernicus Climate Change Service information [2020]
- `GIA_Caron_stats_05.nc`: Glacial Isostatic Adjustment estimates from Caron, L., Ivins, E. R., Larour, E., Adhikari, S., Nilsson, J., & Blewitt, G. (2018). GIA Model Statistics for GRACE Hydrology, Cryosphere, and Ocean Science. Geophysical Research Letters, 45(5), 2203–2212. https://doi.org/10.1002/2017GL076644. The source data have been re-gridded onto a 0.5 degree grid. 

- `GMSL_TPJAOS_5.0_199209_202106.txt`: Global Mean Sea Level Trend from Integrated Multi-Mission Ocean Altimeters TOPEX/Poseidon, Jason-1, OSTM/Jason-2, and Jason-3 Version 5.1 [Data set]. NASA Physical Oceanography DAAC. https://doi.org/10.5067/GMSLM-TJ151. This altimetry dataset uses the methods as described in Beckley, B. D., Callahan, P. S., Hancock, D. W., Mitchum, G. T., & Ray, R. D. (2017). On the “Cal-Mode” Correction to TOPEX Satellite Altimetry and Its Effect on the Global Mean Sea Level Time Series. Journal of Geophysical Research: Oceans, 122(11), 8371–8384. https://doi.org/10.1002/2017JC013090

- `GMSL_ensembles.nc`: Ensemble GMSL reconstruction from tide-gauges based on Frederikse, T., Landerer, F., Caron, L., Adhikari, S., Parkes, D., Humphrey, V. W., Dangendorf, S., Hogarth, P., Zanna, L., Cheng, L., & Wu, Y.-H. (2020). The causes of sea-level rise since 1900. Nature, 584(7821), 393–397. https://doi.org/10.1038/s41586-020-2591-3

- `US_tg_monthly.xlsx`: Tide gauge observations from the NOAA tide gauge network
- `enso_correction.mat`: GMSL correction for ENSO/PDO from Hamlington, B. D., Frederikse, T., Nerem, R. S., Fasullo, J. T., & Adhikari, S. (2020). Investigating the Acceleration of Regional Sea‐level Rise During the Satellite Altimeter Era. Geophysical Research Letters. https://doi.org/10.1029/2019GL086528

- `grd_1992_2020.nc`: Seafloor deformation due to contemporary GRD effects based on Frederikse, T., Landerer, F., Caron, L., Adhikari, S., Parkes, D., Humphrey, V. W., Dangendorf, S., Hogarth, P., Zanna, L., Cheng, L., & Wu, Y.-H. (2020). The causes of sea-level rise since 1900. Nature, 584(7821), 393–397. https://doi.org/10.1038/s41586-020-2591-3

- `region_mask.nc`: Mask with all regions

- `ClimIdx`: Map with climate indices (NAO, PDO, MEI) used to remove internal variability. All the indices come from NOAA [Physical Sciences Laboratory (PSL)](https://psl.noaa.gov/data/climateindices/) and [NOAA Climate Prediction Centre (CPC)](https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml)

- `NCA5_projections` Contains the NCA5 projections for each scenario (Low, IntLow, Int, IntHigh, and High). For each scenario, the GMSL projections, projections at tide-gauge locations and on a 1-degree grid are provided.  

#### Output data
- `TR_global_projections.nc`: GMSL projections, trajectory, and observations for the report
- `TR_regional_projections.nc`: Regional observations, projections and trajectories for the report
- `TR_local_projections.nc`: Local observations, projections and trajectories for the report
- `TR_gridded_projections.nc`: Gridded projections for the report

To read the NetCDF files, many free software is available, including [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html) and [Panoply](https://www.giss.nasa.gov/tools/panoply/). NetCDF packages are available for Julia and Python as well to directly import the data. 

### GMT
This directory contains the [GMT](https://www.generic-mapping-tools.org/) scripts to make Figures 1.2, 2.1, 2.2, 2.6, and A.1.2 from the report. To generate the figures, make sure GMT is installed and run the Shell script in each directory. 