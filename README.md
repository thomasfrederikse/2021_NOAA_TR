# Code and data for Section 2 of the Interagency report: Global and Regional Sea Level Rise Scenarios for the United States: Updated Mean Projections and Extreme Water Level Probabilities Along U.S. Coastlines
This repository contains almost all code and data needed to compute the trajectories, projections, and observations for the Interagency report: Global and Regional Sea Level Rise Scenarios for the United States: Updated Mean Projections and Extreme Water Level Probabilities Along U.S. Coastlines.

## Authors
- William V. Sweet, NOAA National Ocean Service 
- Benjamin D. Hamlington, NASA Jet Propulsion Laboratory 
- Robert E. Kopp, Rutgers University 
- Christopher P. Weaver, U.S. Environmental Protection Agency 
- Patrick L. Barnard, U.S. Geological Survey 
- Michael Craghan, U.S. Environmental Protection Agency 
- Gregory Dusek, NOAA National Ocean Service 
- Thomas Frederikse, NASA Jet Propulsion Laboratory 
- Gregory Garner, Rutgers University 
- Ayesha S. Genz, University of Hawai‘i at Mānoa, Cooperative Institute for Marine and Atmospheric Research 
- John P. Krasting, NOAA Geophysical Fluid Dynamics Laboratory 
- Eric Larour, NASA Jet Propulsion Laboratory 
- Doug Marcy, NOAA National Ocean Service 
- John J. Marra, NOAA National Centers for Environmental Information 
- Jayantha Obeysekera, Florida International University 
- Mark Osler, NOAA National Ocean Service 
- Matthew Pendleton, Lynker 
- Daniel Roman, NOAA National Ocean Service 
- Lauren Schmied, FEMA Risk Management Directorate 
- William C. Veatch, U.S. Army Corps of Engineers 
- Kathleen D. White, U.S. Department of Defense 
- Casey Zuzak, FEMA Risk Management Directorate

## Contents
This data and code set contains the following directories:
### Results
The `Results` folder contains the resulting projections, trajectories and observations from the report.

- `TR_global_projections.nc`: GMSL projections, trajectory, and observations 
- `TR_regional_projections.nc`: Regional observations, projections and trajectories 
- `TR_local_projections.nc`: Local observations, projections and trajectories
- `TR_gridded_projections.nc`: Gridded projections

These files are in the NetCDF forrmat. To read the NetCDF files, many free software packages are available, including [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html) and [Panoply](https://www.giss.nasa.gov/tools/panoply/). Free NetCDF packages are available to directly import the data into [Julia](https://github.com/Alexander-Barth/NCDatasets.jl) and [Python](https://unidata.github.io/netcdf4-python/) code.

### Code
The `Code` folder contains all the computer code used to read and analyze the observations and the projections, and to generate the trajectories. 

To run this code, you need [Julia](https://julialang.org/). The code requires the Julia packages `CSV`, `Interpolations`, `JSON`,  `LoopVectorization`,  `MAT`,  `NCDatasets`, `NetCDF`, `Plots`, `XLSX`, `LinearAlgebra`, and `Statistics`. They can be installed by pressing `]` at the Julia REPL and typing:
```
add CSV Interpolations JSON LoopVectorization MAT NCDatasets NetCDF Plots XLSX LinearAlgebra Statistics
```

This program also requires [Hector](http://segal.ubi.pt/hector/). Hector needs to be installed or compiled. In the file `Hector.jl` update the path to the Hector executable on lines 30 and 104. 

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
The `Data` directory contains the input data sets used during the computations. Please appropriately cite the source data if used. It contains the following:

Directories:
- `ClimIdx`: Map with climate indices (NAO, PDO, MEI) used to remove internal variability. All the indices come from NOAA [Physical Sciences Laboratory (PSL)](https://psl.noaa.gov/data/climateindices/) and [NOAA Climate Prediction Centre (CPC)](https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml)
- `NCA5_projections` Contains the NCA5 projections for each scenario (Low, IntLow, Int, IntHigh, and High). For each scenario, the GMSL projections, projections at tide-gauge locations and on a 1-degree grid are provided.  

Files:
- `basin_codes.nc`: Map with basin codes. from Eric Leuliette/NOAA. Data provided by the NOAA Laboratory for Satellite Altimetry.
- `CDS_monthly_1993_2020.nc`: Monthly-mean sea level (1993-2020) from gridded altimetry. Obtained from [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-sea-level-global). This dataset contains modified Copernicus Climate Change Service information [2020]
- `enso_correction.mat`: GMSL correction for ENSO/PDO from Hamlington, B. D., Frederikse, T., Nerem, R. S., Fasullo, J. T., & Adhikari, S. (2020). Investigating the Acceleration of Regional Sea‐level Rise During the Satellite Altimeter Era. Geophysical Research Letters. https://doi.org/10.1029/2019GL086528
- `filelist_psmsl.txt`: List with PSMSL file names and PSMSL IDs. Obtained from the Permanent Service for Mean Sea Level ([PSMSL](http://www.psmsl.org/)), 2021, Retrieved 29 Nov 2021. Simon J. Holgate, Andrew Matthews, Philip L. Woodworth, Lesley J. Rickards, Mark E. Tamisiea, Elizabeth Bradshaw, Peter R. Foden, Kathleen M. Gordon, Svetlana Jevrejeva, and Jeff Pugh (2013) New Data Systems and Products at the Permanent Service for Mean Sea Level. Journal of Coastal Research: Volume 29, Issue 3: pp. 493 – 504. https://doi.org/:10.2112/JCOASTRES-D-12-00175.1.
- `GEBCO_bathymetry_05.nc`: Bathymetry map of the global oceans from the General Bathymetric Chart of the Oceans ([GEBCO](https://www.gebco.net/)). Source: GEBCO Compilation Group (2021) GEBCO 2021 Grid (`doi:10.5285/c6612cbe-50b3-0cff-e053-6c86abc09f8f`) The source data have been re-gridded onto a 0.5 degree grid. 
- `GIA_Caron_stats_05.nc`: Glacial Isostatic Adjustment estimates from Caron, L., Ivins, E. R., Larour, E., Adhikari, S., Nilsson, J., & Blewitt, G. (2018). GIA Model Statistics for GRACE Hydrology, Cryosphere, and Ocean Science. Geophysical Research Letters, 45(5), 2203–2212. <https://doi.org/10.1002/2017GL076644>. The source data have been re-gridded onto a 0.5 degree grid. 
- `global_timeseries_measures.nc`: Time series of estimated 20th-century GMSL and its components, based on Frederikse, T., Landerer, F., Caron, L., Adhikari, S., Parkes, D., Humphrey, V. W., Dangendorf, S., Hogarth, P., Zanna, L., Cheng, L., & Wu, Y.-H. (2020). The causes of sea-level rise since 1900. Nature, 584(7821), 393–397. https://doi.org/10.1038/s41586-020-2591-3
- `GMSL_ensembles.nc`: Ensemble GMSL reconstruction from tide-gauges based on Frederikse, T., Landerer, F., Caron, L., Adhikari, S., Parkes, D., Humphrey, V. W., Dangendorf, S., Hogarth, P., Zanna, L., Cheng, L., & Wu, Y.-H. (2020). The causes of sea-level rise since 1900. Nature, 584(7821), 393–397. https://doi.org/10.1038/s41586-020-2591-3
- `GMSL_TPJAOS_5.0_199209_202106.txt`: Global Mean Sea Level Trend from Integrated Multi-Mission Ocean Altimeters TOPEX/Poseidon, Jason-1, OSTM/Jason-2, and Jason-3 Version 5.1 [Data set]. NASA Physical Oceanography DAAC. https://doi.org/10.5067/GMSLM-TJ151. This altimetry dataset uses the methods as described in Beckley, B. D., Callahan, P. S., Hancock, D. W., Mitchum, G. T., & Ray, R. D. (2017). On the “Cal-Mode” Correction to TOPEX Satellite Altimetry and Its Effect on the Global Mean Sea Level Time Series. Journal of Geophysical Research: Oceans, 122(11), 8371–8384. https://doi.org/10.1002/2017JC013090
- `grd_1992_2020.nc`: Seafloor deformation due to contemporary GRD effects based on Frederikse, T., Landerer, F., Caron, L., Adhikari, S., Parkes, D., Humphrey, V. W., Dangendorf, S., Hogarth, P., Zanna, L., Cheng, L., & Wu, Y.-H. (2020). The causes of sea-level rise since 1900. Nature, 584(7821), 393–397. https://doi.org/10.1038/s41586-020-2591-3
- `region_mask.nc`: Mask with the definition of all regions.
- `US_tg_monthly.xlsx`: Tide gauge observations from the NOAA tide gauge network

### GMT
This directory contains the [GMT](https://www.generic-mapping-tools.org/) scripts to make Figures 1.2, 2.1, 2.2, 2.6, and A.1.2 from the report. To generate the figures, make sure GMT is installed and run the Shell script in each directory. 