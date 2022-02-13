# -------------------------------------------------------------------------------------------------
# Interagency technical report "Global and Regional Sea Level Rise Scenarios for the United States: 
# Updated Mean Projections and Extreme Water Level Probabilities Along U.S. Coastlines" 
# 
# This script calls all the routines to compute and analyze the projections, trajectories and 
# observations. All results are stored in netCDF files. This script also stores the results for
# some of the figures in the report as text files, which can be read by the GMT scripts used
# to make these figures.
# -------------------------------------------------------------------------------------------------

dir_code = homedir()*"/Projects/2021_NOAA_TR/Code/"
include(dir_code*"ConvertNCA5ToGrid.jl")
include(dir_code*"Masks.jl")
include(dir_code*"GlobalProjections.jl")
include(dir_code*"RegionalProjections.jl")
include(dir_code*"LocalProjections.jl")
include(dir_code*"GriddedProjections.jl")
include(dir_code*"SaveFigureData.jl")

function main()
    # Define some projection settings
    settings = DefSettings()

    # Routines to do some preprocessing. Not needed to run: output has been provided
    # Currently commented out:
    Masks.CreateMask(settings)                       # Create a region mask for each of the regions
    #ConvertNCA5ToGrid.RunConvertNCA5ToGrid(settings) # Covert the NCA5 projections to single GMSL, tide-gauge and gridded projections
    
    # Create the projections
    GlobalProjections.RunGlobalProjections(settings)
    RegionalProjections.RunRegionalProjections(settings)
    LocalProjections.RunLocalProjections(settings)
    GriddedProjections.RunGriddedProjections(settings)

    # Save the data for the report figures
    SaveFigureData.RunSaveFigureData(settings)
    return nothing
end

function DefSettings()
    settings = Dict()
    settings["dir_project"] = homedir()*"/Projects/2021_NOAA_TR/"

    # File names
    settings["fn_region_mask"] = settings["dir_project"] * "Data/region_mask.nc"
    settings["fn_basin_codes"] = settings["dir_project"] * "Data/basin_codes.nc"
    settings["fn_bathymetry"]  = settings["dir_project"] * "Data/GEBCO_bathymetry_05.nc"

    # Observed sea level
    settings["fn_tg_data"]             = settings["dir_project"] * "Data/US_tg_monthly.xlsx" 
    settings["fn_gmsl_20c_ensembles"]  = settings["dir_project"] * "Data/GMSL_ensembles.nc" 
    settings["fn_gmsl_20c_mean"]       = settings["dir_project"] * "Data/global_timeseries_measures.nc"
    settings["fn_gmsl_altimetry_GSFC"] = settings["dir_project"] * "Data/GMSL_TPJAOS_5.0_199209_202106.txt"
    settings["fn_GIA"]                 = settings["dir_project"] * "Data/GIA_Caron_stats_05.nc"
    settings["fn_altimetry"]           = settings["dir_project"] * "Data/CDS_monthly_1993_2020.nc"
    settings["fn_GRD"]                 = settings["dir_project"] * "Data/grd_1992_2020.nc"

    # Climate indices
    settings["fn_gmsl_ENSO"] = settings["dir_project"]*"Data/ClimIdx/enso_correction.mat"
    settings["fn_NAO"] = settings["dir_project"]*"Data/ClimIdx/nao_pc_monthly.txt" 
    settings["fn_PDO"] = settings["dir_project"]*"Data/ClimIdx/pdo.txt" 
    settings["fn_MEI_1"] = settings["dir_project"]*"Data/ClimIdx/MEIext.data" 
    settings["fn_MEI_2"] = settings["dir_project"]*"Data/ClimIdx/meiv2.data" 
    
    # NCA5 projections
    settings["dir_NCA5_raw"] = homedir() * "/Data/NCA5/"
    settings["dir_NCA5"] =settings["dir_project"]*"Data/NCA5_projections/"

    # Files to save the Technical Report projections
    settings["fn_proj_glb"] = settings["dir_project"]*"Results/TR_global_projections.nc" 
    settings["fn_proj_reg"] = settings["dir_project"]*"Results/TR_regional_projections.nc" 
    settings["fn_proj_lcl"] = settings["dir_project"]*"Results/TR_local_projections.nc" 
    settings["fn_proj_gri"] = settings["dir_project"]*"Results/TR_gridded_projections.nc" 

    # Scenarios
    settings["NCA5_scenarios"]  = ["Low","IntLow","Int","IntHigh","High"]
    settings["regions"]         = ["USA","EC", "SE", "GCE",  "GCW",  "SWC","NWC","PAC","CAR","ALN","ALS"] # V1.0 names
    settings["regions"]         = ["USA","NE", "SE", "EGOM", "WGOM", "SW", "NW", "PAC","CAR","NAL","SAL"] # V1.1 names
    settings["processes"]       = ["AIS","GIS","glaciers","landwaterstorage","oceandynamics","total","verticallandmotion"]
    settings["processes_clim"]  = ["AIS","GIS","glaciers","landwaterstorage","oceandynamics"]

    # Trajectory years
    settings["years"] = [1900:2150...] 
    settings["years_trajectory_global"] = [1970:2100...] 
    settings["years_trajectory"] = [1970:2050...] 
    settings["years_tg"] = [1920:2020...]
    settings["years_baseline"] = [2000]

    # Figure directories
    settings["dir_gmt"] = settings["dir_project"]*"GMT/" 
    settings["dir_fig_1_gmsl_usa"] = settings["dir_gmt"]*"fig_1_gmsl_usa/" 
    settings["dir_fig_2_regional"] = settings["dir_gmt"]*"fig_2_regional/" 
    settings["dir_fig_3_map"] = settings["dir_gmt"]*"fig_3_map/" 
    settings["dir_fig_4_gmsl"] = settings["dir_gmt"]*"fig_4_gmsl/" 
    settings["dir_fig_5_divergence"] = settings["dir_gmt"]*"fig_5_divergence/" 
    return settings
end

main()