# Run all the code for the NOAA TR
dir_code = homedir()*"/Projects/2021_NOAA_TR/Code/"
include(dir_code*"ConvertNCA5ToGrid.jl")
include(dir_code*"Masks.jl")
include(dir_code*"ComputeRegionalObs.jl")
include(dir_code*"GlobalProjections.jl")
include(dir_code*"RegionalProjections.jl")
include(dir_code*"LocalProjections.jl")
include(dir_code*"SaveFigureData.jl")

function main()
    settings = DefSettings()
    ConvertNCA5ToGrid.RunConvertNCA5ToGrid(settings) # Only needed once
    ComputeRegionalObs.RunComputeRegionalObs(settings)
    GlobalProjections.RunGlobalProjections(settings)
    RegionalProjections.RunRegionalProjections(settings)
    LocalProjections.RunLocalProjections(settings)
    SaveFigureData.RunSaveFigureData(settings)
    return nothing
end

function DefSettings()
    settings = Dict()
    settings["dir_project"] = homedir()*"/Projects/2021_NOAA_TR/"

    # File names
    settings["fn_region_mask"] = settings["dir_project"]*"Data/region_mask.nc"
    settings["fn_basin_codes"] = settings["dir_project"]*"Data/basin_codes.nc"
    settings["fn_bathymetry"] = settings["dir_project"]*"Data/GEBCO_bathymetry_05.nc"

    # Observed sea level
    settings["fn_tg_data"] = settings["dir_project"]*"Data/US_tg_monthly.xlsx" 
    settings["fn_GMSL_20c"] = settings["dir_project"]* "Data/GMSL_ensembles.nc" 
    settings["fn_GMSL_GSFC"] = settings["dir_project"]*"Data/GMSL_TPJAOS_5.0_199209_202106.txt"

    # Altimetry and GMSL
    settings["fn_gmsl_20c"] = settings["dir_project"]*"Data/global_timeseries_measures.nc"
    settings["fn_GIA"] = settings["dir_project"]*"Data/GIA_Caron_stats_05.nc"
    settings["fn_altimetry"] = settings["dir_project"]*"Data/CDS_monthly_1993_2019.nc"
    settings["fn_GRD"] = settings["dir_project"]*"Data/grd_1992_2020.nc"

    # Indices
    settings["fn_NAO"] = settings["dir_project"]*"Data/ClimIdx/nao_pc_monthly.txt" 
    settings["fn_PDO"] = settings["dir_project"]*"Data/ClimIdx/pdo.txt" 
    settings["fn_MEI_1"] = settings["dir_project"]*"Data/ClimIdx/MEIext.data" 
    settings["fn_MEI_2"] = settings["dir_project"]*"Data/ClimIdx/meiv2.data" 
    settings["fn_GMSL_ENSO"] = settings["dir_project"]*"Data/enso_correction.mat"
    
    # NCA5 projections
    settings["dir_NCA5_raw"] = homedir() * "/Data/NCA5/"
    settings["dir_NCA5"] =settings["dir_project"]*"Data/NCA5_projections/"


    # settings["fn_GMSL_NCA5"] = settings["dir_project"]*"Data/NCA5_RSL_projections_vlm_gmsl.nc" 
    # settings["fn_grid_NCA5"] = settings["dir_project"]*"Data/NCA5_RSL_projections_vlm_grid.nc" 
    # Files to save the NOAA projections
    settings["fn_proj_glb"] = settings["dir_project"]*"Data/NOAA_TR_global_projections.nc" 
    settings["fn_proj_reg"] = settings["dir_project"]*"Data/NOAA_TR_regional_projections.nc" 
    settings["fn_proj_lcl"] = settings["dir_project"]*"Data/NOAA_TR_local_projections.nc" 
    settings["fn_regional_obs"] = settings["dir_project"]* "Data/NOAA_TR_regional_obs.nc"   

    # Scenarios
    settings["NCA5_scenarios"]  = ["Low","IntLow","Int","IntHigh","High"]
    settings["regions"]         = ["USA","EC","SE","GCE","GCW","SWC","NWC","PAC","CAR","ALN","ALS"]
    settings["processes"]       = ["AIS","GIS","glaciers","landwaterstorage","oceandynamics","total","verticallandmotion"]
    settings["processes_clim"]  = ["AIS","GIS","glaciers","landwaterstorage","oceandynamics"]

    # Trajectory years
    settings["years_trajectory"] = [1970:2200...]
    settings["years_tg"] = [1920:2020...]
    settings["years_baseline"] = [1999:2001...]
    settings["years"] = [1900:2200...]

    # Figure directories
    settings["dir_gmt"] = settings["dir_project"]*"GMT/" 
    settings["dir_fig_1_gmsl_usa"] = settings["dir_gmt"]*"fig_1_gmsl_usa/" 
    settings["dir_fig_2_regional"] = settings["dir_gmt"]*"fig_2_regional/" 
    settings["dir_fig_3_map"] = settings["dir_gmt"]*"fig_3_map/" 
    settings["dir_fig_4_gmsl"] = settings["dir_gmt"]*"fig_4_gmsl/" 
    settings["dir_fig_5_divergence"] = settings["dir_gmt"]*"fig_5_divergence/" 
    return settings
end