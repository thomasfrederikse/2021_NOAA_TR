# ---------------------------------------------------------------
# Read all the NCA5 files and reformat them to a file 
# thats usable throughout the report. All heights in mm
# 
# For each scenario, the following files are made:
#   NCA5_scn_gmsl.nc: Global-mean sea level for each scenario
#   NCA5_scn_grid.nc: Local sea level on a 1° grid
#   NCA5_scn_tg.nc: Local sea level at PSMSL tide-gauge locations
# This script only has to be run once.
# ---------------------------------------------------------------
module ConvertNCA5ToGrid

using NetCDF
using DelimitedFiles
using Statistics
using NCDatasets
using Interpolations
using DelimitedFiles

function RunConvertNCA5ToGrid(settings)
    println("\nProcessing NCA5 scenarios...")
    ConvertGMSL(settings)
    println("  Processing gridded and tide-gauge scenarios...")
    for scn ∈ 1:length(settings["NCA5_scenarios"])
        ConvertIndivScn(scn,settings)
    end
    println("  Processing gridded and tide-gauge scenarios done")
    println("Processing NCA5 scenarios done\n")
    return nothing
end

function ConvertGMSL(settings)
    # -------------------------------------------
    # GMSL is in different file
    # Read the truncated file with GMSL values
    # Save as netcdf
    # Currently, I don't have the GMSL components
    # -------------------------------------------
    println("  Processing GMSL done")
    fn_in = settings["dir_NCA5_raw"]*"NCA5_GMSL.csv"
    data_raw = readdlm(fn_in,',',skipstart=1)
    traw = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120,2130,2140, 2150, 2200, 2250, 2300]
    gmsl = Dict()
    gmsl["years"] = traw
    for (idx,scn) ∈ enumerate(settings["NCA5_scenarios"])
        gmsl[scn] = 10 .* convert.(Int,data_raw[3*(idx-1)+1:3*idx,7:end]')
        @. gmsl[scn] = gmsl[scn][:,[2,1,3]]
    end
    # Save as netcdf
    for (idx,scn) ∈ enumerate(settings["NCA5_scenarios"])
        fn_gmsl_out = settings["dir_NCA5"] * "NCA5_"*scn*"_gmsl.nc"
        fh = Dataset(fn_gmsl_out,"c")
        defDim(fh,"years", length(gmsl["years"]))
        defDim(fh,"percentiles",3)
        defVar(fh,"years",trunc.(Int16,gmsl["years"]),("years",),deflatelevel=6)
        defVar(fh,"percentiles",trunc.(Int16,[17,50,83]),("percentiles",),deflatelevel=6)
        defVar(fh,"total",trunc.(Int16,gmsl[scn]),("years","percentiles",),deflatelevel=6)
        close(fh)
    end
    println("  Processing GMSL done")
end

function ConvertIndivScn(scn,settings)
    # -------------------------------------------
    # Read all processes for individual scenarios
    # Save lat/lon grid and tide-gauge values in
    # individual files
    # -------------------------------------------
    println("   Scenario "*settings["NCA5_scenarios"][scn]*"...")
    scn_grid = Dict()
    scn_tg   = Dict()
    scenario_target = [30,50,100,150,200,250]

    # Read file info
    fn = settings["dir_NCA5_raw"] * "gmsl"*lpad(scenario_target[1],3,"0")*"/AIS_gmsl"*lpad(scenario_target[1],3,"0")*".nc"
    quantiles_raw = ncread(fn,"quantiles")
    idx_low = findfirst(quantiles_raw.>0.1699)
    idx_med = findfirst(quantiles_raw.>0.4999)
    idx_high = findfirst(quantiles_raw.>0.8299)
    scn_grid["quantiles"] = quantiles_raw[[idx_low,idx_med,idx_high]]
    scn_tg["quantiles"] = quantiles_raw[[idx_low,idx_med,idx_high]]
    scn_grid["years"] = ncread(fn,"years")
    scn_tg["years"] = ncread(fn,"years")
    acc_loc_grid = ncread(fn,"locations") .> 4000
    acc_loc_tg = ncread(fn,"locations") .< 4000
    latmat = reshape(ncread(fn,"lat")[acc_loc_grid],(360,181))
    lonmat = reshape(ncread(fn,"lon")[acc_loc_grid],(360,181))
    @. lonmat = mod(lonmat,360)
    scn_grid["lon"] = lonmat[:,1]
    scn_grid["lat"] = reverse(latmat[1,:])

    scn_tg["loc"] = ncread(fn,"locations")[acc_loc_tg]
    scn_tg["lon"] = mod.(ncread(fn,"lon")[acc_loc_tg],360)
    scn_tg["lat"] = ncread(fn,"lat")[acc_loc_tg]

    # Read individual fields
    for process ∈ settings["processes"]
        println("   Reading "*process)
        fn = settings["dir_NCA5_raw"] * "gmsl"*lpad(scenario_target[scn],3,"0")*"/"*process*"_gmsl"*lpad(scenario_target[scn],3,"0")*".nc"
        if (process == "landwaterstorage") || (process == "verticallandmotion")
            scn_grid[process] = reshape(ncread(fn,"sea_level_change",start=[findfirst(acc_loc_grid),2,1],count=[-1,-1,-1])[:,:,[idx_low,idx_med,idx_high]],(360,181,length(scn_grid["years"]),3));
            scn_tg[process] = ncread(fn,"sea_level_change",start=[1,2,1],count=[sum(acc_loc_tg),-1,-1])[:,:,[idx_low,idx_med,idx_high]];
        else
            scn_grid[process] = reshape(ncread(fn,"sea_level_change",start=[findfirst(acc_loc_grid),1,1],count=[-1,-1,-1])[:,:,[idx_low,idx_med,idx_high]],(360,181,length(scn_grid["years"]),3));
            scn_tg[process]  = ncread(fn,"sea_level_change",start=[1,1,1],count=[sum(acc_loc_tg),-1,-1])[:,:,[idx_low,idx_med,idx_high]];
        end
        scn_grid[process] = convert.(Float32,scn_grid[process]);
        scn_tg[process] = convert.(Float32,scn_tg[process]);
        reverse!(scn_grid[process],dims=2);
    end

    # Estimate climate-driven sea-level changes (ice, sterodynamic, land water storage, but not VLM)
    # Assume errors are independent, so they're added in quadrature
    println("   Summing climate-related fields...")

    scn_grid["total_climate"] = zeros(Float32,size(scn_grid["total"]))
    scn_tg["total_climate"] = zeros(Float32,size(scn_tg["total"]))
    for process ∈ settings["processes_clim"]
        @. scn_grid["total_climate"][:,:,:,2] += scn_grid[process][:,:,:,2]
        @. scn_grid["total_climate"][:,:,:,1] += (abs((scn_grid[process][:,:,:,1] -  scn_grid[process][:,:,:,2]))^2)
        @. scn_grid["total_climate"][:,:,:,3] += (abs((scn_grid[process][:,:,:,3] -  scn_grid[process][:,:,:,2]))^2)
        @. scn_tg["total_climate"][:,:,2] += scn_tg[process][:,:,2]
        @. scn_tg["total_climate"][:,:,1] += (abs((scn_tg[process][:,:,1] -  scn_tg[process][:,:,2]))^2)
        @. scn_tg["total_climate"][:,:,3] += (abs((scn_tg[process][:,:,3] -  scn_tg[process][:,:,2]))^2)
    end
    @. scn_grid["total_climate"][:,:,:,1] = scn_grid["total_climate"][:,:,:,2] - sqrt(@. scn_grid["total_climate"][:,:,:,1])
    @. scn_grid["total_climate"][:,:,:,3] = scn_grid["total_climate"][:,:,:,2] + sqrt(@. scn_grid["total_climate"][:,:,:,3])
    @. scn_tg["total_climate"][:,:,1] = scn_tg["total_climate"][:,:,2] - sqrt(@. scn_tg["total_climate"][:,:,1])
    @. scn_tg["total_climate"][:,:,3] = scn_tg["total_climate"][:,:,2] + sqrt(@. scn_tg["total_climate"][:,:,3])

    # Fix NaNs
    @. scn_grid["total_climate"][scn_grid["total_climate"]<-32000] = -32000 
    @. scn_grid["total_climate"][scn_grid["total_climate"]>32000] = -32000 
    @. scn_tg["total_climate"][scn_tg["total_climate"]<-32000] = -32000 
    @. scn_tg["total_climate"][scn_tg["total_climate"]>32000] = -32000 
    for process ∈ settings["processes"]
        @. scn_grid[process][scn_grid[process]<-32000] = -32000 
        @. scn_grid[process][scn_grid[process]>32000] = -32000 
        @. scn_tg[process][scn_tg[process]<-32000] = -32000 
        @. scn_tg[process][scn_tg[process]>32000] = -32000 
    end

    # Save scenario file
    println("   Saving")
    fn_grid_out = settings["dir_NCA5"] * "NCA5_"*settings["NCA5_scenarios"][scn]*"_grid.nc"
    fn_tg_out   = settings["dir_NCA5"] * "NCA5_"*settings["NCA5_scenarios"][scn]*"_tg.nc"

    fh = Dataset(fn_grid_out,"c")
    defDim(fh,"years", length(scn_grid["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"lon", length(scn_grid["lon"]))
    defDim(fh,"lat", length(scn_grid["lat"]))
    defVar(fh,"lon",scn_grid["lon"],("lon",),deflatelevel=6)
    defVar(fh,"lat",scn_grid["lat"],("lat",),deflatelevel=6)
    defVar(fh,"years",trunc.(Int16,scn_grid["years"]),("years",),deflatelevel=6)
    defVar(fh,"percentiles",trunc.(Int16,[17,50,83]),("percentiles",),deflatelevel=6)
    for process ∈ settings["processes"]
        defVar(fh,process,trunc.(Int16,scn_grid[process]),("lon","lat","years","percentiles",),deflatelevel=6)
    end
    defVar(fh,"total_climate",trunc.(Int16,scn_grid["total_climate"]),("lon","lat","years","percentiles",),deflatelevel=6)
    close(fh)

    fh = Dataset(fn_tg_out,"c")
    defDim(fh,"years", length(scn_tg["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"loc", length(scn_tg["loc"]))
    defVar(fh,"loc",scn_tg["loc"],("loc",),deflatelevel=6)
    defVar(fh,"lon",scn_tg["lon"],("loc",),deflatelevel=6)
    defVar(fh,"lat",scn_tg["lat"],("loc",),deflatelevel=6)
    defVar(fh,"years",trunc.(Int16,scn_tg["years"]),("years",),deflatelevel=6)
    defVar(fh,"percentiles",trunc.(Int16,[17,50,83]),("percentiles",),deflatelevel=6)
    for process ∈ settings["processes"]
        defVar(fh,process,trunc.(Int16,scn_tg[process]),("loc","years","percentiles",),deflatelevel=6)
    end
    defVar(fh,"total_climate",trunc.(Int16,scn_tg["total_climate"]),("loc","years","percentiles",),deflatelevel=6)
    close(fh)
    println("   Scenario "*settings["NCA5_scenarios"][scn]*" done")
    return nothing
end

end