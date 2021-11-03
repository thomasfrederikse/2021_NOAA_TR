# ---------------------------------------------
# Read all the NCA5 files and extract 
# the grids and tide-gauge locations
# ---------------------------------------------
using NetCDF
using DelimitedFiles
using Statistics
using NCDatasets
using Interpolations

function main()
    settings = def_settings()
    for scn ∈ 1:length(settings["scenario_name"])
        reformat_scn(scn,settings)
    end
    return nothing
end

function def_settings()
    settings = Dict()
    settings["dir_NCA5"] = homedir()*"/Data/NCA5/"
    settings["dir_proj"] = homedir()*"/Projects/2021_NOAA_TR/Data/NCA5_grids/"
    settings["scenario_name"] = ["Low","IntLow","Int","IntHigh","High","Extreme"]
    settings["scenario_target"] = [30,50,100,150,200,250]
    settings["processes"] = ["AIS","GIS","glaciers","landwaterstorage","oceandynamics","total","verticallandmotion"]
    settings["processes_clim"] = ["AIS","GIS","glaciers","landwaterstorage","oceandynamics"]
    return settings
end

function reformat_scn(scn,settings)
    println("Scenario "*settings["scenario_name"][scn])
    scn_grid = Dict()
    scn_tg   = Dict()

    # Prepare
    fn = settings["dir_NCA5"] * "gmsl"*lpad(settings["scenario_target"][1],3,"0")*"/AIS_gmsl"*lpad(settings["scenario_target"][1],3,"0")*".nc"
    quantiles_raw = ncread(fn,"quantiles")
    idx_low = findfirst(quantiles_raw.>0.1699)
    idx_med = findfirst(quantiles_raw.>0.4999)
    idx_high = findfirst(quantiles_raw.>0.9499)
    scn_grid["quantiles"] = [0.17,0.50,0.83]
    scn_tg["quantiles"] = [0.17,0.50,0.83]
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

    # Read fields
    for process ∈ settings["processes"]
        println("  "*process)
        fn = settings["dir_NCA5"] * "gmsl"*lpad(settings["scenario_target"][scn],3,"0")*"/"*process*"_gmsl"*lpad(settings["scenario_target"][scn],3,"0")*".nc"
        if (process == "landwaterstorage") || (process == "verticallandmotion")
            scn_grid[process] = reshape(ncread(fn,"sea_level_change",start=[findfirst(acc_loc),2,1],count=[-1,-1,-1])[:,:,[idx_low,idx_med,idx_high]],(360,181,length(scn_grid["years"]),3));
            scn_tg[process] = ncread(fn,"sea_level_change",start=[1,2,1],count=[sum(acc_loc_tg),-1,-1])[:,:,[idx_low,idx_med,idx_high]];
        else
            scn_grid[process] = reshape(ncread(fn,"sea_level_change",start=[findfirst(acc_loc),1,1],count=[-1,-1,-1])[:,:,[idx_low,idx_med,idx_high]],(360,181,length(scn_grid["years"]),3));
            scn_tg[process]  = ncread(fn,"sea_level_change",start=[1,1,1],count=[sum(acc_loc_tg),-1,-1])[:,:,[idx_low,idx_med,idx_high]];
        end
        scn_grid[process] = convert.(Float32,scn_grid[process]);
        scn_tg[process] = convert.(Float32,scn_tg[process]);
        reverse!(scn_grid[process],dims=2);
    end

    # Estimate climate-related scenarios
    println("  total_climate")

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

    # Fix NaN
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
    println("  saving")
    fn_grid_out = settings["dir_proj"] * "NCA5_"*settings["scenario_name"][scn]*"_grid.nc"
    fn_tg_out = settings["dir_proj"] * "NCA5_"*settings["scenario_name"][scn]*"_tg.nc"

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

    return nothing
end