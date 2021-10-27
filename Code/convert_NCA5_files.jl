# ---------------------------------------------
# Read all the NCA5 files and extract the grids
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
    scn_data = Dict()

    # Prepare
    fn = settings["dir_NCA5"] * "gmsl"*lpad(settings["scenario_target"][1],3,"0")*"/AIS_gmsl"*lpad(settings["scenario_target"][1],3,"0")*".nc"
    scn_data["years"] = ncread(fn,"years")
    quantiles_raw = ncread(fn,"quantiles")
    idx_low = findfirst(quantiles_raw.>0.1699)
    idx_med = findfirst(quantiles_raw.>0.4999)
    idx_high = findfirst(quantiles_raw.>0.9499)
    scn_data["quantiles"] = [0.17,0.50,0.83]
    acc_loc = ncread(fn,"locations") .> 4000
    latmat = reshape(ncread(fn,"lat")[acc_loc],(360,181))
    lonmat = reshape(ncread(fn,"lon")[acc_loc],(360,181))
    @. lonmat = mod(lonmat,360)
    scn_data["lon"] = lonmat[:,1]
    scn_data["lat"] = reverse(latmat[1,:])

    # Read fields
    for process ∈ settings["processes"]
        println("  "*process)
        fn = settings["dir_NCA5"] * "gmsl"*lpad(settings["scenario_target"][scn],3,"0")*"/"*process*"_gmsl"*lpad(settings["scenario_target"][scn],3,"0")*".nc"
        if (process == "landwaterstorage") || (process == "verticallandmotion")
            scn_data[process] = reshape(ncread(fn,"sea_level_change",start=[findfirst(acc_loc),2,1],count=[-1,-1,-1])[:,:,[idx_low,idx_med,idx_high]],(360,181,length(years),3));
        else
            scn_data[process] = reshape(ncread(fn,"sea_level_change",start=[findfirst(acc_loc),1,1],count=[-1,-1,-1])[:,:,[idx_low,idx_med,idx_high]],(360,181,length(years),3));
        end
        scn_data[process] = convert.(Float32,scn_data[process])
        reverse!(scn_data[process],dims=2)
    end

    # Estimate climate-related scenarios
    println("  total_climate")

    scn_data["total_climate"] = zeros(Float32,size(scn_data["total"]))
    for process ∈ settings["processes_clim"]
        @. scn_data["total_climate"][:,:,:,2] += scn_data[process][:,:,:,2]
        @. scn_data["total_climate"][:,:,:,1] += (abs((scn_data[process][:,:,:,1] -  scn_data[process][:,:,:,2]))^2)
        @. scn_data["total_climate"][:,:,:,3] += (abs((scn_data[process][:,:,:,3] -  scn_data[process][:,:,:,2]))^2)
    end
    @. scn_data["total_climate"][:,:,:,1] = scn_data["total_climate"][:,:,:,2] - sqrt(@. scn_data["total_climate"][:,:,:,1])
    @. scn_data["total_climate"][:,:,:,3] = scn_data["total_climate"][:,:,:,2] + sqrt(@. scn_data["total_climate"][:,:,:,3])
    
    # Save scenario file
    println("  saving")
    @. scn_data["total_climate"][scn_data["total_climate"]<-32000] = -32000 
    @. scn_data["total_climate"][scn_data["total_climate"]>32000] = -32000 
    for process ∈ settings["processes"]
        @. scn_data[process][scn_data[process]<-32000] = -32000 
        @. scn_data[process][scn_data[process]>32000] = -32000 
    end
    fn_out = settings["dir_proj"] * "NCA5_"*settings["scenario_name"][scn]*"_grid.nc"

    fh = Dataset(fn_out,"c")
    defDim(fh,"years", length(scn_data["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"lon", length(scn_data["lon"]))
    defDim(fh,"lat", length(scn_data["lat"]))
    defVar(fh,"lon",scn_data["lon"],("lon",),deflatelevel=6)
    defVar(fh,"lat",scn_data["lat"],("lat",),deflatelevel=6)
    defVar(fh,"years",trunc.(Int16,scn_data["years"]),("years",),deflatelevel=6)
    defVar(fh,"percentiles",trunc.(Int16,[17,50,83]),("percentiles",),deflatelevel=6)
    for process ∈ settings["processes"]
        defVar(fh,process,trunc.(Int16,scn_data[process]),("lon","lat","years","percentiles",),deflatelevel=6)
    end
    defVar(fh,"total_climate",trunc.(Int16,scn_data["total_climate"]),("lon","lat","years","percentiles",),deflatelevel=6)
    close(fh)
    return nothing
end