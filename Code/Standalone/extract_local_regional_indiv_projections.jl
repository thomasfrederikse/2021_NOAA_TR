dir_code = homedir()*"/Projects/2021_NOAA_TR/Code/"
include(dir_code*"Masks.jl")
include(dir_code*"ComputeRegionalObs.jl")
include(dir_code*"RegionalProjections.jl")
using NetCDF
using DelimitedFiles
using Plots
using Statistics
using NCDatasets
using MAT
using Interpolations
using CSV

# Read individual component projections and get local/regional projections
function ProcessComponentProjections(settings)
    processes = ["AIS","GIS","glaciers","landwaterstorage","oceandynamics","total","verticallandmotion","total_climate"]

    # Get TG locations
    tg_monthly = ComputeRegionalObs.ReadTGData(settings)
    mask = Masks.ReadMask(settings) 

    # Allocate arrays
    projections_local = zeros(Float32,length(tg_monthly["station_names"]),length(settings["NCA5_scenarios"]),length(processes),29,3);
    projections_region = zeros(Float32,length(settings["regions"]),length(settings["NCA5_scenarios"]),length(processes),29,3);

    # Masks on scenario grid and sampling points for TG stations
    fn = settings["dir_component_grid"]*"NCA5_Int_grid.nc"
    lon = ncread(fn,"lon")
    lat = ncread(fn,"lat")
    years = ncread(fn,"years")
    sample_tg_idx = zeros(Int32,length(tg_monthly["station_names"]),2)
    for station ∈ 1:length(tg_monthly["station_names"])
        sample_tg_idx[station,1] = argmin(abs.(lon.-tg_monthly["station_coords"][station,1]))
        sample_tg_idx[station,2] = argmin(abs.(lat.-tg_monthly["station_coords"][station,2]))
    end

    mask_int = zeros(Bool,length(lon),length(lat),length(settings["regions"]));
    wgt_int = zeros(Float32,length(lon),length(lat),length(settings["regions"]));
    area = RegionalProjections.ComputeGridArea(lon,lat)

    for (idx,region) ∈ enumerate(settings["regions"])
        mask_int[:,:,idx] .= LinearInterpolation((mask["ϕ"], mask["θ"]),convert.(Float32,mask[region]),extrapolation_bc = Line())(lon,lat) .> 0.5;
        wgt_int_lcl = zeros(Float32,length(lon),length(lat))
        wgt_int_lcl[mask_int[:,:,idx]] .= area[mask_int[:,:,idx]] ./ sum(area[mask_int[:,:,idx]])
        wgt_int[:,:,idx] = wgt_int_lcl
    end

    for (scn_idx,scenario) ∈ enumerate(settings["NCA5_scenarios"]) 
        fn = settings["dir_component_grid"]*"NCA5_"*scenario*"_grid.nc"
        for (prc_idx,process) ∈ enumerate(processes) 
            data_lcl = ncread(fn,process);
            for station ∈ 1:length(tg_monthly["station_names"])
                @. projections_local[station,scn_idx,prc_idx,:,:] = data_lcl[sample_tg_idx[station,1],sample_tg_idx[station,2],:,:]
            end
            for (region_idx,region) ∈ enumerate(settings["regions"])
                 projections_region[region_idx,scn_idx,prc_idx,:,:] =  sum(data_lcl .* reshape(wgt_int[:,:,region_idx],(length(lon),length(lat),1,1)),dims=(1,2))[1,1,:,:]
            end
        end
    end

    # Save

    fn_lcl = settings["dir_component_grid"]*"projections_tg_locs.nc"
    fh = Dataset(fn_lcl,"c")
    defDim(fh,"tg", length(tg_monthly["station_names"]))
    defDim(fh,"scenario", length(settings["NCA5_scenarios"]))
    defDim(fh,"process", length(processes))
    defDim(fh,"year", length(years))
    defDim(fh,"percentiles", 3)


    defVar(fh,"tg",tg_monthly["station_names"],("tg",),deflatelevel=5)
    defVar(fh,"lon",tg_monthly["station_coords"][:,1],("tg",),deflatelevel=5)
    defVar(fh,"lat",tg_monthly["station_coords"][:,2],("tg",),deflatelevel=5)
    defVar(fh,"scenario",settings["NCA5_scenarios"],("scenario",),deflatelevel=5)
    defVar(fh,"process",processes,("process",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)
    defVar(fh,"years",years,("years",),deflatelevel=5)

    defVar(fh,"projection_local",projections_local,("tg","scenario","process","year","percentiles"),deflatelevel=5)
    close(fh)

    fn_reg = settings["dir_component_grid"]*"projections_regions.nc"
    fh = Dataset(fn_reg,"c")
    defDim(fh,"region", length(settings["regions"]))
    defDim(fh,"scenario", length(settings["NCA5_scenarios"]))
    defDim(fh,"process", length(processes))
    defDim(fh,"year", length(years))
    defDim(fh,"percentiles", 3)

    defVar(fh,"region",settings["regions"],("region",),deflatelevel=5)
    defVar(fh,"scenario",settings["NCA5_scenarios"],("scenario",),deflatelevel=5)
    defVar(fh,"process",processes,("process",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)
    defVar(fh,"years",years,("years",),deflatelevel=5)
    defVar(fh,"projection_regional",projections_region,("region","scenario","process","year","percentiles"),deflatelevel=5)
    close(fh)
    end











