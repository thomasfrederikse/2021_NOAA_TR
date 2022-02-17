# -------------------------------------------
# LocalProjections:
# For each individual tide-gauge location:
# - Read local observations
# - Compute local extrapolation 
# - Read local NCA5 projections
# - Put everything on a year 2000 baseline
# Save all to NOAA_TR_local_projections.nc
# -------------------------------------------

module LocalProjections

using NetCDF
using DelimitedFiles
using Plots
using Statistics
using NCDatasets
using MAT
using Interpolations
using CSV

dir_code = homedir()*"/Projects/2021_NOAA_TR/Code/"
include(dir_code*"Masks.jl")
include(dir_code*"Hector.jl")
include(dir_code*"RegionalProjections.jl")
include(dir_code*"ProcessObservations.jl")

function RunLocalProjections(settings)
    println("\nLocal observations, trajectories and scenarios...")
    local_obs  = ProcessObservations.ReadLocalObs(settings)
    ExtrapolateLocalTrajectory!(local_obs,settings)
    NCA5_local,years_NCA5,pct_NCA5 = ReadLocalProjections(local_obs,settings)
    correct_baselines!(local_obs,NCA5_local, settings)
    save_data(local_obs,NCA5_local, years_NCA5, pct_NCA5,settings)
    println("Local observations, trajectories and scenarios done\n")
end

function ExtrapolateLocalTrajectory!(local_obs,settings)
    # Compute and extrapolate the trajectory in observed sea level
    println("  Extrapolating tide-gauge records to compute trajectory...")
    # Prepare design matrices
    y_acc = findall(in(settings["years_trajectory"]),settings["years_tg"])
    amat = ones(length(y_acc),3)
    amat[:,2] = settings["years_tg"][y_acc] .- mean(settings["years_tg"][y_acc])
    amat[:,3] = (settings["years_tg"][y_acc] .- mean(settings["years_tg"][y_acc])).^2
    amat_extend = ones(length(settings["years_trajectory"]),3)
    amat_extend[:,2] .= settings["years_trajectory"] .- mean(settings["years_tg"][y_acc])
    amat_extend[:,3] = (settings["years_trajectory"] .- mean(settings["years_tg"][y_acc])).^2

    # Do the extrapolation
    local_obs["obs_length"] = zeros(Int32,length(local_obs["name"])); # Number of observation-years used to estimate trajectory
    local_obs["obs_lt_30"] = zeros(Int32,length(local_obs["name"])); # Number of observation-years used to estimate trajectory
    local_obs["rsl_trend"] = zeros(Float32,length(local_obs["name"]),3);
    local_obs["rsl_accel"] = zeros(Float32,length(local_obs["name"]),3);
    local_obs["rsl_trajectory"] = zeros(Float32,length(local_obs["name"]),length(settings["years"]),3);
    for tg in 1:length(local_obs["name"])
        local_obs["obs_length"][tg] = sum(isfinite.(local_obs["rsl"][tg,y_acc]))
        local_obs["obs_length"][tg] >= 30 ? local_obs["obs_lt_30"][tg] = 1 : local_obs["obs_lt_30"][tg] = 0
        trend_file = Hector.EstTrend(settings["years_tg"][y_acc],local_obs["rsl"][tg,y_acc];accel=true,model="Powerlaw",SA=false,SSA=false,monthly=false,tref=mean(settings["years_tg"][y_acc]))
        μ_sol = [trend_file["bias"],trend_file["trend"],trend_file["accel"]]
        σ_sol= [trend_file["bias_sigma"],trend_file["trend_sigma"],trend_file["accel_sigma"]]
        sol_arr = reshape(μ_sol,(1,3)) .+ randn(5000,3) .* reshape(σ_sol,(1,3))
        trajectory_arr = zeros(Float32,length(settings["years_trajectory"]),5000)
        [trajectory_arr[:,i] = amat_extend * @views sol_arr[i,:] for i ∈ 1:5000]
        trajectory_stats = zeros(Float32,length(settings["years_trajectory"]),3)
        [trajectory_stats[t,:] = quantile((@views trajectory_arr[t,:]),[0.17,0.50,0.83]) for t ∈ 1:length(settings["years_trajectory"])]
        local_obs["rsl_trajectory"][tg,:,:] = LinearInterpolation((settings["years_trajectory"],[1.0f0:3.0f0...]),trajectory_stats,extrapolation_bc=NaN32)[settings["years"],[1.0f0:3.0f0...]]
        local_obs["rsl_trend"][tg,:] = [μ_sol[2]-σ_sol[2],μ_sol[2],μ_sol[2]+σ_sol[2]]
        local_obs["rsl_accel"][tg,:] = 2 .*[μ_sol[3]-σ_sol[3],μ_sol[3],μ_sol[3]+σ_sol[3]]
    end
    return nothing
end

function ReadLocalProjections(local_obs,settings)
    # Locations and time
    fn_grid = settings["dir_NCA5"]*"NCA5_Low_grid.nc"
    lon_NCA5 = ncread(fn_grid,"lon")
    lat_NCA5 = ncread(fn_grid,"lat")
    years_NCA5 = convert.(Float32,ncread(fn_grid,"years",start=[1],count=[-1]))
    pct_NCA5 = convert.(Float32,ncread(fn_grid,"percentiles"))

    fh_tg = Dataset(settings["dir_NCA5"]*"NCA5_Low_tg.nc","r")
    psmsl_loc = fh_tg["loc"][:]
    close(fh_tg)

    NCA5_loc_is_present = zeros(Bool,length(local_obs["name"]))
    NCA5_loc_tg = zeros(Int32,length(local_obs["name"]))
    NCA5_loc_grid = zeros(Int,length(local_obs["name"]),2)
    for tg in 1:length(local_obs["name"])
        loc_idx = findfirst(==(local_obs["psmsl_id"][tg]),psmsl_loc)
        if ~isnothing(loc_idx)
            # Tide gauge is in PSMSL/local database: use tide-gauge projection
            NCA5_loc_is_present[tg] = true
            NCA5_loc_tg[tg] = loc_idx
            NCA5_loc_grid[tg,1] = argmin(@. abs(local_obs["coords"][tg,1]-lon_NCA5))
            NCA5_loc_grid[tg,2] = argmin(@. abs(local_obs["coords"][tg,2]-lat_NCA5))
        else
            # Tide gauge is not in PSMSL/local database: sample location from the grid
            NCA5_loc_is_present[tg] = false
            NCA5_loc_tg[tg] = -1
            NCA5_loc_grid[tg,1] = argmin(@. abs(local_obs["coords"][tg,1]-lon_NCA5))
            NCA5_loc_grid[tg,2] = argmin(@. abs(local_obs["coords"][tg,2]-lat_NCA5))
        end
    end

    # Create data sctructure
    NCA5_local = Array{Dict}(undef,length(local_obs["name"]))
    for tg in 1:length(local_obs["name"])
        NCA5_local[tg] = Dict()
        NCA5_local[tg]["Loc_is_in_proj"] = NCA5_loc_is_present[tg]
        for scn in settings["NCA5_scenarios"]
            NCA5_local[tg][scn] = Dict()
        end
    end

    # Read data
    for scenario in settings["NCA5_scenarios"]
        println("   Scenario "*scenario*"...")
        fn_grid = settings["dir_NCA5"]*"NCA5_"*scenario*"_grid.nc"
        fn_tg = settings["dir_NCA5"]*"NCA5_"*scenario*"_tg.nc"
        for prc in settings["processes"]
            NCA5_prc_grid = convert.(Float32,ncread(fn_grid,prc));
            NCA5_prc_tg = convert.(Float32,ncread(fn_tg,prc));
            for tg in 1:length(local_obs["name"])
                if NCA5_loc_is_present[tg] # Use TG scenario
                    NCA5_local[tg][scenario][prc] = NCA5_prc_tg[NCA5_loc_tg[tg],:,:]
                else # Use gridded scenario
                    NCA5_local[tg][scenario][prc] = NCA5_prc_grid[NCA5_loc_grid[tg,1],NCA5_loc_grid[tg,2],:,:]
                end
            end
        end
    end
    return NCA5_local,years_NCA5,pct_NCA5
end

function correct_baselines!(local_obs,NCA5_local, settings)
    println("  Correcting baselines for trajectory, observations and scenarios...")
    # Correct baseline: baseline for everything is trajectory value in 2000
    traj_baseline = findall(in(settings["years_baseline"]),settings["years"])
    traj_nca5_start = findfirst(settings["years"].==2005)

    for tg in 1:length(local_obs["name"])
        # 1. Baseline trajectory
        local_obs["rsl_trajectory"][tg,:,:] .-= local_obs["rsl_trajectory"][tg,traj_baseline,2]
        
        # 2. Baseline observations 
        acc_obs = isfinite.(local_obs["rsl"][tg,:])
        yr_int = intersect(settings["years_tg"][acc_obs],settings["years_trajectory"])
        obs_mn = mean(local_obs["rsl"][tg,(findall(in(yr_int),settings["years_tg"]))])
        traj_mn = mean(local_obs["rsl_trajectory"][tg,(findall(in(yr_int),settings["years"])),2])
        local_obs["rsl"][tg,:] = local_obs["rsl"][tg,:] .- obs_mn .+ traj_mn

        # 3. Baseline scenarios
        for scenario in settings["NCA5_scenarios"]
            NCA5_local[tg][scenario]["total"] .+= local_obs["rsl_trajectory"][tg,traj_nca5_start,2]
            NCA5_local[tg][scenario]["verticallandmotion"] .+= local_obs["rsl_trajectory"][tg,traj_nca5_start,2]
        end
    end
    return nothing
end

function save_data(local_obs,NCA5_local, years_NCA5, pct_NCA5,settings)
    println("  Saving data...");
    scn_proj = Dict()
    for scenario ∈ settings["NCA5_scenarios"]
        scn_proj[scenario] = Dict()
        for prc in settings["processes"]
            scn_proj[scenario][prc] = Array{Float32}(undef,length(local_obs["name"]),length(settings["years"]),3)
            for tg in 1:length(local_obs["name"])
                scn_proj[scenario][prc][tg,:,:] = LinearInterpolation((years_NCA5, pct_NCA5), NCA5_local[tg][scenario][prc],extrapolation_bc=NaN32)(settings["years"],pct_NCA5)
            end
        end
    end

    rsl_obs = zeros(Float32,length(local_obs["name"]),length(settings["years"])) .* NaN32
    acc_obs = findall(in(settings["years_tg"]),settings["years"])
    rsl_obs[:,acc_obs] = local_obs["rsl"]
   
    fh = Dataset(settings["fn_proj_lcl"],"c")
    fh.attrib["title"] = "Local projections"
    fh.attrib["description"] = "Local projections, trajectory, and observations at tide-gauge locations for the interagency report: Global and Regional Sea Level Rise Scenarios for the United States: Updated Mean Projections and Extreme Water Level Probabilities Along U.S. Coastlines"
    fh.attrib["processes"] = "AIS: Antarctic Ice Sheet, GIS: Greenland, glaciers: Glaciers and Ice Caps, landwaterstorage: Liquid water storage changes on land, oceandynamics: Ocean dynamics and global thermosteric expansion, verticallandmotion: Vertical land motion, total: All processes combined."

    defDim(fh,"years", length(settings["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"tg", length(local_obs["name"]))

    defVar(fh,"tg",local_obs["name"],("tg",),deflatelevel=5, attrib = Dict("description" => "Name of the tide gauge", "units" => "-"))
    defVar(fh,"lon",local_obs["coords"][:,1],("tg",),deflatelevel=5, attrib = Dict("description" => "Longitude of the tide gauge", "units" => "degrees East"))
    defVar(fh,"lat",local_obs["coords"][:,2],("tg",),deflatelevel=5, attrib = Dict("description" => "Latitude of the tide gauge", "units" => "degrees North"))
    defVar(fh,"years",settings["years"],("years",),deflatelevel=5, attrib = Dict("units" => "years"))
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5, attrib = Dict("description" => "The percentile of the projection. The 50th percentile shows the median projection, and the 17th and 83rd show the upper- and lower bound on the 1 sigma level.", "units" => "-"))

    defVar(fh,"PSMSL_id",convert.(Int32,local_obs["psmsl_id"]),("tg",),deflatelevel=5, attrib = Dict("description" => "PSMSL ID to which the tide gauge corresponds. Stations that are not in PSMSL are flagged -1."))
    defVar(fh,"QC_flag",convert.(Int8,local_obs["qc_flag"]),("tg",),deflatelevel=5, attrib = Dict("description" => "Quality control flag. Value of 1 means possible issue. Currently manual process."))
    # Write trajectory
    defVar(fh,"rsl_obs",rsl_obs,("tg","years"),deflatelevel=5, attrib = Dict("description" => "Observed local sea-level changes, relative to year 2000 (tide-gauge observations).", "units" => "mm"))
    defVar(fh,"number_of_observation_years_for_trajectory_estimation",local_obs["obs_length"],("tg",),deflatelevel=5, attrib = Dict("description" => "Number of observation-years used to estimate the trajectory.", "units" => "-"))
    defVar(fh,"number_of_observation_years_for_trajectory_estimation_geq_30",local_obs["obs_lt_30"],("tg",),deflatelevel=5, attrib = Dict("description" => "Flag if number of observation-years used for the trajectory exceeds 30. 1: 30 years or more, 0: less than 30 years.", "units" => "-"))
    defVar(fh,"rsl_trajectory",local_obs["rsl_trajectory"],("tg","years","percentiles"),deflatelevel=5, attrib = Dict("description" => "Estimated local sea-level trajectory, relative to year 2000.", "units" => "mm"))
    defVar(fh,"rsl_trend",local_obs["rsl_trend"],("tg","percentiles"),deflatelevel=5, attrib = Dict("description" => "Linear trend of estimated trajectory", "units" => "mm yr-1"))
    defVar(fh,"rsl_accel",local_obs["rsl_accel"],("tg","percentiles"),deflatelevel=5, attrib = Dict("description" => "Acceleration (half quadratic) of estimated trajectory", "units" => "mm yr-2"))
    # Write scenarios
    for scn ∈ settings["NCA5_scenarios"]
        for prc in settings["processes"]
            defVar(fh,"rsl_"*prc*"_"*scn,scn_proj[scn][prc],("tg","years","percentiles",),deflatelevel=5, attrib = Dict("description" => "Projected local sea level for scenario "*scn*" and process "*prc*", relative to year 2000.", "units" => "mm"))
        end
    end
    close(fh)
    return nothing
end

function example_plots(local_obs,NCA5_local, years_NCA5, pct_NCA5,settings)
    # Make a simple plot of all local time series 
    for tg in 1:length(local_obs["name"])
        scn_cl = cgrad(:viridis, 3, categorical = true)
        plot(settings["years_tg"],local_obs["rsl"][tg,:],color=:red,linewidth=2)
        plot!(settings["years"],local_obs["rsl_trajectory"][tg,:,2],color=:blue,linewidth=2)
        plot!(years_NCA5,NCA5_local[tg]["Low"]["total"][:,2],color=scn_cl[1],linewidth=2)
        plot!(years_NCA5,NCA5_local[tg]["Int"]["total"][:,2],color=scn_cl[2],linewidth=2)
        plot!(years_NCA5,NCA5_local[tg]["High"]["total"][:,2],color=scn_cl[3],linewidth=2,legend=false)
        xlims!((1970,2050))
        ylims!((-200,1000))
        title!(local_obs["name"][tg])
        savefig(homedir()*"/Scratch/Plots/"*string(tg)*".png")
    end
end

end