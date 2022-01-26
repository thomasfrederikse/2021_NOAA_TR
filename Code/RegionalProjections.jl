# --------------------------------------------
# Process regional projections and observations
# - Read TG estimates
# - Estimate trend + acceleration using Hector
# - Extrapolate data from tide gauges
# - Read NCA5 projections
# --------------------------------------------

module RegionalProjections
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
include(dir_code*"ProcessObservations.jl")
function RunRegionalProjections(settings)
    println("\nRegional observations, trajectories and scenarios...")
    local_obs  = ProcessObservations.ReadLocalObs(settings)
    region_obs = ProcessObservations.ComputeRegionObs(local_obs,settings)
    ComputeRegionalTrajectory!(region_obs,settings)
    NCA5_regional = ReadNCA5Regional(settings)
    correct_baselines!(region_obs,NCA5_regional,settings)
    SaveRegionData(region_obs,NCA5_regional,settings)
    println("Regional observations, trajectories and scenarios done\n")
    return nothing
end

function ComputeRegionalTrajectory!(region_obs,settings)
    # -----------------------------------------------
    # Compute the trajectory for regional sea level 
    # and extrapolate
    # -----------------------------------------------

    println("  Computing trajectories...")
    years_estimate = intersect(settings["years_trajectory"],settings["years_tg"])
    y_acc = findall(in(years_estimate),settings["years_tg"])

    amat = ones(length(y_acc),3)
    amat[:,2] = settings["years_tg"][y_acc] .- mean(settings["years_tg"][y_acc])
    amat[:,3] = (settings["years_tg"][y_acc] .- mean(settings["years_tg"][y_acc])).^2

    amat_extend = ones(length(settings["years_trajectory"]),3)
    amat_extend[:,2] .= settings["years_trajectory"] .- mean(settings["years_tg"][y_acc])
    amat_extend[:,3] = (settings["years_trajectory"] .- mean(settings["years_tg"][y_acc])).^2
    for region in settings["regions"]
        trend_file = Hector.EstTrend(settings["years_tg"][y_acc],region_obs[region]["rsl"][y_acc];accel=true,model="Powerlaw",SA=false,SSA=false,monthly=false,tref=mean(settings["years_tg"][y_acc]))
        μ_sol = [trend_file["bias"],trend_file["trend"],trend_file["accel"]]
        σ_sol= [trend_file["bias_sigma"],trend_file["trend_sigma"],trend_file["accel_sigma"]]
        sol_arr = reshape(μ_sol,(1,3)) .+ randn(5000,3) .* reshape(σ_sol,(1,3))
        trajectory_arr = zeros(Float32,length(settings["years_trajectory"]),5000)
        [trajectory_arr[:,i] = amat_extend * @views sol_arr[i,:] for i ∈ 1:5000]
        trajectory_stats = zeros(Float32,length(settings["years_trajectory"]),3)
        [trajectory_stats[t,:] = quantile((@views trajectory_arr[t,:]),[0.17,0.50,0.83]) for t ∈ 1:length(settings["years_trajectory"])]
        region_obs[region]["rsl_trajectory"] = trajectory_stats
        region_obs[region]["trend"] = [μ_sol[2]-σ_sol[2],μ_sol[2],μ_sol[2]+σ_sol[2]]
        region_obs[region]["accel"] = 2 .*[μ_sol[3]-σ_sol[3],μ_sol[3],μ_sol[3]+σ_sol[3]]
    end
    return nothing
end

function ReadNCA5Regional(settings)
    println("  Reading NCA5 scenarios...")
    
    # Read mask and grid information
    mask = Masks.ReadMask(settings) 
    area = ComputeGridArea(mask["ϕ"],mask["θ"])

    fn = settings["dir_NCA5"]*"NCA5_Low_grid.nc"
    lon_NCA5 = ncread(fn,"lon")
    lat_NCA5 = ncread(fn,"lat")
    years_NCA5 = convert.(Float32,ncread(fn,"years",start=[1],count=[-1]))
    pct_NCA5 = convert.(Float32,ncread(fn,"percentiles"))

    # Read and interpolate the data
    # Structure NCA5_projections[scenario][process][region]
    # Store in dictionary
    NCA5_regional = Dict()
    NCA5_regional["years"] = years_NCA5
    NCA5_regional["percentiles"] = pct_NCA5
    for scn in settings["NCA5_scenarios"]
        println("   Scenario "*scn*"...")
        NCA5_regional[scn] = Dict()
        fn = settings["dir_NCA5"]*"NCA5_"*scn*"_grid.nc"
        for prc in settings["processes"]
            NCA5_regional[scn][prc] = Dict()
            NCA5_prc_raw = convert.(Float32,ncread(fn,prc,start=[1,1,1,1],count=[-1,-1,-1,-1]));
            NCA5_prc_int = LinearInterpolation((lon_NCA5, lat_NCA5,years_NCA5,pct_NCA5),NCA5_prc_raw,extrapolation_bc = Line())(mask["ϕ"],mask["θ"],years_NCA5,pct_NCA5);
            for region in settings["regions"]
                NCA5_regional[scn][prc][region] = (sum((@. NCA5_prc_int*area*mask[region]),dims=(1,2)) / sum(area .* mask[region]))[1,1,:,:]
            end
        end
    end
    return NCA5_regional
end

function ComputeGridArea(lon,lat)
    gridsize = abs(lat[2]-lat[1])
    area = @. deg2rad(gridsize) * (sind(lat+gridsize/2)-sind(lat-gridsize/2)) * 6371000^2
    return repeat(area',size(lon,1))
end

function correct_baselines!(region_obs,NCA5_regional,settings)
    # Correct baseline: baseline for everything is trajectory value in 2000
    traj_baseline = findall(in(settings["years_baseline"]),settings["years_trajectory"])
    traj_nca5_start = findfirst(settings["years_trajectory"].==2005)

    # TG and altimetry average
    traj_tg_bl = findall(in(intersect(settings["years_trajectory"],settings["years_tg"])),settings["years_trajectory"])
    obs_tg_bl  = findall(in(intersect(settings["years_trajectory"],settings["years_tg"])),settings["years_tg"])

    obs_alt_bl = findall(isfinite,region_obs["USA"]["alt"])
    traj_alt_bl = findall(in(intersect(settings["years_trajectory"],settings["years_tg"][obs_alt_bl])),settings["years_trajectory"])
    for region in settings["regions"]
        region_obs[region]["rsl_trajectory"] .-= mean(region_obs[region]["rsl_trajectory"][traj_baseline,2])
        region_obs[region]["rsl"] = region_obs[region]["rsl"] .- mean(filter(isfinite,region_obs[region]["rsl"][obs_tg_bl])) .+ mean(region_obs[region]["rsl_trajectory"][traj_tg_bl,2])
        region_obs[region]["alt"] = region_obs[region]["alt"] .- mean(region_obs[region]["alt"][obs_alt_bl]) .+ mean(region_obs[region]["rsl_trajectory"][traj_alt_bl,2])
        for scenario in settings["NCA5_scenarios"]
            NCA5_regional[scenario]["total"][region] .+= region_obs[region]["rsl_trajectory"][traj_nca5_start,2]
            NCA5_regional[scenario]["verticallandmotion"][region] .+= region_obs[region]["rsl_trajectory"][traj_nca5_start,2]
        end
    end
    return nothing
end

function SaveRegionData(region_obs,NCA5_regional,settings)
    println("  Saving...")
    fh = Dataset(settings["fn_proj_reg"],"c")
    defDim(fh,"years", length(settings["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"region", length(settings["regions"]))

    defVar(fh,"region",settings["regions"],("region",),deflatelevel=5)
    defVar(fh,"years",settings["years"],("years",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)

    # Write observations
    rsl_obs_array = zeros(Float32,length(settings["regions"]),length(settings["years"])) .* NaN32
    alt_obs_array = zeros(Float32,length(settings["regions"]),length(settings["years"])) .* NaN32

    for (region_idx,region) ∈ enumerate(settings["regions"])
        rsl_obs_array[region_idx,:] = LinearInterpolation((convert.(Float32,settings["years_tg"])),region_obs[region]["rsl"],extrapolation_bc=NaN32)(settings["years"])
        alt_obs_array[region_idx,:] = LinearInterpolation((convert.(Float32,settings["years_tg"])),region_obs[region]["alt"],extrapolation_bc=NaN32)(settings["years"])
    end
    defVar(fh,"rsl_obs",rsl_obs_array,("region","years"),deflatelevel=5)
    defVar(fh,"alt_obs",alt_obs_array,("region","years"),deflatelevel=5)

    # Write trajectory
    trajectory = zeros(Float32,length(settings["regions"]),length(settings["years"]),3) .* NaN32
    [trajectory[region_idx,:,:] = (LinearInterpolation((settings["years_trajectory"],[1.0f0:3.0f0...]),region_obs[region]["rsl_trajectory"],extrapolation_bc=NaN32)[settings["years"],[1.0f0:3.0f0...]]) for (region_idx,region) ∈ enumerate(settings["regions"])]
    defVar(fh,"rsl_trajectory",trajectory,("region","years","percentiles"),deflatelevel=5)

    # Write trends and accelerations
    trends = zeros(Float32,length(settings["regions"]),3) .* NaN32
    accels = zeros(Float32,length(settings["regions"]),3) .* NaN32
    [trends[region_idx,:] = region_obs[region]["trend"] for (region_idx,region) ∈ enumerate(settings["regions"])]
    [accels[region_idx,:] = region_obs[region]["accel"] for (region_idx,region) ∈ enumerate(settings["regions"])]
    defVar(fh,"rsl_trend",trends,("region","percentiles"),deflatelevel=5)
    defVar(fh,"rsl_accel",accels,("region","percentiles"),deflatelevel=5)

    # Write projections
    scn_array = zeros(Float32,length(settings["regions"]),length(settings["years"]),3) .* NaN32
    for prc in settings["processes"] 
        for scn ∈ settings["NCA5_scenarios"]
            for (region_idx,region) ∈ enumerate(settings["regions"])
                scn_array[region_idx,:,:] = (LinearInterpolation((NCA5_regional["years"],[1.0f0:3.0f0...]),NCA5_regional[scn][prc][region],extrapolation_bc=NaN32)[settings["years"],[1.0f0:3.0f0...]])
            end
            defVar(fh,"rsl_"*prc*"_"*scn,scn_array,("region","years","percentiles",),deflatelevel=5)
        end
    end
    close(fh)
end

end