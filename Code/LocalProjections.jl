module LocalProjections
# -------------------------------------------
# Compute local extrapolation and projections
# for individual tide-gauge locations
# -------------------------------------------
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
include(dir_code*"ComputeRegionalObs.jl")
include(dir_code*"RegionalProjections.jl")

function RunLocalProjections(settings)
    local_projections = ReadTGData(settings)
    ExtrapolateLocalTrajectory!(local_projections,settings)
    local_NCA5, years_NCA5, pct_NCA5 = ReadLocalProjections(local_projections,settings)
    save_data(local_projections,local_NCA5, years_NCA5, pct_NCA5,settings)
end

function ReadTGData(settings)
    tg_monthly = ComputeRegionalObs.ReadTGData(settings)
    tg_annual  = ComputeRegionalObs.ConvertMonthlyToAnnual(tg_monthly,settings)
    local_projections = Array{Dict}(undef,length(tg_annual["station_names"]))
    for tg ∈ eachindex(local_projections)
        local_projections[tg] = Dict()
        local_projections[tg]["station_name"] = tg_annual["station_names"][tg]
        local_projections[tg]["station_coords"] = tg_annual["station_coords"][tg,:]
        local_projections[tg]["η_tg"] = LinearInterpolation(tg_annual["years"], tg_annual["rsl"][tg,:],extrapolation_bc=NaN32)(settings["years"])
    end
    return local_projections
end

function ExtrapolateLocalTrajectory!(local_projections,settings)
    println("  Extrapolating tide-gauge records to compute trajectory...")
    # Prepare design matrices
    years_estimate = intersect(settings["years_trajectory"],settings["years_tg"])
    y_acc = findall(in(years_estimate),settings["years"])
    amat = ones(length(y_acc),3)
    amat[:,2] = settings["years"][y_acc] .- mean(settings["years"][y_acc])
    amat[:,3] = (settings["years"][y_acc] .- mean(settings["years"][y_acc])).^2
    amat_extend = ones(length(settings["years_trajectory"]),3)
    amat_extend[:,2] .= settings["years_trajectory"] .- mean(settings["years"][y_acc])
    amat_extend[:,3] = (settings["years_trajectory"] .- mean(settings["years"][y_acc])).^2

    # Do the extrapolation
    for tg ∈ eachindex(local_projections)
        trend_file = Hector.EstTrend(settings["years"][y_acc],local_projections[tg]["η_tg"][y_acc];accel=true,model="Powerlaw",SA=false,SSA=false,monthly=false)
        μ_sol = [trend_file["bias"],trend_file["trend"],trend_file["accel"]]
        σ_sol= [trend_file["bias_sigma"],trend_file["trend_sigma"],trend_file["accel_sigma"]]
        sol_arr = reshape(μ_sol,(1,3)) .+ randn(5000,3) .* reshape(σ_sol,(1,3))
        trajectory_arr = zeros(Float32,length(settings["years_trajectory"]),5000)
        [trajectory_arr[:,i] = amat_extend * @views sol_arr[i,:] for i ∈ 1:5000]
        trajectory_stats = zeros(Float32,length(settings["years_trajectory"]),3)
        [trajectory_stats[t,:] = quantile((@views trajectory_arr[t,:]),[0.17,0.50,0.83]) for t ∈ 1:length(settings["years_trajectory"])]
        local_projections[tg]["η_trajectory"] = LinearInterpolation((settings["years_trajectory"],[1.0f0:3.0f0...]),trajectory_stats,extrapolation_bc=NaN32)[settings["years"],[1.0f0:3.0f0...]]
        local_projections[tg]["trend"] = [μ_sol[2]-σ_sol[2],μ_sol[2],μ_sol[2]+σ_sol[2]]
        local_projections[tg]["accel"] = 2 .*[μ_sol[3]-σ_sol[3],μ_sol[3],μ_sol[3]+σ_sol[3]]
    end
    return nothing
end

function ReadLocalProjections(local_projections,settings)
    # Locations and time
    fn = settings["dir_NCA5"]*"NCA5_Low_grid.nc"
    lon_NCA5 = ncread(fn,"lon")
    lat_NCA5 = ncread(fn,"lat")
    years_NCA5 = convert.(Float32,ncread(fn,"years",start=[1],count=[9]))
    pct_NCA5 = convert.(Float32,ncread(fn,"percentiles"))

    # Find nearest grid cell
    NCA5_loc = zeros(Int,length(local_projections),2)
    for tg in eachindex(local_projections)
        NCA5_loc[tg,1] = argmin(@. abs(local_projections[tg]["station_coords"][1]-lon_NCA5))
        NCA5_loc[tg,2] = argmin(@. abs(local_projections[tg]["station_coords"][2]-lat_NCA5))
    end

    # Create data sctructure
    local_NCA5 = Array{Dict}(undef,length(tg_annual["station_names"]))
    for tg in eachindex(local_projections)
        local_NCA5[tg] = Dict()
        for scn in settings["NCA5_scenarios"]
            local_NCA5[tg][scn] = Dict()
        end
    end

    # Read data
    for scn in settings["NCA5_scenarios"]
        println("   Scenario "*scn*"...")
        fn = settings["dir_NCA5"]*"NCA5_"*scn*"_grid.nc"
        for prc in settings["processes"]
            NCA5_prc = convert.(Float32,ncread(fn,prc,start=[1,1,1,1],count=[-1,-1,9,-1]));
            for tg in eachindex(local_projections)
                local_NCA5[tg][scn][prc] = NCA5_prc[NCA5_loc[tg,1],NCA5_loc[tg,2],:,:]
            end
        end
    end

    # For each scenario, match 2020 value with trajectory for:
    #  total
    #  vlm
    trajectory_idx = findfirst(settings["years"] .== years_NCA5[1])
    for tg in eachindex(local_projections)
        traj_value = local_projections[tg]["η_trajectory"][trajectory_idx,2]
        for scn in settings["NCA5_scenarios"]
            diff_value = local_NCA5[tg][scn]["total"][1,2] - traj_value
            local_NCA5[tg][scn]["total"] .-= diff_value
            local_NCA5[tg][scn]["verticallandmotion"] .-= diff_value
        end
    end
    return local_NCA5,years_NCA5,pct_NCA5
end

function save_data(local_projections,local_NCA5, years_NCA5, pct_NCA5,settings)
    println("  Saving data...")
    station_names = Array{String}(undef,length(local_projections))
    station_coords = Array{Float32}(undef,length(local_projections),2)
    η_tg = Array{Float32}(undef,length(local_projections),length(settings["years"]))
    η_trajectory = Array{Float32}(undef,length(local_projections),length(settings["years"]),3)

    trends = Array{Float32}(undef,length(local_projections),3)
    accels = Array{Float32}(undef,length(local_projections),3)

    scn_proj = Dict()
    for scn ∈ settings["NCA5_scenarios"]
        scn_proj[scn] = Dict()
        for prc in settings["processes"]
            scn_proj[scn][prc] = Array{Float32}(undef,length(local_projections),length(settings["years"]),3)
        end
    end

    for tg ∈ eachindex(local_projections)
        station_names[tg] = local_projections[tg]["station_name"]
        station_coords[tg,:] = local_projections[tg]["station_coords"]
        η_tg[tg,:] = local_projections[tg]["η_tg"]
        η_trajectory[tg,:,:] = local_projections[tg]["η_trajectory"]
        for scn ∈ settings["NCA5_scenarios"]
            for prc in settings["processes"]
                scn_proj[scn][prc][tg,:,:] = LinearInterpolation((years_NCA5, pct_NCA5), local_NCA5[tg][scn][prc],extrapolation_bc=NaN32)(settings["years"],pct_NCA5)
            end
        end
        trends[tg,:] = local_projections[tg]["trend"]
        accels[tg,:] = local_projections[tg]["accel"]
    end
    fh = Dataset(settings["fn_proj_lcl"],"c")
    defDim(fh,"years", length(settings["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"tg", length(local_projections))

    defVar(fh,"tg",station_names,("tg",),deflatelevel=5)
    defVar(fh,"lon",station_coords[:,1],("tg",),deflatelevel=5)
    defVar(fh,"lat",station_coords[:,2],("tg",),deflatelevel=5)
    defVar(fh,"years",settings["years"],("years",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)

    # Write trajectory
    defVar(fh,"MSL_Observed",η_tg,("tg","years"),deflatelevel=5)
    defVar(fh,"MSL_Trajectory",η_trajectory,("tg","years","percentiles"),deflatelevel=5)
    defVar(fh,"MSL_trend",trends,("tg","percentiles"),deflatelevel=5)
    defVar(fh,"MSL_accel",accels,("tg","percentiles"),deflatelevel=5)
    # Write scenarios
    for scn ∈ settings["NCA5_scenarios"]
        for prc in settings["processes"]
            defVar(fh,"MSL_"*prc*"_"*scn,scn_proj[scn][prc],("tg","years","percentiles",),deflatelevel=5)
        end
    end
    close(fh)
    return nothing
end

end
