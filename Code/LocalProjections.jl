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
    tg_monthly = ComputeRegionalObs.ReadTGData(settings)
    tg_annual = ComputeRegionalObs.ConvertMonthlyToAnnual(tg_monthly,settings)

    local_projections = Array{Dict}(undef,length(tg_annual["station_names"]))
    for tg ∈ eachindex(local_projections)
        local_projections[tg] = Dict()
        local_projections[tg]["station_name"] = tg_annual["station_names"][tg]
        local_projections[tg]["station_coords"] = tg_annual["station_coords"][tg,:]
        local_projections[tg]["η_tg"] = LinearInterpolation(tg_annual["years"], tg_annual["rsl"][tg,:],extrapolation_bc=NaN32)(settings["years"])
    end

    # 1. Extrapolate
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
    end

    # 2. Read projections
    println("  Reading NCA5 projections...")
    ϕ,θ,years,percentiles,NCA_grid = RegionalProjections.ReadNCA5(settings)
    baseline_idx = settings["years"] .== years[1]
    for tg ∈ eachindex(local_projections)
        ϕ_idx = argmin(@. abs(local_projections[tg]["station_coords"][1]-ϕ))
        θ_idx = argmin(@. abs(local_projections[tg]["station_coords"][2]-θ))
        local_projections[tg]["η_projection"] = Dict()
        for scenario ∈ settings["NCA5_scenarios"]
            local_projections[tg]["η_projection"][scenario] = LinearInterpolation((years,[1.0f0:3.0f0...]),NCA_grid[scenario][ϕ_idx,θ_idx,:,:],extrapolation_bc=NaN32)[settings["years"],[1.0f0:3.0f0...]]
            local_projections[tg]["η_projection"][scenario] = local_projections[tg]["η_projection"][scenario] .- local_projections[tg]["η_projection"][scenario][baseline_idx,2] .+ local_projections[tg]["η_trajectory"][baseline_idx,2]
        end
    end

    # 3. Save
    println("  Saving data...")
    station_names = Array{String}(undef,length(local_projections))
    station_coords = Array{Float32}(undef,length(local_projections),2)
    η_tg = Array{Float32}(undef,length(local_projections),length(settings["years"]))
    η_trajectory = Array{Float32}(undef,length(local_projections),length(settings["years"]),3)
    scn_proj = Dict()
    [scn_proj[scenario] = Array{Float32}(undef,length(local_projections),length(settings["years"]),3) for scenario ∈ settings["NCA5_scenarios"]]
    for tg ∈ eachindex(local_projections)
        station_names[tg] = local_projections[tg]["station_name"]
        station_coords[tg,:] = local_projections[tg]["station_coords"]
        η_tg[tg,:] = local_projections[tg]["η_tg"]
        η_trajectory[tg,:,:] = local_projections[tg]["η_trajectory"]
        [scn_proj[scenario][tg,:,:] = local_projections[tg]["η_projection"][scenario] for scenario ∈ settings["NCA5_scenarios"]]
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
    # Write scenarios
    for scenario ∈ settings["NCA5_scenarios"]
        defVar(fh,"MSL_"*scenario,scn_proj[scenario],("tg","years","percentiles",),deflatelevel=5)
    end
    close(fh)
end






