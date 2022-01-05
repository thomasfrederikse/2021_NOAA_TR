# -------------------------------------------
# Process global projections and observations
# - Read altimetry observation (GSFC/MEASURES)
# - Read TG reconstruction ensemble (F2020)
# - Estimate trend + acceleration using Hector
# - Extrapolate data from tide gauges
# - Read NCA5 projections
# --------------------------------------------
module GlobalProjections
using NetCDF
using DelimitedFiles
using Statistics
using NCDatasets
using MAT
using Interpolations
using CSV
dir_code = homedir()*"/Projects/2021_NOAA_TR/Code/"
include(dir_code*"Hector.jl")
function RunGlobalProjections(settings)
    println("\nGlobal-mean observations, trajectories, and scenarios...")
    GMSL_20c = read_GMSL_20c(settings)
    GMSL_Trajectory = compute_GMSL_trajectory(GMSL_20c,settings)
    GMSL_Altimetry  = read_GMSL_Altimetry(settings)
    GMSL_NCA5       = read_NCA5(settings)
    correct_baselines!(GMSL_20c,GMSL_Trajectory,GMSL_Altimetry,GMSL_NCA5,settings)
    save_data(GMSL_20c,GMSL_Trajectory,GMSL_Altimetry,GMSL_NCA5,settings)
    println("Global-mean observations, trajectories, and scenarios done\n")
    return nothing
end

function read_GMSL_20c(settings)
    println("  Reading tide-gauge GMSL...")
    # Read F2020 GMSL
    fh = Dataset(settings["fn_gmsl_20c_ensembles"],"r")
    t = fh["time"][:]
    Λ = fh["likelihood"][:]
    η = fh["GMSL"][:]

    # Adjust to baseline
    idx_baseline = findall(in(settings["years_baseline"]),t)
    η .-= mean(η[idx_baseline,:],dims=1)
    # Mean and standard deviation
    η_mean = zeros(length(t),3)
    sidx = zeros(Int,length(Λ))
    𝚲    = zeros(Float64,length(Λ))
    Η    = zeros(Float64,length(Λ))
    for ts ∈ 1:length(t)
        sortperm!(sidx, η[ts,:])
        cumsum!(𝚲,Λ[sidx])
        @. Η = @views η[ts,sidx];
        η_mean[ts,1] = Η[findfirst(>(0.17),𝚲)]
        η_mean[ts,2] = Η[findfirst(>(0.5),𝚲)]
        η_mean[ts,3] = Η[findfirst(>(0.83),𝚲)]
    end
    GMSL_20c = Dict()
    GMSL_20c["years"] = t
    GMSL_20c["rsl"] = η_mean
    GMSL_20c["likelihood"] = Λ
    GMSL_20c["rsl_ensemble"] = η
    return GMSL_20c
end

function compute_GMSL_trajectory(GMSL_20c,settings)
    # -----------------------------------------------
    # Compute the trajectory for GMSL and extrapolate
    # -----------------------------------------------
    println("  Computing GMSL trajectories...")
    t_acc = findall(in(settings["years_trajectory_global"]),GMSL_20c["years"])

    # Compute trend and acceleration in each of the ensemble members of F2020 GMSL
    amat = ones(length(t_acc),3)
    amat[:,2] = GMSL_20c["years"][t_acc] .- mean(GMSL_20c["years"][t_acc])
    @. amat[:,3] = 0.5 * (GMSL_20c["years"][t_acc] - $mean(GMSL_20c["years"][t_acc]))^2
    amat_tr = transpose(amat)
    amat_sq = inv(amat_tr*amat)
    sol_arr = zeros(length(GMSL_20c["likelihood"]),3)
    amat_expand = ones(length(settings["years_trajectory_global"]),3)
    amat_expand[:,2] = settings["years_trajectory_global"] .- mean(GMSL_20c["years"][t_acc])
    @. amat_expand[:,3] = 0.5 * (settings["years_trajectory_global"] - $mean(GMSL_20c["years"][t_acc]))^2
    η_projected = zeros(length(settings["years_trajectory_global"]),length(GMSL_20c["likelihood"]))

    # Generate autocorrelated noise to account for uncertainties in the trajectory due to
    # internal variability and add this to the ensemble members
    noise = Hector.GenerateNoise(GMSL_20c["years"][t_acc],GMSL_20c["rsl"][t_acc,2],length(GMSL_20c["likelihood"]),length(GMSL_20c["years"]);accel=true,model="Powerlaw");
    GMSL_20c["rsl_ensemble"] .+= noise

    # Compute trend and acceleration for each ens member, and compute trajectory
    # time series: design matrix * solution
    for λ ∈ 1:length(GMSL_20c["likelihood"])
        sol_arr[λ,:] = amat_sq * (amat_tr*GMSL_20c["rsl_ensemble"][t_acc,λ])
        η_projected[:,λ] = amat_expand * sol_arr[λ,:]
    end

    # Statistics: compute mean and confidence intervals
    # of the extrapolated trajectory
    projection_mean = zeros(length(settings["years_trajectory_global"]),3)
    sidx = zeros(Int,length(GMSL_20c["likelihood"]))
    𝚲    = zeros(Float64,length(GMSL_20c["likelihood"]))
    Η    = zeros(Float64,length(GMSL_20c["likelihood"]))
    for ts ∈ 1:length(settings["years_trajectory_global"])
        sortperm!(sidx, η_projected[ts,:])
        cumsum!(𝚲,GMSL_20c["likelihood"][sidx])
        @. Η = @views η_projected[ts,sidx];
        projection_mean[ts,1] = Η[findfirst(>(0.17),𝚲)]
        projection_mean[ts,2] = Η[findfirst(>(0.5),𝚲)]
        projection_mean[ts,3] = Η[findfirst(>(0.83),𝚲)]
    end
    GMSL_Trajectory = Dict()
    GMSL_Trajectory["rsl"] = projection_mean
    GMSL_Trajectory["likelihood"] = GMSL_20c["likelihood"]

    # Statistics: get trend and acceleration of trajectory
    GMSL_Trajectory["trend"] = zeros(Float32,3)
    sidx = sortperm(sol_arr[:,2])
    𝚲 = cumsum(GMSL_20c["likelihood"][sidx])
    Η = sol_arr[sidx,2];
    GMSL_Trajectory["trend"][1] = Η[findfirst(>(0.17),𝚲)]
    GMSL_Trajectory["trend"][2] = Η[findfirst(>(0.5),𝚲)]
    GMSL_Trajectory["trend"][3] = Η[findfirst(>(0.83),𝚲)]

    GMSL_Trajectory["accel"] = zeros(Float32,3)
    sidx = sortperm(sol_arr[:,3])
    𝚲 = cumsum(GMSL_20c["likelihood"][sidx])
    Η = sol_arr[sidx,3];
    GMSL_Trajectory["accel"][1] = Η[findfirst(>(0.17),𝚲)]
    GMSL_Trajectory["accel"][2] = Η[findfirst(>(0.5),𝚲)]
    GMSL_Trajectory["accel"][3] = Η[findfirst(>(0.83),𝚲)]
    return GMSL_Trajectory
end

function read_GMSL_Altimetry(settings)
    println("  Reading altimetry GMSL...")
    # Read Altimetry GMSL from Goddard and apply corrections for
    # internal variability due to ENSO/PDO and take annual-means
    GMSL_Altimetry = Dict()
    GMSL_Altimetry["years"] = [1993:2020...]
      # Read data
    GSFC_raw = readdlm(settings["fn_gmsl_altimetry_GSFC"],skipstart=48)[:,[3,12]]
    ENSO_raw = matread(settings["fn_gmsl_ENSO"])
    # To monthly-mean GMSL 
    tbounds = zeros(length(ENSO_raw["time"]),2)
    tbounds[:,1] = [1993.0:1/12:2021-1/12...]
    tbounds[:,2] = tbounds[:,1] .+ 1/12
    GSFC_monthly = zeros(Float32,length(ENSO_raw["time"]))
    [GSFC_monthly[t] = mean(GSFC_raw[(GSFC_raw[:,1] .>= tbounds[t,1]) .& (GSFC_raw[:,1] .<= tbounds[t,2]),2]) - ENSO_raw["gmsl"][1,t] for t ∈ 1:length(ENSO_raw["time"])]
     GMSL_Altimetry["alt"] =mean(reshape(GSFC_monthly,(12,:)),dims=1)[1,:] # Monthly-mean to annual-mean
    return GMSL_Altimetry
end

function read_NCA5(settings)
    println("  Reading NCA5 scenarios...")
    # Read the GMSL scenarios
    GMSL_NCA5 = Dict()
    GMSL_NCA5["years"] = ncread(settings["dir_NCA5"] * "NCA5_Low_gmsl.nc","years",start=[1],count=[-1])
    for scenario ∈ settings["NCA5_scenarios"]
        fn = settings["dir_NCA5"] * "NCA5_"*scenario*"_gmsl.nc"
        GMSL_NCA5[scenario] = convert.(Float32,ncread(fn,"total",start=[1,1],count=[-1,-1]))    
    end
    return GMSL_NCA5
end

function correct_baselines!(GMSL_20c,GMSL_Trajectory,GMSL_Altimetry,GMSL_NCA5,settings)
    # Correct baseline: baseline for everything is trajectory value in 2000
    traj_baseline = findall(in(settings["years_baseline"]),settings["years_trajectory_global"])
    traj_nca5_start = findfirst(settings["years_trajectory_global"].==2005)

    GMSL_Trajectory["rsl"] .-= mean(GMSL_Trajectory["rsl"][traj_baseline,:])

    traj_bl = findall(in(intersect(settings["years_trajectory_global"],GMSL_20c["years"])),settings["years_trajectory_global"])
    obs_bl  = findall(in(intersect(settings["years_trajectory_global"],GMSL_20c["years"])),GMSL_20c["years"])
    GMSL_20c["rsl"] = GMSL_20c["rsl"] .- mean(GMSL_20c["rsl"][obs_bl,2]) .+ mean(GMSL_Trajectory["rsl"][traj_bl,2])

    traj_bl = findall(in(intersect(settings["years_trajectory_global"],GMSL_Altimetry["years"])),settings["years_trajectory_global"])
    obs_bl  = findall(in(intersect(settings["years_trajectory_global"],GMSL_Altimetry["years"])),GMSL_Altimetry["years"])
    GMSL_Altimetry["alt"] = GMSL_Altimetry["alt"] .- mean(GMSL_Altimetry["alt"][obs_bl]) .+ mean(GMSL_Trajectory["rsl"][traj_bl,2])

    for scenario in settings["NCA5_scenarios"]
        GMSL_NCA5[scenario] .+= GMSL_Trajectory["rsl"][traj_nca5_start,2]
    end
    return nothing
end

function save_data(GMSL_20c,GMSL_Trajectory,GMSL_Altimetry,GMSL_NCA5,settings)
    println("  Saving data...")
    fh = Dataset(settings["fn_proj_glb"],"c")
    defDim(fh,"years", length(settings["years"]))
    defDim(fh,"percentiles",3)
    defVar(fh,"years",settings["years"],("years",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)

    # Write observations and trajectory
    defVar(fh,"GMSL_20c",Float32,("years","percentiles"),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,GMSL_20c["years"]),[1.0f0:3.0f0...]),GMSL_20c["rsl"],extrapolation_bc=NaN32)(settings["years"],[1.0f0:3.0f0...])
    defVar(fh,"GMSL_Altimetry",Float32,("years",),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,GMSL_Altimetry["years"])),GMSL_Altimetry["alt"],extrapolation_bc=NaN32)(settings["years"])
    defVar(fh,"GMSL_Trajectory",Float32,("years","percentiles"),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,settings["years_trajectory_global"]),[1.0f0:3.0f0...]),GMSL_Trajectory["rsl"],extrapolation_bc=NaN32)(settings["years"],[1.0f0:3.0f0...])
    
    defVar(fh,"GMSL_trend",GMSL_Trajectory["trend"],("percentiles",),deflatelevel=5)
    defVar(fh,"GMSL_accel",GMSL_Trajectory["accel"],("percentiles",),deflatelevel=5)

    # Write scenarios
    for scenario ∈ settings["NCA5_scenarios"]
        defVar(fh,"GMSL_"*scenario,Float32,("years","percentiles"),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,GMSL_NCA5["years"]),[1.0f0:3.0f0...]),GMSL_NCA5[scenario],extrapolation_bc=NaN32)(settings["years"],[1.0f0:3.0f0...])
    end
    close(fh)
    return nothing
end
end



