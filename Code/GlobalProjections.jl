# -------------------------------------------
# Process global projections and observations
# - Read altimetry observation (GSFC/MEASURES)
# - Read TG reconstruction ensemble (F2020)
# - Estimate trend + acceleration using Hector
# - Extrapolate data from tide gauges
# - Read AR6 projections
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
    println("Global-mean observations and scenarios...")
    GMSL_20c = read_GMSL_20c(settings)
    GMSL_Trajectory = compute_GMSL_trajectory(GMSL_20c,settings)
    GMSL_Altimetry = read_GMSL_Altimetry(settings)
    GMSL_NCA5 = read_NCA5(GMSL_Trajectory,settings)
    save_data(GMSL_20c,GMSL_Trajectory,GMSL_Altimetry,GMSL_NCA5,settings)
    println("Global-mean observations and scenarios done")
    return nothing
end

function read_GMSL_20c(settings)
    println("  Reading tide-gauge GMSL...")
    # Read F2020 GMSL
    fh = Dataset(settings["fn_GMSL_20c"],"r")
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
    GMSL_20c["η_mean"] = η_mean
    GMSL_20c["Λ"] = Λ
    GMSL_20c["η"] = η
    return GMSL_20c
end

function compute_GMSL_trajectory(GMSL_20c,settings)
    println("  Computing GMSL trajectories...")
    t_acc = findall(in(settings["years_trajectory"]),GMSL_20c["years"])
    # Compute and extrapolate trend and acceleration
    amat = ones(length(t_acc),3)
    amat[:,2] = GMSL_20c["years"][t_acc] .- mean(GMSL_20c["years"][t_acc])
    @. amat[:,3] = (GMSL_20c["years"][t_acc] - $mean(GMSL_20c["years"][t_acc]))^2
    amat_tr = transpose(amat)
    amat_sq = inv(amat_tr*amat)
    sol_arr = zeros(length(GMSL_20c["Λ"]),3)
    amat_expand = ones(length(settings["years_trajectory"]),3)
    amat_expand[:,2] = settings["years_trajectory"] .- mean(GMSL_20c["years"][t_acc])
    @. amat_expand[:,3] = (settings["years_trajectory"] - $mean(GMSL_20c["years"][t_acc]))^2
    η_projected = zeros(length(settings["years_trajectory"]),length(GMSL_20c["Λ"]))
    noise = Hector.GenerateNoise(GMSL_20c["years"][t_acc],GMSL_20c["η_mean"][t_acc,2],length(GMSL_20c["Λ"]),length(GMSL_20c["years"]);accel=true,model="Powerlaw");
    GMSL_20c["η"] .+= noise
    for λ ∈ 1:length(GMSL_20c["Λ"])
        sol_arr[λ,:] = amat_sq * (amat_tr*GMSL_20c["η"][t_acc,λ])
        η_projected[:,λ] = amat_expand * sol_arr[λ,:]
    end
    # Perturb estimated trend
    # Projection mean and std
    projection_mean = zeros(length(settings["years_trajectory"]),3)
    sidx = zeros(Int,length(GMSL_20c["Λ"]))
    𝚲    = zeros(Float64,length(GMSL_20c["Λ"]))
    Η    = zeros(Float64,length(GMSL_20c["Λ"]))
    for ts ∈ 1:length(settings["years_trajectory"])
        sortperm!(sidx, η_projected[ts,:])
        cumsum!(𝚲,GMSL_20c["Λ"][sidx])
        @. Η = @views η_projected[ts,sidx];
        projection_mean[ts,1] = Η[findfirst(>(0.05),𝚲)]
        projection_mean[ts,2] = Η[findfirst(>(0.5),𝚲)]
        projection_mean[ts,3] = Η[findfirst(>(0.95),𝚲)]
    end
    GMSL_Trajectory = Dict()
    GMSL_Trajectory["years"] = settings["years_trajectory"]
    GMSL_Trajectory["η_mean"] = projection_mean
    GMSL_Trajectory["Λ"] = GMSL_20c["Λ"]
    GMSL_Trajectory["η"] = η_projected
    return GMSL_Trajectory
end

function read_GMSL_Altimetry(settings)
    println("  Reading altimetry GMSL...")
    # Read Altimetry GMSL from Goddard and apply corrections for
    # internal variability due to ENSO/PDO and take annual-means
    GMSL_Altimetry = Dict()
    GMSL_Altimetry["years"] = [1993:2020...]
      # Read data
    GSFC_raw = readdlm(settings["fn_GMSL_GSFC"],skipstart=48)[:,[3,12]]
    ENSO_raw = matread(settings["fn_GMSL_ENSO"])
    # To monthly-mean GMSL 
    tbounds = zeros(length(ENSO_raw["time"]),2)
    tbounds[:,1] = [1993.0:1/12:2021-1/12...]
    tbounds[:,2] = tbounds[:,1] .+ 1/12
    GSFC_monthly = zeros(Float32,length(ENSO_raw["time"]))
    [GSFC_monthly[t] = mean(GSFC_raw[(GSFC_raw[:,1] .>= tbounds[t,1]) .& (GSFC_raw[:,1] .<= tbounds[t,2]),2]) - ENSO_raw["gmsl"][1,t] for t ∈ 1:length(ENSO_raw["time"])]
     GMSL_Altimetry["η_mean"] =mean(reshape(GSFC_monthly,(12,:)),dims=1)[1,:] # Monthly-mean to annual-mean
    # Remove baseline
    idx_baseline = findall(in(settings["years_baseline"]),GMSL_Altimetry["years"])
    GMSL_Altimetry["η_mean"] .-= mean(GMSL_Altimetry["η_mean"][idx_baseline])
    return GMSL_Altimetry
end

function read_NCA5(GMSL_Trajectory,settings)
    println("  Reading NCA5 scenarios...")
    # Read the GMSL scenarios
    GMSL_NCA5 = Dict()
    GMSL_NCA5["years"] = ncread(settings["fn_GMSL_NCA5"],"time",start=[1],count=[9])
    [GMSL_NCA5[scenario] = 10.0f0 .* ncread(settings["fn_GMSL_NCA5"],scenario,start=[1,1],count=[9,-1]) for scenario ∈ settings["NCA5_scenarios"]]
    # Tie trajectory to scenario in 2020
    traj_idx = findall(GMSL_Trajectory["years"] .== 2020)
    [GMSL_NCA5[scenario] = @. GMSL_NCA5[scenario] - GMSL_NCA5[scenario][1,2] + GMSL_Trajectory["η_mean"][traj_idx,2] for scenario ∈ settings["NCA5_scenarios"]]
    return GMSL_NCA5
end

function save_data(GMSL_20c,GMSL_Trajectory,GMSL_Altimetry,GMSL_NCA5,settings)
    println("  Saving data...")
    fh = Dataset(settings["fn_proj_glb"],"c")
    defDim(fh,"years", length(settings["years"]))
    defDim(fh,"percentiles",3)
    defVar(fh,"years",settings["years"],("years",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)
    # Write observations and trajectory
    defVar(fh,"GMSL_20c",Float32,("years","percentiles"),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,GMSL_20c["years"]),[1.0f0:3.0f0...]),GMSL_20c["η_mean"],extrapolation_bc=NaN32)(settings["years"],[1.0f0:3.0f0...])
    defVar(fh,"GMSL_Altimetry",Float32,("years",),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,GMSL_Altimetry["years"])),GMSL_Altimetry["η_mean"],extrapolation_bc=NaN32)(settings["years"])
    defVar(fh,"GMSL_Trajectory",Float32,("years","percentiles"),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,GMSL_Trajectory["years"]),[1.0f0:3.0f0...]),GMSL_Trajectory["η_mean"],extrapolation_bc=NaN32)(settings["years"],[1.0f0:3.0f0...])
    # Write scenarios
    for scenario ∈ settings["NCA5_scenarios"]
        defVar(fh,"GMSL_"*scenario,Float32,("years","percentiles"),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,GMSL_NCA5["years"]),[1.0f0:3.0f0...]),GMSL_NCA5[scenario],extrapolation_bc=NaN32)(settings["years"],[1.0f0:3.0f0...])
    end
    close(fh)
    return nothing
end
end



