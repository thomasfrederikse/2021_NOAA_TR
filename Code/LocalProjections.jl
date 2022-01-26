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
    fn = settings["dir_NCA5"]*"NCA5_Low_grid.nc"
    lon_NCA5 = ncread(fn,"lon")
    lat_NCA5 = ncread(fn,"lat")
    years_NCA5 = convert.(Float32,ncread(fn,"years",start=[1],count=[-1]))
    pct_NCA5 = convert.(Float32,ncread(fn,"percentiles"))

    # Find nearest grid cell
    NCA5_loc = zeros(Int,length(local_obs["name"]),2)
    for tg in 1:length(local_obs["name"])
        NCA5_loc[tg,1] = argmin(@. abs(local_obs["coords"][tg,1]-lon_NCA5))
        NCA5_loc[tg,2] = argmin(@. abs(local_obs["coords"][tg,2]-lat_NCA5))
    end

    # Create data sctructure
    NCA5_local = Array{Dict}(undef,length(local_obs["name"]))
    for tg in 1:length(local_obs["name"])
        NCA5_local[tg] = Dict()
        for scn in settings["NCA5_scenarios"]
            NCA5_local[tg][scn] = Dict()
        end
    end

    # Read data
    for scenario in settings["NCA5_scenarios"]
        println("   Scenario "*scenario*"...")
        fn = settings["dir_NCA5"]*"NCA5_"*scenario*"_grid.nc"
        for prc in settings["processes"]
            NCA5_prc = convert.(Float32,ncread(fn,prc,start=[1,1,1,1],count=[-1,-1,-1,-1]));
            for tg in 1:length(local_obs["name"])
                NCA5_local[tg][scenario][prc] = NCA5_prc[NCA5_loc[tg,1],NCA5_loc[tg,2],:,:]
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
    defDim(fh,"years", length(settings["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"tg", length(local_obs["name"]))

    defVar(fh,"tg",local_obs["name"],("tg",),deflatelevel=5)
    defVar(fh,"lon",local_obs["coords"][:,1],("tg",),deflatelevel=5)
    defVar(fh,"lat",local_obs["coords"][:,2],("tg",),deflatelevel=5)
    defVar(fh,"years",settings["years"],("years",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)

    # Write trajectory
    defVar(fh,"rsl_obs",rsl_obs,("tg","years"),deflatelevel=5)
    defVar(fh,"number_of_observation_years_for_trajectory_estimation",local_obs["obs_length"],("tg",),deflatelevel=5)
    defVar(fh,"number_of_observation_years_for_trajectory_estimation_geq_30",local_obs["obs_lt_30"],("tg",),deflatelevel=5)
    defVar(fh,"rsl_trajectory",local_obs["rsl_trajectory"],("tg","years","percentiles"),deflatelevel=5)
    defVar(fh,"rsl_trend",local_obs["rsl_trend"],("tg","percentiles"),deflatelevel=5)
    defVar(fh,"rsl_accel",local_obs["rsl_accel"],("tg","percentiles"),deflatelevel=5)
    # Write scenarios
    for scn ∈ settings["NCA5_scenarios"]
        for prc in settings["processes"]
            defVar(fh,"rsl_"*prc*"_"*scn,scn_proj[scn][prc],("tg","years","percentiles",),deflatelevel=5)
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