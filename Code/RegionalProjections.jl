module RegionalProjections
# --------------------------------------------
# Process regiona; projections and observations
# - Read TG estimates
# - Estimate trend + acceleration using Hector
# - Extrapolate data from tide gauges
# - Read NCA5 projections
# --------------------------------------------
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

function RunRegionalProjections(settings)
    println("Regional trajectories and scenarios...")
    region_data = ReadVirstats(settings)
    ComputeRegionalTrajectory!(region_data,settings)
    NCA_projections = ReadNCA5Regional(region_data,settings)
    SaveData(region_data,NCA_projections,settings)
    println("Regional trajectories and scenarios done")
    return nothing
end

function ReadVirstats(settings)
    println("  Reading regional TG curves...")
    region_data = Dict()
    for region ∈ settings["regions"]
        region_data[region]= Dict()
        region_data[region]["η_tg"] = ncread(settings["fn_regional_obs"],region*"_tg")
    end
    return region_data
end

function ComputeRegionalTrajectory!(region_data,settings)
    println("  Computing trajectories...")
    years_estimate = intersect(settings["years_trajectory"],settings["years_tg"])
    y_acc = findall(in(years_estimate),settings["years"])

    amat = ones(length(y_acc),3)
    amat[:,2] = settings["years"][y_acc] .- mean(settings["years"][y_acc])
    amat[:,3] = (settings["years"][y_acc] .- mean(settings["years"][y_acc])).^2

    amat_extend = ones(length(settings["years_trajectory"]),3)
    amat_extend[:,2] .= settings["years_trajectory"] .- mean(settings["years"][y_acc])
    amat_extend[:,3] = (settings["years_trajectory"] .- mean(settings["years"][y_acc])).^2
    for region ∈ settings["regions"]
        trend_file = Hector.EstTrend(settings["years"][y_acc],region_data[region]["η_tg"][y_acc];accel=true,model="Powerlaw",SA=false,SSA=false,monthly=false)
        μ_sol = [trend_file["bias"],trend_file["trend"],trend_file["accel"]]
        σ_sol= [trend_file["bias_sigma"],trend_file["trend_sigma"],trend_file["accel_sigma"]]
        sol_arr = reshape(μ_sol,(1,3)) .+ randn(5000,3) .* reshape(σ_sol,(1,3))
        trajectory_arr = zeros(Float32,length(settings["years_trajectory"]),5000)
        [trajectory_arr[:,i] = amat_extend * @views sol_arr[i,:] for i ∈ 1:5000]
        trajectory_stats = zeros(Float32,length(settings["years_trajectory"]),3)
        [trajectory_stats[t,:] = quantile((@views trajectory_arr[t,:]),[0.17,0.50,0.83]) for t ∈ 1:length(settings["years_trajectory"])]
        region_data[region]["η_trajectory"] = trajectory_stats
    end
    return nothing
end

function ReadNCA5Regional(region_data,settings)
    println("  Reading NCA5 scenarios...")
    # Read mask
    mask = Masks.ReadMask(settings) 
    area = ComputeGridArea(mask["ϕ"],mask["θ"])
    ϕ,θ,years,percentiles,NCA_grid = ReadNCA5(settings)

    NCA_projections = Dict()
    NCA_projections["years"] = years
    [NCA_projections[region] = Dict() for region ∈ settings["regions"]]
    trajectory_idx = settings["years_trajectory"] .== years[1]
    for scenario ∈ settings["NCA5_scenarios"]
        NCA_int = LinearInterpolation((ϕ, θ,years,percentiles),NCA_grid[scenario],extrapolation_bc = Line())(mask["ϕ"],mask["θ"],years,percentiles);
        for region ∈ settings["regions"]
            NCA_projections[region][scenario] = (sum((@. NCA_int * area*mask[region]),dims=(1,2)) / sum(area .* mask[region]))[1,1,:,:]
            NCA_projections[region][scenario] = NCA_projections[region][scenario] .- NCA_projections[region][scenario][1,2] .+ region_data[region]["η_trajectory"][trajectory_idx,2]
        end
    end
    return NCA_projections
end

function ReadNCA5(settings)
    # Read projections
    ϕ = convert.(Float32,ncread(settings["fn_grid_NCA5"],"lon"))
    θ = convert.(Float32,ncread(settings["fn_grid_NCA5"],"lat"))
    years = convert.(Float32,ncread(settings["fn_grid_NCA5"],"time"))[1:9]
    percentiles = convert.(Float32,ncread(settings["fn_grid_NCA5"],"percentiles"))
    NCA_grid = Dict()
    [NCA_grid[scn] = 10 .* convert.(Float32,ncread(settings["fn_grid_NCA5"],scn,start=[1,1,1,1],count=[-1,-1,9,-1])) for scn ∈ settings["NCA5_scenarios"]]
    return ϕ,θ,years,percentiles,NCA_grid
end

function ComputeGridArea(lon,lat)
    gridsize = abs(lat[2]-lat[1])
    area = @. deg2rad(gridsize) * (sind(lat+gridsize/2)-sind(lat-gridsize/2)) * 6371000^2
    return repeat(area',size(lon,1))
end


function SaveData(region_data,NCA_projections,settings)
    println("  Saving...")
    fh = Dataset(settings["fn_proj_reg"],"c")
    defDim(fh,"years", length(settings["years"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"region", length(settings["regions"]))

    defVar(fh,"region",settings["regions"],("region",),deflatelevel=5)
    defVar(fh,"years",settings["years"],("years",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)

    # Write trajectory
    trajectory = zeros(Float32,length(settings["regions"]),length(settings["years"]),3) .* NaN32
    [trajectory[region_idx,:,:] = (LinearInterpolation((settings["years_trajectory"],[1.0f0:3.0f0...]),region_data[region]["η_trajectory"],extrapolation_bc=NaN32)[settings["years"],[1.0f0:3.0f0...]]) for (region_idx,region) ∈ enumerate(settings["regions"])]
    defVar(fh,"MSL_Trajectory",trajectory,("region","years","percentiles"),deflatelevel=5)

    # Write scenarios
    scen_array = zeros(Float32,length(settings["regions"]),length(settings["years"]),3) .* NaN32
    for scenario ∈ settings["NCA5_scenarios"]
        for (region_idx,region) ∈ enumerate(settings["regions"])
            scen_array[region_idx,:,:] = (LinearInterpolation((NCA_projections["years"],[1.0f0:3.0f0...]),NCA_projections[region][scenario],extrapolation_bc=NaN32)[settings["years"],[1.0f0:3.0f0...]])
        end
        defVar(fh,"MSL_"*scenario,scen_array,("region","years","percentiles",),deflatelevel=5)
    end
    close(fh)
end

end