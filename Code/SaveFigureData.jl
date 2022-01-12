# ---------------------------------------------------
# Save the data for the figures in Chapter 2
# The figures themselves are created 
# from GMT scripts in the ./GMT/ directory
# ---------------------------------------------------
module SaveFigureData

using NetCDF
using DelimitedFiles
using Statistics
function RunSaveFigureData(settings)
    println("Saving figure data...")
    fig_1_gmsl_usa(settings)
    fig_2_regional(settings)
    fig_4_gmsl(settings)
    fig_5_divergence(settings)
    return
end

function fig_1_gmsl_usa(settings)
    # Comparison of observed global-mean (Panel a) and Contiguous US (Panel b), Observed sea level, trajectory and scenarios
    # Panel a: GMSL
    GMSL_20c = ncread(settings["fn_proj_glb"],"GMSL_20c")
    GMSL_Altimetry = ncread(settings["fn_proj_glb"],"GMSL_Altimetry")
    GMSL_Trajectory = ncread(settings["fn_proj_glb"],"GMSL_Trajectory")
    save_ts_unc_gmt(settings["years"],GMSL_20c./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_20c.txt")
    save_ts_gmt(settings["years"],GMSL_Altimetry./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_Altimetry.txt")
    save_ts_unc_gmt(settings["years"],GMSL_Trajectory./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_Trajectory.txt")
    for scenario ∈ settings["NCA5_scenarios"]
        GMSL_scn = ncread(settings["fn_proj_glb"],"GMSL_"*scenario)
        save_ts_unc_gmt(settings["years"],GMSL_scn./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_"*scenario*".txt")
    end

    # Panel b: Contiguous US
    USA_rsl = ncread(settings["fn_proj_reg"],"rsl_obs")[1,:]
    USA_alt = ncread(settings["fn_proj_reg"],"alt_obs")[1,:]
    USA_trajectory = ncread(settings["fn_proj_reg"],"rsl_trajectory")[1,:,:]
    save_ts_gmt(settings["years"],USA_rsl./1000,settings["dir_fig_1_gmsl_usa"]*"USA_20c.txt")
    save_ts_gmt(settings["years"],USA_alt./1000,settings["dir_fig_1_gmsl_usa"]*"USA_Altimetry.txt")
    save_ts_unc_gmt(settings["years"],USA_trajectory./1000,settings["dir_fig_1_gmsl_usa"]*"USA_Trajectory.txt")

    for scenario ∈ settings["NCA5_scenarios"]
        USA_scn = ncread(settings["fn_proj_reg"],"rsl_total_"*scenario)[1,:,:]
        save_ts_unc_gmt(settings["years"],USA_scn./1000,settings["dir_fig_1_gmsl_usa"]*"USA_"*scenario*".txt")
    end
    return nothing
end

function fig_2_regional(settings)
    # Plot observations, trajectory and scenarios for each individual region
    years = ncread(settings["fn_proj_glb"],"years")
    reg_rsl = ncread(settings["fn_proj_reg"],"rsl_obs")
    reg_alt = ncread(settings["fn_proj_reg"],"alt_obs")
    reg_trajectory = ncread(settings["fn_proj_reg"],"rsl_trajectory")
    for (region_idx,region) ∈ enumerate(settings["regions"])
        save_ts_gmt(years,reg_rsl[region_idx,:]./1000,settings["dir_fig_2_regional"]*region*"_20c.txt")
        save_ts_gmt(years,reg_alt[region_idx,:]./1000,settings["dir_fig_2_regional"]*region*"_Altimetry.txt")
        save_ts_unc_gmt(years,reg_trajectory[region_idx,:,:]./1000,settings["dir_fig_2_regional"]*region*"_Trajectory.txt")
        for scenario ∈ settings["NCA5_scenarios"]
            reg_scn = ncread(settings["fn_proj_reg"],"rsl_total_"*scenario)[region_idx,:,:]
            save_ts_unc_gmt(years,reg_scn./1000,settings["dir_fig_2_regional"]*region*"_"*scenario*".txt")
        end
    end
    return nothing
end

function fig_4_gmsl(settings)
    # Observed global sea level change and its components (Panel a) and GMSL versus contiguous US sea level (Panel b)
    t = [1900:2018...]
    glac = zeros(length(t),3)
    GrIS = zeros(length(t),3)
    AIS = zeros(length(t),3)
    tws = zeros(length(t),3)
    steric = zeros(length(t),3)
    budget = zeros(length(t),3)
    GMSL = zeros(length(t),3)
    read_gmsl_var!("glac",glac,settings)
    read_gmsl_var!("GrIS",GrIS,settings)
    read_gmsl_var!("AIS",AIS,settings)
    read_gmsl_var!("tws",tws,settings)
    read_gmsl_var!("sum_of_contrib_processes",budget,settings)
    read_gmsl_var!("sum_of_contrib_processes",budget,settings)
    read_gmsl_var2!("global_average_sea_level_change",GMSL,settings)
    read_gmsl_var2!("global_average_thermosteric_sea_level_change",steric,settings)
    USA_20c = ncread(settings["fn_proj_reg"],"rsl_obs")[1,:]

    save_ts_unc_gmt(t,0.001.*(glac.-glac[101,2]),settings["dir_fig_4_gmsl"]*"glac.txt")
    save_ts_unc_gmt(t,0.001.*(GrIS.-GrIS[101,2]),settings["dir_fig_4_gmsl"]*"GrIS.txt")
    save_ts_unc_gmt(t,0.001.*(AIS.-AIS[101,2]),settings["dir_fig_4_gmsl"]*"AIS.txt")
    save_ts_unc_gmt(t,0.001.*(tws.-tws[101,2]),settings["dir_fig_4_gmsl"]*"tws.txt")
    save_ts_unc_gmt(t,0.001.*(steric.-steric[101,2]),settings["dir_fig_4_gmsl"]*"steric.txt")
    save_ts_unc_gmt(t,0.001.*(budget.-budget[101,2]),settings["dir_fig_4_gmsl"]*"budget.txt")
    save_ts_unc_gmt(t,0.001.*(GMSL.-GMSL[101,2]) ,settings["dir_fig_4_gmsl"]*"GMSL.txt")
    save_ts_gmt(settings["years"],0.001.*(USA_20c.-USA_20c[101]),settings["dir_fig_4_gmsl"]*"USA_MSL.txt")
    return nothing
end

function fig_5_divergence(settings)
    # Plot where each scenario diverges from the trajectory
    # This script uses 2σ instead of 1σ in the rest of the figures

    # Observed and trajectory time series
    GMSL_20c = ncread(settings["fn_proj_glb"],"GMSL_20c")
    GMSL_Trajectory = ncread(settings["fn_proj_glb"],"GMSL_Trajectory")
    save_ts_unc_gmt(settings["years"],0.001.*(GMSL_20c),settings["dir_fig_5_divergence"]*"GMSL_20c.txt")
    save_ts_unc_gmt(settings["years"],0.001.*(GMSL_Trajectory),settings["dir_fig_5_divergence"]*"GMSL_Trajectory.txt")

    # Panel a: relative to trajectory
    for scenario ∈ settings["NCA5_scenarios"]
        # For each scenario:
        # Read and save time series
        # Determine where this scenario starts diverging
        # from trajectory (Both above and below the trajectory)
        GMSL_scn = ncread(settings["fn_proj_glb"],"GMSL_"*scenario)
        save_ts_unc_gmt(settings["years"],GMSL_scn./1000,settings["dir_fig_5_divergence"]*"GMSL_"*scenario*".txt")
        diff_mean = @. (GMSL_scn[:,2] - GMSL_Trajectory[:,2])
        diff_err = @. 2 * sqrt((0.5*(GMSL_scn[:,3] - GMSL_scn[:,1]))^2 + (0.5*(GMSL_Trajectory[:,3] - GMSL_Trajectory[:,1]))^2)
        scn_above_traj = findfirst(@. (diff_mean - diff_err) >= 0)
        scn_below_traj = findfirst(@. (diff_mean + diff_err) <= 0)
        if !isnothing(scn_above_traj)
            p_plot = [settings["years"][scn_above_traj] 0.001*GMSL_scn[scn_above_traj,2] 1]
            l_plot = [settings["years"][scn_above_traj] 0.001*GMSL_scn[scn_above_traj,2] ; settings["years"][scn_above_traj] -0.05]
            writedlm(settings["dir_fig_5_divergence"]*"div_p_"*scenario*".txt",p_plot,";")
            writedlm(settings["dir_fig_5_divergence"]*"div_l_"*scenario*".txt",l_plot,";")
        elseif !isnothing(scn_below_traj)
            p_plot = [settings["years"][scn_below_traj] 0.001*GMSL_scn[scn_below_traj,2] 1]
            l_plot = [settings["years"][scn_below_traj] 0.001*GMSL_scn[scn_below_traj,2] ; settings["years"][scn_below_traj] -0.05 ]
            writedlm(settings["dir_fig_5_divergence"]*"div_p_"*scenario*".txt",p_plot,";")
            writedlm(settings["dir_fig_5_divergence"]*"div_l_"*scenario*".txt",l_plot,";")
        end
    end

    # Panel b: relative to Intermediate
    GMSL_int = ncread(settings["fn_proj_glb"],"GMSL_Int")
    for scenario ∈ settings["NCA5_scenarios"]
        # For each scenario:
        # Read and save time series
        # Determine where this scenario starts diverging
        # from trajectory (Both above and below the trajectory)
        GMSL_scn = ncread(settings["fn_proj_glb"],"GMSL_"*scenario)
        diff_mean = @. (GMSL_scn[:,2] - GMSL_int[:,2])
        diff_err = @. 2 * sqrt((0.5*(GMSL_scn[:,3] - GMSL_scn[:,1]))^2 + (0.5*(GMSL_int[:,3] - GMSL_int[:,1]))^2)
        scn_above_traj = findfirst(@. (diff_mean - diff_err) >= 0)
        scn_below_traj = findfirst(@. (diff_mean + diff_err) <= 0)
        if !isnothing(scn_above_traj)
            p_plot = [settings["years"][scn_above_traj] 0.001*GMSL_scn[scn_above_traj,2] 1]
            l_plot = [settings["years"][scn_above_traj] 0.001*GMSL_scn[scn_above_traj,2] ; settings["years"][scn_above_traj] -0.05]
            writedlm(settings["dir_fig_5_divergence"]*"div_p_Int_"*scenario*".txt",p_plot,";")
            writedlm(settings["dir_fig_5_divergence"]*"div_l_Int_"*scenario*".txt",l_plot,";")
        elseif !isnothing(scn_below_traj)
            p_plot = [settings["years"][scn_below_traj] 0.001*GMSL_scn[scn_below_traj,2] 1]
            l_plot = [settings["years"][scn_below_traj] 0.001*GMSL_scn[scn_below_traj,2] ; settings["years"][scn_below_traj] -0.05 ]
            writedlm(settings["dir_fig_5_divergence"]*"div_p_Int_"*scenario*".txt",p_plot,";")
            writedlm(settings["dir_fig_5_divergence"]*"div_l_Int_"*scenario*".txt",l_plot,";")
        end
    end
    return nothing
end

# -----------------
# General functions
# -----------------
function read_gmsl_var!(varname,store,settings)
    store[:,1] = ncread(settings["fn_gmsl_20c_mean"],varname*"_lower")
    store[:,2] = ncread(settings["fn_gmsl_20c_mean"],varname*"_mean")
    store[:,3] = ncread(settings["fn_gmsl_20c_mean"],varname*"_upper")
    return nothing
end

function read_gmsl_var2!(varname,store,settings)
    store[:,1] = ncread(settings["fn_gmsl_20c_mean"],varname*"_lower")
    store[:,2] = ncread(settings["fn_gmsl_20c_mean"],varname)
    store[:,3] = ncread(settings["fn_gmsl_20c_mean"],varname*"_upper")
    return nothing
end

function save_ts_unc_gmt(years,ts,fn)
    acc_idx = isfinite.(ts[:,2])
    sdata = [years[acc_idx] ts[acc_idx,2] ts[acc_idx,1] ts[acc_idx,3]]
    writedlm(fn,sdata,";")
    return nothing
end

function save_ts_gmt(years,ts,fn)
    acc_idx = isfinite.(ts)
    sdata = [years[acc_idx] ts[acc_idx]]
    writedlm(fn,sdata,";")
    return nothing
end

end
