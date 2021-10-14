using Plots
using NetCDF
using DelimitedFiles
# Make simple plots and save for GMT
function simple_plots(settings)

    return
end

function figure_1(settings)
    years = ncread(settings["fn_proj_glb"],"years")
    # 1a GMSL
    GMSL_20c = ncread(settings["fn_proj_glb"],"GMSL_20c")
    GMSL_Altimetry = ncread(settings["fn_proj_glb"],"GMSL_Altimetry")
    GMSL_Trajectory = ncread(settings["fn_proj_glb"],"GMSL_Trajectory")
    # Save time series
    settings["dir_fig_2_regional"]
    save_ts_unc_gmt(years,GMSL_20c./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_20c.txt")
    save_ts_gmt(years,GMSL_Altimetry./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_Altimetry.txt")
    save_ts_unc_gmt(years,GMSL_Trajectory./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_Trajectory.txt")
    plot(years,GMSL_20c,color=:black)
    plot!(years,GMSL_Trajectory,color=:red,legend=false)
    for scenario ∈ settings["NCA5_scenarios"]
        GMSL_scn = ncread(settings["fn_proj_glb"],"GMSL_"*scenario)
        plot!(years,GMSL_scn[:,2],color=:green,legend=false)
        save_ts_unc_gmt(years,GMSL_scn./1000,settings["dir_fig_1_gmsl_usa"]*"GMSL_"*scenario*".txt")
    end
    plot!(years,GMSL_Altimetry,color=:blue)
    plot!()

    # 1b USA
    USA_20c = ncread(settings["fn_regional_obs"],"USA_tg")
    USA_rsl_alt = ncread(settings["fn_regional_obs"],"USA_rsl_alt")
    USA_trajectory = ncread(settings["fn_proj_reg"],"MSL_Trajectory")[1,:,:]
    save_ts_gmt(years,USA_20c./1000,settings["dir_fig_1_gmsl_usa"]*"USA_20c.txt")
    save_ts_gmt(years,USA_rsl_alt./1000,settings["dir_fig_1_gmsl_usa"]*"USA_Altimetry.txt")
    save_ts_unc_gmt(years,USA_trajectory./1000,settings["dir_fig_1_gmsl_usa"]*"USA_Trajectory.txt")

    plot(years,USA_20c,color=:black)
    plot!(years,USA_trajectory,color=:red,legend=false)
    for scenario ∈ settings["NCA5_scenarios"]
        USA_scn = ncread(settings["fn_proj_reg"],"MSL_"*scenario)[1,:,:]
        plot!(years,USA_scn[:,2],color=:green,legend=false)
        save_ts_unc_gmt(years,USA_scn./1000,settings["dir_fig_1_gmsl_usa"]*"USA_"*scenario*".txt")
    end
    plot!(years,USA_rsl_alt,color=:blue,legend=false)
end

function figure_2(settings)
    plot_array = Any[]
    years = ncread(settings["fn_proj_glb"],"years")

    for (region_idx,region) ∈ enumerate(settings["regions"][2:end])
        reg_20c = ncread(settings["fn_regional_obs"],region*"_tg")
        reg_rsl_alt = ncread(settings["fn_regional_obs"],region*"_rsl_alt")
        reg_trajectory = ncread(settings["fn_proj_reg"],"MSL_Trajectory")[region_idx+1,:,:]
        save_ts_gmt(years,reg_20c./1000,settings["dir_fig_2_regional"]*region*"_20c.txt")
        save_ts_gmt(years,reg_rsl_alt./1000,settings["dir_fig_2_regional"]*region*"_Altimetry.txt")
        save_ts_unc_gmt(years,reg_trajectory./1000,settings["dir_fig_2_regional"]*region*"_Trajectory.txt")
        p = plot(years,reg_20c,color=:black)
        plot!(years,reg_trajectory,color=:red,legend=false)
        for scenario ∈ settings["NCA5_scenarios"]
            reg_scn = ncread(settings["fn_proj_reg"],"MSL_"*scenario)[region_idx+1,:,:]
            save_ts_unc_gmt(years,reg_scn./1000,settings["dir_fig_2_regional"]*region*"_"*scenario*".txt")
            plot!(years,reg_scn[:,2],color=:green,legend=false)
        end
        plot!(years,reg_rsl_alt,color=:blue,legend=false)
        title!(region)
        push!(plot_array,p)
    end
    plot(plot_array..., layout = (5,2),size=(1000,1200))
    savefig("extrap_regr_alt.png")
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
