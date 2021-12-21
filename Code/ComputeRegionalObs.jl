# --------------------------------------------------------------------
# Run a simple virtual-station average for US tide-gauge observations
# Read altimetry observations
# Provide annual-mean reconstruction
# --------------------------------------------------------------------
module ComputeRegionalObs
using NCDatasets
using NetCDF
using XLSX
using Dates
using Statistics
using Interpolations
using DelimitedFiles
dir_code = homedir()*"/Projects/2021_NOAA_TR/Code/"
include(dir_code*"Masks.jl")

function RunComputeRegionalObs(settings)
    println("\nRegional observations...")
    tg_monthly = ReadTGData(settings)
    tg_annual = ConvertMonthlyToAnnual(tg_monthly,settings)
    region_data = MapToRegions(tg_annual,settings)
    MergeVirstat!(region_data,settings)
    RemoveClimVar!(region_data,settings)
    ReadAltimetry!(region_data,settings)
    SaveData(region_data,settings)
    SaveGridGMT(region_data, settings)
    println("Regional observations done\n")
    return nothing
end

function ReadTGData(settings)
    # --------------------------------------------------
    # Read the monthly tide-gauge data provided by Billy
    # Returns monthly-mean sea level in mm
    # --------------------------------------------------
    println("  Reading tide-gauge data...")
    tg_monthly = Dict()
    data_raw = XLSX.readxlsx(settings["fn_tg_data"])["monthly MSL"][:]
    tg_monthly["station_names"] = Vector{String}(data_raw[1,2:122])
    tg_monthly["station_coords"] = reverse!(Array{Float32}(data_raw[2:3,2:122]'),dims=2)
    @. tg_monthly["station_coords"][tg_monthly["station_coords"][:,1]<0,1] += 360  # Longitudes should cover [0-360]
    tval = data_raw[4:end,1]
    tvec = zeros(Int,length(tval),2)
    for t ∈ 1:length(tval)
        tvec[t,1] = Year(tval[t]).value
        tvec[t,2] = Month(tval[t]).value
    end
    tg_monthly["time"] = tval
    tg_monthly["tvec"] = tvec
    tg_monthly["tval"] = tg_monthly["tvec"][:,1] + (tg_monthly["tvec"][:,2]./12 .- 1/24)

    rsl_raw = data_raw[4:end,2:122]
    @. rsl_raw[rsl_raw=="NaN"] = NaN
    tg_monthly["rsl"] = Array{Float32}(rsl_raw')
    tg_monthly["rsl"] .*= 1000.0f0

    # Remove Seward because of jump
    stat_acc = ones(Bool,size(tg_monthly["rsl"],1))
    stat_acc[107] = false
    tg_monthly["station_names"] = tg_monthly["station_names"][stat_acc]
    tg_monthly["rsl"] = tg_monthly["rsl"][stat_acc,:]
    tg_monthly["station_coords"] = tg_monthly["station_coords"][stat_acc,:]
    return tg_monthly
end

function ConvertMonthlyToAnnual(tg_monthly,settings)
    # ------------------------------------------------------------
    # Compute annual-mean sea level from monthly-mean sea level
    # - First remove mean seasonal cycle, which would give trouble
    #   for years with missing months
    # - For years with at least 10 months with data, compute
    #   annual mean. Otherwise, NaN32
    # ------------------------------------------------------------
    println("  Converting tide-gauge data to annual-means...")
    tg_annual = Dict()
    tg_annual["years"] = settings["years_tg"]
    tg_annual["station_names"] = tg_monthly["station_names"]
    tg_annual["station_coords"] = tg_monthly["station_coords"]
    tg_annual["rsl"] = zeros(Float32,size(tg_monthly["rsl"],1),length(tg_annual["years"])) .* NaN32
    for stat ∈ 1:size(tg_monthly["rsl"],1)
        rsl_annual_lcl = zeros(length(tg_annual["years"])) .* NaN32
        rsl_monthly_lcl = zeros(length(tg_monthly["tval"])) .* NaN32
        # Remove seasonal cycle
        acc_idx = @. isfinite(tg_monthly["rsl"][stat,:])
        amat = ones(sum(acc_idx),6)
        amat[:,2] = tg_monthly["tval"][acc_idx] .- mean(tg_monthly["tval"][acc_idx])
        amat[:,3] = @. sin(2*π*tg_monthly["tval"][acc_idx])
        amat[:,4] = @. cos(2*π*tg_monthly["tval"][acc_idx])
        amat[:,5] = @. sin(4*π*tg_monthly["tval"][acc_idx])
        amat[:,6] = @. cos(4*π*tg_monthly["tval"][acc_idx])
        sol = amat\tg_monthly["rsl"][stat,acc_idx]    
        sol[2] = 0
        rsl_monthly_lcl[acc_idx] = tg_monthly["rsl"][stat,acc_idx] - amat*sol
        rsl_monthly_lcl = reshape(rsl_monthly_lcl,(12,length(tg_annual["years"])))
        for yr ∈ 1:length(tg_annual["years"])
            if sum(isnan.(rsl_monthly_lcl[:,yr])) < 3
                rsl_annual_lcl[yr] = mean(filter(isfinite,rsl_monthly_lcl[:,yr]))
            end
        end
        tg_annual["rsl"][stat,:] = rsl_annual_lcl
    end
    return tg_annual
end

function MapToRegions(tg_annual,settings)
    # --------------------------------------------------------------
    # Map the individual tide-gauge records to one of the 10 regions
    # --------------------------------------------------------------
    println("  Mapping stations onto regions...")
    mask = Masks.ReadMask(settings)
    region_num = zeros(Bool,length(settings["regions"]),size(tg_annual["station_coords"],1))
    
    # Alaska stations
    in_alaska = zeros(Bool,size(tg_annual["station_coords"],1))
    for stat ∈ 1:size(tg_annual["station_coords"],1)
        strip(tg_annual["station_names"][stat][end-2:end]) == "AK" ? in_alaska[stat] = true : nothing
    end
    ALN_list = ["Unalaska, AK","Prudhoe Bay, AK","Nome, AK","Adak Island, AK ","Unalaska, AK"]

    for (region_idx,region) ∈ enumerate(settings["regions"])
        if (region!="ALN") & (region!="ALS")
            for stat ∈ 1:size(tg_annual["station_coords"],1)
                ϕ_idx = argmin(@. abs(tg_annual["station_coords"][stat,1] - mask["ϕ"]))
                θ_idx = argmin(@. abs(tg_annual["station_coords"][stat,2] - mask["θ"]))
                if mask[region*"_unc"][ϕ_idx,θ_idx]
                    region_num[region_idx,stat] = true
                end
            end
        elseif region=="ALN"
            for stat ∈ 1:size(tg_annual["station_coords"],1)
                tg_annual["station_names"][stat] in ALN_list ? region_num[region_idx,stat] = true : nothing
            end
        elseif region=="ALS"
            for stat ∈ 1:size(tg_annual["station_coords"],1)
                if in_alaska[stat] & !region_num[region_idx-1,stat]
                    region_num[region_idx,stat] = true
                end
            end
        end
    end
    region_data = Dict()
    for (region_idx,region) ∈ enumerate(settings["regions"])
        region_data[region] = Dict()
        region_data[region]["station_coords"] = tg_annual["station_coords"][region_num[region_idx,:],:]
        region_data[region]["station_names"] = tg_annual["station_names"][region_num[region_idx,:]]
        region_data[region]["rsl"] = tg_annual["rsl"][region_num[region_idx,:],:]
    end
    region_data["years"] = tg_annual["years"]
    return region_data
end

function MergeVirstat!(region_data,settings)
    # ------------------------------------------------------------------------------------------------
    # Merge the tide gauges in each region to an average curve using the
    # Dangendorf et al. (2017) virtual station method:
    # Iteratively find the two nearest stations, subtract local mean sea level
    # over the years where both stations have data from both stations. Then take
    # the average RSL of both stations. Remove these two stations from the list,
    # but add a new one with the averaged sea-level curve halfway both stations.
    # Repeat this procedure until only one virtual station is left. The resulting
    # curve is the estimate of regional sea level.
    #
    # For more details, see:
    #  Dangendorf, S., Marcos, M., Wöppelmann, G., Conrad, C. P., Frederikse, T., & Riva, R. (2017). 
    #  Reassessment of 20th century global mean sea level rise. Proceedings of the National Academy 
    #  of Sciences, 114(23), 5946–5951. https://doi.org/10.1073/pnas.1616007114
    # ------------------------------------------------------------------------------------------------
    println("  Merge stations into regional means using virtual-station approach...")
    idx_baseline = findall(in(settings["years_baseline"]),settings["years_tg"])
    for region ∈ settings["regions"]
        rsl_reg = copy(region_data[region]["rsl"])
        loc_reg = copy(region_data[region]["station_coords"])
        while size(loc_reg,1) > 1
            # Find two nearest stations
            merge_idx, ϕ_mid, θ_mid = FindStatsToMerge(rsl_reg,loc_reg)
            merge_data = rsl_reg[[merge_idx[1],merge_idx[2]],:]
            # Subtract the mean over the overlap period
            ovl = @. isfinite(merge_data[1,:]) & isfinite(merge_data[2,:])
            merge_data[1,:] .-= mean(merge_data[1,ovl])
            merge_data[2,:] .-= mean(merge_data[2,ovl])
            merged_data = MergeTwoStats(merge_data)
            # Replace stations by virtaul station
            keep_idx = ones(Bool,size(rsl_reg,1))
            keep_idx[[merge_idx[1],merge_idx[2]]] .= false
            rsl_reg = rsl_reg[keep_idx,:]
            loc_reg = loc_reg[keep_idx,:]
            rsl_reg = vcat(rsl_reg, merged_data')
            loc_reg = vcat(loc_reg, [ϕ_mid, θ_mid]')
        end
        region_data[region]["rsl_virstat"] = rsl_reg[:]
        region_data[region]["rsl_virstat"] .-= mean(region_data[region]["rsl_virstat"][idx_baseline])
    end
    return nothing
end


function ReadAltimetry!(region_data,settings)
    # ------------------------------------------------------------------
    # Read satellite altimetry and add GIA and GRD corrections to obtain
    # relative sea level, which is comparable to tide-gauge observations
    # ------------------------------------------------------------------
    println("  Reading altimetry...")
    mask = Masks.ReadMask(settings) 
    area = ComputeGridArea(mask["ϕ"],mask["θ"])
    clim_indices = ReadClimIndices(settings)

    # Read altimetry data
    ssh_grid = convert.(Float32,ncread(settings["fn_altimetry"],"ssh"))
    @. ssh_grid[ssh_grid > 31000] = 0
    @. ssh_grid .*= 0.05
    alt_time = convert.(Float32,ncread(settings["fn_altimetry"],"time"))
    alt_year = [1993:2019...]
    # From high frquency data to annual-means
    ssh_grid_year = zeros(Float32,size(ssh_grid,1),size(ssh_grid,2),length(alt_year));
    for (idx,yr) ∈ enumerate(alt_year)
        acc_idx = floor.(alt_time) .== yr
        ssh_grid_year[:,:,idx] = mean(ssh_grid[:,:,acc_idx],dims=3)
    end

    # Take regional average and remove NAO, PDO and ENSO effects using
    # ordinary least squares
    amat = ones(size(alt_year,1),5);
    amat[:,2] = alt_year .- mean(alt_year);
    amat[:,3] = LinearInterpolation(clim_indices[1]["years"], clim_indices[1]["index"],extrapolation_bc=Line())(alt_year)
    amat[:,4] = LinearInterpolation(clim_indices[2]["years"], clim_indices[2]["index"],extrapolation_bc=Line())(alt_year)
    amat[:,5] = LinearInterpolation(clim_indices[3]["years"], clim_indices[3]["index"],extrapolation_bc=Line())(alt_year)
    for region ∈ settings["regions"]
        ssh_reg_raw = sum( @.(ssh_grid_year * mask[region] * area),dims=(1,2))[:] ./ sum(mask[region] .* area)
        sol = amat\ssh_reg_raw
        sol[2] = 0
        region_data[region]["gsl_alt"] = ssh_reg_raw - amat*sol
    end

    # 3 Add GIA and contemporary GRD
    GIA = ComputeGIA(mask,area,settings) # GIA from Caron et al. (2019)
    GRD = ComputeGRD(mask,area,alt_year,settings) # Contemporary GRD effects from Frederikse e al. (2020)
    for (region_idx,region) ∈ enumerate(settings["regions"])
        region_data[region]["rsl_alt"] = region_data[region]["gsl_alt"] - (GIA[region_idx] .* (alt_year .- mean(alt_year))) - GRD[:,region_idx]
    end

    # Adjust baseline of all observations
    idx_baseline = findall(in(settings["years_baseline"]),alt_year)
    for region ∈ settings["regions"]
        region_data[region]["rsl_alt"] .-= mean(region_data[region]["rsl_alt"][idx_baseline])
        region_data[region]["gsl_alt"] .-= mean(region_data[region]["gsl_alt"][idx_baseline])
    end
    region_data["alt_year"] = alt_year
    return nothing
end

function ComputeGIA(mask,area,settings)
    # --------------------------------------------------------------------------------------
    # Read GIA solid-Earth deformation from:
    # Caron, L., Ivins, E. R., Larour, E., Adhikari, S., Nilsson, J., & Blewitt, G. (2018). 
    # GIA Model Statistics for GRACE Hydrology, Cryosphere, and Ocean Science. 
    # Geophysical Research Letters, 45(5), 2203–2212. https://doi.org/10.1002/2017GL076644
    # --------------------------------------------------------------------------------------
    println("  Reading GIA...")
    GIA = zeros(Float32,length(settings["regions"]))
    GIA_rad = ncread(settings["fn_GIA"],"rad_mean")
    for (region_idx,region) ∈ enumerate(settings["regions"])
        GIA[region_idx] = sum(@.(GIA_rad * mask[region] * area)) ./ sum(mask[region] .* area)
    end
    return GIA
end

function ComputeGRD(mask,area,years,settings)
    # --------------------------------------------------------------------------------------
    # Read GRD solid-Earth deformation from:
    # Frederikse, T., Landerer, F., Caron, L., Adhikari, S., Parkes, D., Humphrey, V. W., 
    # Dangendorf, S., Hogarth, P., Zanna, L., Cheng, L., & Wu, Y.-H. (2020). 
    # The causes of sea-level rise since 1900. 
    # Nature, 584(7821), 393–397. https://doi.org/10.1038/s41586-020-2591-3
    # --------------------------------------------------------------------------------------
    println("  Reading contemporary GRD...")
    GRD = zeros(Float32,length(years),length(settings["regions"]))
    GRD_rad = ncread(settings["fn_GRD"],"rad") .* 0.05
    GRD_years = ncread(settings["fn_GRD"],"time")
    for (region_idx,region) ∈ enumerate(settings["regions"])
        GRD_lcl = sum( @.(GRD_rad * mask[region] * area),dims=(1,2))[:] ./ sum(mask[region] .* area)
        GRD[:,region_idx] = LinearInterpolation(GRD_years, GRD_lcl,extrapolation_bc=Line())(years)
    end
    return GRD
end

function RemoveClimVar!(region_data,settings)
    # -------------------------------------------------------------------
    # For each region remove the effects of climate indices (NAO/PDO/MEI)
    # from observed sea level using OLS. 
    # -------------------------------------------------------------------
    clim_indices = ReadClimIndices(settings)
    for region ∈ settings["regions"]
        region_data[region]["rsl_virstat_novar"] = RemoveVariability(settings["years_tg"],region_data[region]["rsl_virstat"],clim_indices,settings)
    end
    return nothing
end

function RemoveVariability(years,tseries,clim_indices,settings)
    acc_idx = @. isfinite(tseries)
    amat = ones(sum(acc_idx),5);
    amat[:,2] = years[acc_idx].- mean(years[acc_idx]);
    amat[:,3] = LinearInterpolation(clim_indices[1]["years"], clim_indices[1]["index"],extrapolation_bc=Line())(years[acc_idx])
    amat[:,4] = LinearInterpolation(clim_indices[2]["years"], clim_indices[2]["index"],extrapolation_bc=Line())(years[acc_idx])
    amat[:,5] = LinearInterpolation(clim_indices[3]["years"], clim_indices[3]["index"],extrapolation_bc=Line())(years[acc_idx])
    sol = amat\tseries[acc_idx]
    sol[2] = 0
    tseries_var = zeros(length(tseries)) .* NaN32
    tseries_var[acc_idx] = amat*sol
    tseries_novar = tseries - tseries_var
    idx_baseline = findall(in(settings["years_baseline"]),years)
    tseries_var .-= mean(tseries_var[idx_baseline])
    tseries_novar .-= mean(tseries_novar[idx_baseline])
    return tseries_novar
end

function ReadClimIndices(settings)
    # ---------------------------------------------------------------
    # Read climate indices used for removing natural variability
    # All indices come from NOAA Physical Sciences Laboratory (PSL)
    # https://psl.noaa.gov/data/climateindices/
    # and NOAA Climate Prediction Centre (CPC)
    # https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml
    # ---------------------------------------------------------------
    clim_indices = Array{Dict}(undef,3)
    [clim_indices[i] = Dict() for i ∈ 1:3]

    # NAO
    clim_indices[1]["name"] = ["NAO"]
    clim_indices[1]["years"] = [1899:2020...]
    clim_indices[1]["index"] =  mean(reshape(readdlm(settings["fn_NAO"],skipstart=1)[1:122,2:end]'[:],(12,:)),dims=1)[1,:]

    # PDO
    clim_indices[2]["name"] = ["PDO"]
    clim_indices[2]["years"] = [1854:2020...]
    clim_indices[2]["index"] =  mean(reshape(readdlm(settings["fn_PDO"],skipstart=2)[1:167,2:end]'[:],(12,:)),dims=1)[1,:]

    # MEI
    MEI_1 = readdlm(settings["fn_MEI_1"])[1:end,2:13]'[:]
    MEI_1_time = [1871+1/24:1/12:2005+23/24...]
    MEI_2 = readdlm(settings["fn_MEI_2"],skipstart=1)[1:42,2:13]'[:]
    MEI_2_time = [1979+1/24:1/12:2020+23/24...]
    MEI_time = union(MEI_1_time,MEI_2_time)
    MEI = zeros(length(MEI_time))
    MEI[end-503:end] = MEI_2
    idx_1 = findall(in(MEI_2_time),MEI_1_time) 
    idx_2 = findall(in(MEI_1_time),MEI_2_time) 
    idx_3 = findall(!in(MEI_2_time),MEI_1_time) 
    MEI_1 = MEI_1 .- mean(MEI_1[idx_1]) .+ mean(MEI_2[idx_2])
    MEI[1:end-504] = MEI_1[idx_3]
    clim_indices[3]["name"] = ["MEI"]
    clim_indices[3]["years"] = [1871:2020...]
    clim_indices[3]["index"] =  mean(reshape(MEI,(12,:)),dims=1)[1,:]
    return clim_indices
end

function SaveData(region_data,settings)
    println("  Saving...")
    # Save the time series
    fh = Dataset(settings["fn_regional_obs"],"c")
    defDim(fh,"years", length(settings["years"]))
    defVar(fh,"years",convert.(Float32,settings["years"]),("years",),deflatelevel=5)
    for region ∈ settings["regions"]     
        defVar(fh,region*"_tg",Float32,("years",),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,settings["years_tg"])),region_data[region]["rsl_virstat_novar"],extrapolation_bc=NaN32)(settings["years"])
        defVar(fh,region*"_rsl_alt",Float32,("years",),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,region_data["alt_year"] )),region_data[region]["rsl_alt"],extrapolation_bc=NaN32)(settings["years"])
        defVar(fh,region*"_gsl_alt",Float32,("years",),deflatelevel=5)[:] = LinearInterpolation((convert.(Float32,region_data["alt_year"] )),region_data[region]["gsl_alt"],extrapolation_bc=NaN32)(settings["years"])
    end
    close(fh)
    return nothing
end

function SaveGridGMT(region_data, settings)
    for region ∈ settings["regions"]     
        writedlm(settings["dir_fig_3_map"]*region*".txt",region_data[region]["station_coords"],";")
    end
    mask = Masks.ReadMask(settings) 
    mask_num = zeros(Int16,size(mask["EC"]))
    mask_num[mask["EC"]] .= 1
    mask_num[mask["SE"]] .= 2
    mask_num[mask["GCE"]] .= 3
    mask_num[mask["GCW"]] .= 4
    mask_num[mask["SWC"]] .= 5
    mask_num[mask["NWC"]] .= 6
    mask_num[mask["PAC"]] .= 7
    mask_num[mask["CAR"]] .= 8
    mask_num[mask["ALN"]] .= 9
    mask_num[mask["ALS"]] .= 10

    return nothing
end

### General functions
function MergeTwoStats(merge_data)
    merged_data = zeros(Float32,size(merge_data,2))
    for t ∈ 1:length(merged_data)
        acc = isfinite.(merge_data[:,t])
        merged_data[t] = ifelse(sum(acc)>0,mean(merge_data[acc,t]),NaN)            
    end
    return merged_data
end

function FindStatsToMerge(rsl_reg,loc_reg)
    dst_arr = ones(Float32,size(loc_reg,1),size(loc_reg,1)) * 100000
    for i ∈ 1:size(loc_reg,1), j ∈ i+1:size(loc_reg,1)
        if sum(@. isfinite(rsl_reg[i,:] * rsl_reg[j,:])) > 10
            dst_arr[i,j] = ComputeDistance(loc_reg[i,1],loc_reg[i,2],loc_reg[j,1],loc_reg[j,2])
        end
    end
    merge_idx = argmin(dst_arr)
    ϕ_mid, θ_mid = ComputeMidpoint(loc_reg[merge_idx[1],1],loc_reg[merge_idx[1],2],loc_reg[merge_idx[2],1],loc_reg[merge_idx[2],2])
    return merge_idx, ϕ_mid, θ_mid
end

function ComputeDistance(ϕ1,θ1,ϕ2,θ2)
    dist = 2*6371*asin(sqrt(sind(0.5*(θ2-θ1))^2+cosd(θ1)*cosd(θ2)*sind(0.5*(ϕ1-ϕ2))^2))
    return dist
end

function ComputeMidpoint(ϕ1,θ1,ϕ2,θ2)
    Bx = cos(deg2rad(θ2)) * cos(deg2rad(ϕ2 - ϕ1))
    By = cos(deg2rad(θ2)) * sin(deg2rad(ϕ2 - ϕ1))
    ϕ_mid = rad2deg(deg2rad(ϕ1)+atan(By, cos(deg2rad(θ1)) + Bx))
    θ_mid = rad2deg(atan(sin(deg2rad(θ1))+sin(deg2rad(θ2)), sqrt((cos(deg2rad(θ1))+Bx)^2 + By^2)))
    return ϕ_mid, θ_mid
end

function ComputeGridArea(lon,lat)
    gridsize = abs(lat[2]-lat[1])
    area = @. deg2rad(gridsize) * (sind(lat+gridsize/2)-sind(lat-gridsize/2)) * 6371000^2
    return repeat(area',size(lon,1))
end


end