# --------------------------------------------------------------------
# ProcessObservations: all functions to read and post-process
# observations from altimetry and tide gauges
#   Run a virtual-station average for US tide-gauge observations
#   Read altimetry observations
#   Return annual-mean reconstructions
# These routines are called from LocalProjections and 
# RegionalProjections
# --------------------------------------------------------------------

module ProcessObservations
using NCDatasets
using NetCDF
using Dates
using XLSX
using DelimitedFiles
using Statistics
using MAT
using Interpolations
using CSV
dir_code = homedir()*"/Projects/2021_NOAA_TR/Code/"
include(dir_code*"Masks.jl")
include(dir_code*"Hector.jl")

function ReadLocalObs(settings)
    # -------------------------------------------------------------------------------------------
    # Read tide-gauge data and return list of local annual-mean tide-gauge data, 
    # the station names and coordinates
    #
    # To compute annual-mean sea level from monthly-mean sea level:
    # - First remove mean seasonal cycle, which would give trouble for years with missing months
    # - For years with at least 10 months with data, compute annual mean. Otherwise, NaN32
    # -------------------------------------------------------------------------------------------
    println("  Reading tide-gauge data...")
    # Read data file
    tg_data_raw = XLSX.readxlsx(settings["fn_tg_data"])["monthly MSL"][:]


    
    # Station name, lon,lat coords. Longitude in [0 360] interval
    station_names = Vector{String}(tg_data_raw[1,2:122])
    station_coords = reverse!(Array{Float32}(tg_data_raw[2:3,2:122]'),dims=2)
    @. station_coords[:,1] = mod(station_coords[:,1],360)
    
    # Time steps
    tsteps_raw = tg_data_raw[4:end,1]
    tvec = zeros(Int32,length(tsteps_raw),2)
    for t in 1:length(tsteps_raw)
        tvec[t,1] = Year(tsteps_raw[t]).value
        tvec[t,2] = Month(tsteps_raw[t]).value
    end
    tval = @. tvec[:,1] + (tvec[:,2]/12 - 1/24)

    # Monthly-mean RSL
    rsl_raw = @views tg_data_raw[4:end,2:122]
     @. rsl_raw[rsl_raw=="NaN"] = NaN
    rsl_monthly = Array{Float32}(rsl_raw')
    rsl_monthly .*= 1000.0f0 # From meters to millimeters

    # Remove station Seward because of large vertical datum discontinuity
    stat_acc = ones(Bool,size(rsl_monthly,1))
    stat_acc[107] = false
    station_names = station_names[stat_acc]
    station_coords = station_coords[stat_acc,:]
    rsl_monthly = rsl_monthly[stat_acc,:]

    # From monthly to annual-mean sea level
    local_obs = Dict()
    local_obs["name"]  = station_names
    local_obs["coords"] = station_coords
    local_obs["rsl"]            = zeros(Float32,size(rsl_monthly,1),length(settings["years_tg"])) .* NaN32
    for stat in 1:length(station_names)
        # 1. Remove mean seasonal cycle before computing annual means
        acc_mnth = @. isfinite(rsl_monthly[stat,:])
        amat = ones(Float32,sum(acc_mnth),6)
        amat[:,2] = tval[acc_mnth] .- mean(tval[acc_mnth])
        amat[:,3] = @. sin(2*π*tval[acc_mnth])
        amat[:,4] = @. cos(2*π*tval[acc_mnth])
        amat[:,5] = @. sin(4*π*tval[acc_mnth])
        amat[:,6] = @. cos(4*π*tval[acc_mnth])
        sol = amat\rsl_monthly[stat,acc_mnth]    
        sol[2] = 0
        rsl_monthly[stat,acc_mnth]  -= amat*sol
        rsl_rfmt = reshape(rsl_monthly[stat,:],(12,length(settings["years_tg"])))
        for yr in 1:length(settings["years_tg"])
            sum(isnan.(rsl_rfmt[:,yr])) < 3 ? local_obs["rsl"][stat,yr] = mean(filter(isfinite,rsl_rfmt[:,yr])) : nothing
        end
    end
    local_obs["region"] = MapToRegions(local_obs,settings)

    # Find nearest PSMSL ID (For Nga)
    psmsl_data = readdlm(settings["dir_project"] * "Data/filelist_psmsl.txt" ,';')
    lon_psmsl = mod.(convert.(Float32,psmsl_data[:,3]),360)
    lat_psmsl = convert.(Float32,psmsl_data[:,2])
    # name_psmsl = psmsl_data[:,4]
    id_psmsl = convert.(Int32,psmsl_data[:,1])
    local_obs["psmsl_id"] = zeros(Int,length(local_obs["name"]))
    for tg in 1:length(local_obs["name"])
        dst_arr = (local_obs["coords"][tg,1].-lon_psmsl).^2 + (local_obs["coords"][tg,2].-lat_psmsl).^2
        if minimum(dst_arr) .< 0.001
            local_obs["psmsl_id"][tg] = id_psmsl[argmin(dst_arr)]
            # println(rpad(local_obs["name"][tg],40," ")*" "*name_psmsl[argmin(dst_arr)])
        else
            local_obs["psmsl_id"][tg] = -1
            # println("NOMATCH: "*rpad(local_obs["name"][tg],40," ")*" "*name_psmsl[argmin(dst_arr)])
        end
    end

    # Flag Chesapeake bridge and La Jolla because of possible data problems
    local_obs["qc_flag"] = zeros(Bool,length(local_obs["name"]))
    local_obs["qc_flag"][findfirst(==(256),local_obs["psmsl_id"])] = true
    local_obs["qc_flag"][findfirst(==(1635),local_obs["psmsl_id"])] = true

    return local_obs
end

function MapToRegions(local_obs, settings)
    # --------------------------------------------------------------
    # Map the individual tide-gauge records to one of the 11 regions
    # Return a list with region numbers
    # --------------------------------------------------------------
    println("  Mapping tide-gauge stations onto regions...")
    mask = Masks.ReadMask(settings)
    region_num = zeros(Bool,length(settings["regions"]),size(local_obs["coords"],1))
    
    # Find Alaska stations
    in_alaska = zeros(Bool,size(local_obs["coords"],1))
    for stat ∈ 1:size(local_obs["coords"],1)
        strip(local_obs["name"][stat][end-2:end]) == "AK" ? in_alaska[stat] = true : nothing
    end
    NAL_list = ["Unalaska, AK","Prudhoe Bay, AK","Nome, AK","Adak Island, AK ","Unalaska, AK"]
    # Determine region
    for (region_idx,region) ∈ enumerate(settings["regions"])
        if (region!="NAL") & (region!="SAL")
            for stat ∈ 1:size(local_obs["coords"],1)
                ϕ_idx = argmin(@. abs(local_obs["coords"][stat,1] - mask["ϕ"]))
                θ_idx = argmin(@. abs(local_obs["coords"][stat,2] - mask["θ"]))
                if mask[region*"_unc"][ϕ_idx,θ_idx]
                    region_num[region_idx,stat] = true
                end
            end
        elseif region=="NAL"
            for stat ∈ 1:size(local_obs["coords"],1)
                local_obs["name"][stat] in NAL_list ? region_num[region_idx,stat] = true : nothing
            end
        elseif region=="SAL"
            for stat ∈ 1:size(local_obs["coords"],1)
                in_alaska[stat] & !region_num[region_idx-1,stat] ? region_num[region_idx,stat] = true : nothing
            end
        end
    end
    return region_num
end

function ComputeRegionObs(local_obs,settings)
    # -------------------------------------------------------------
    # Read tide-gauge observations and merge using virtual stations
    # Read altimetry observations using region mask
    # Return both curves per region
    # -------------------------------------------------------------
    println("  Computing regional observations...")
    region_obs = Dict()
    [region_obs[region] = Dict() for region in settings["regions"]]

    ComputeVirstat!(local_obs,region_obs,settings) 
    ReadAltimetry!(region_obs,settings) # Altimetry

    # Contiguous USA curve is average of the six regions that circumference the Contiguous US
    region_obs["USA"]["rsl"] = (region_obs["NE"]["rsl"] + region_obs["SE"]["rsl"] + region_obs["WGOM"]["rsl"] + region_obs["EGOM"]["rsl"] + region_obs["SW"]["rsl"] + region_obs["NW"]["rsl"]) ./6
    
    # Remove natural variability
    clim_indices = ReadClimIndices(settings)
    for region in settings["regions"]
        region_obs[region]["rsl"] = RemoveVariability(settings["years_tg"],region_obs[region]["rsl"],clim_indices)
        region_obs[region]["alt"] = RemoveVariability(settings["years_tg"],region_obs[region]["alt"],clim_indices)
    end
    return region_obs    
end

function ComputeVirstat!(local_obs,region_obs,settings)
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
    for (region_idx, region) in enumerate(settings["regions"])
        # List with locations and rsl data
        rsl_reg   = local_obs["rsl"][local_obs["region"][region_idx,:],:]
        coord_reg = local_obs["coords"][local_obs["region"][region_idx,:],:]
        while size(coord_reg,1) > 1
            # Find two nearest stations
            merge_idx, ϕ_mid, θ_mid = FindStatsToMerge(rsl_reg,coord_reg)
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
            coord_reg = coord_reg[keep_idx,:]
            rsl_reg = vcat(rsl_reg, merged_data')
            coord_reg = vcat(coord_reg, [ϕ_mid, θ_mid]')
        end
        region_obs[region]["rsl"] = rsl_reg[:]
    end
    return nothing
end

function ReadAltimetry!(region_obs,settings)
    # -----------------------------------------------------------------------------------
    # Read satellite altimetry and add GIA and GRD corrections to obtain
    # relative sea level, which is comparable to tide-gauge observations
    # 
    # Altimetry observations ahve been obtained from Copernicus Climate Data Store
    # https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-sea-level-global
    # This dataset contains modified Copernicus Climate Change Service information [2020]
    # -----------------------------------------------------------------------------------
    println("  Reading altimetry...")
    mask = Masks.ReadMask(settings) 
    area = ComputeGridArea(mask["ϕ"],mask["θ"])

    # Read altimetry data
    ssh_grid = convert.(Float32,ncread(settings["fn_altimetry"],"ssh"))
    @. ssh_grid[ssh_grid > 31000] = 0
    @. ssh_grid .*= 0.05
    alt_time = convert.(Float32,ncread(settings["fn_altimetry"],"time"))
    alt_year = [1993:2020...]
    alt_acc  = findall(in(alt_year),settings["years_tg"])

    # From high frquency data to annual-means
    ssh_grid_year = zeros(Float32,size(ssh_grid,1),size(ssh_grid,2),length(alt_year));
    for (idx,yr) ∈ enumerate(alt_year)
        acc_idx = floor.(alt_time) .== yr
        ssh_grid_year[:,:,idx] = mean(ssh_grid[:,:,acc_idx],dims=3)
    end

    # Take regional average 
    for region ∈ settings["regions"]
        region_obs[region]["alt"] = zeros(Float32,length(settings["years_tg"])) .* NaN32
        region_obs[region]["alt"][alt_acc] = sum( @.(ssh_grid_year * mask[region] * area),dims=(1,2))[:] ./ sum(mask[region] .* area)
    end

    # 3 Add GIA and contemporary GRD
    GIA = ComputeGIA(mask,area,settings) # GIA from Caron et al. (2019)
    GRD = ComputeGRD(mask,area,alt_year,settings) # Contemporary GRD effects from Frederikse e al. (2020)
    for (region_idx,region) ∈ enumerate(settings["regions"])
        region_obs[region]["alt"][alt_acc] = region_obs[region]["alt"][alt_acc] - (GIA[region_idx] .* (alt_year .- mean(alt_year))) - GRD[:,region_idx]
    end
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

function RemoveVariability(years,tseries,clim_indices)
    # ---------------------------------------------------
    # Remove the effects of climate indices (NAO/PDO/MEI)
    # from the sea level curves using OLS. 
    # ---------------------------------------------------
    acc_idx = @. isfinite(tseries)
    amat = ones(Float32,sum(acc_idx),5);
    amat[:,2] = years[acc_idx].- mean(years[acc_idx]);
    amat[:,3] = LinearInterpolation(clim_indices[1]["years"], clim_indices[1]["index"],extrapolation_bc=Line())(years[acc_idx])
    amat[:,4] = LinearInterpolation(clim_indices[2]["years"], clim_indices[2]["index"],extrapolation_bc=Line())(years[acc_idx])
    amat[:,5] = LinearInterpolation(clim_indices[3]["years"], clim_indices[3]["index"],extrapolation_bc=Line())(years[acc_idx])
    sol = amat\tseries[acc_idx]
    sol[2] = 0
    tseries_var = zeros(Float32,length(tseries)) .* NaN32
    tseries_var[acc_idx] = amat*sol
    tseries_novar = tseries - tseries_var
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