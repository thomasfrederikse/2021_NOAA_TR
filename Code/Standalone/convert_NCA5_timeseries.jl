# --------------------------------------------
# Read an NCA5 file and convert it to a usable
# NetCDF file
# --------------------------------------------
using NetCDF
using DelimitedFiles
using Plots
using Statistics
using NCDatasets
using Interpolations

function main()
    fn_in = homedir()*"/Projects/2021_NOAA_TR/Data/NCA5_RSL_projections_vlm.csv"
    fn_out_gmsl = homedir()*"/Projects/2021_NOAA_TR/Data/NCA5_RSL_projections_vlm_gmsl.nc"
    fn_out_tg = homedir()*"/Projects/2021_NOAA_TR/Data/NCA5_RSL_projections_vlm_tg.nc"
    fn_out_grid = homedir()*"/Projects/2021_NOAA_TR/Data/NCA5_RSL_projections_vlm_grid.nc"

    scenarios = ["Low","IntLow","Int","IntHigh","High","Extreme"]

    data_raw = readdlm(fn_in,',',skipstart=16)
    traw = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120,2130,2140, 2150, 2200, 2250, 2300]

    gmsl_raw = data_raw[1:18,:]

    gmsl = Dict()
    gmsl["time"] = traw
    for (idx,scn) ∈ enumerate(scenarios)
        gmsl[scn] = convert.(Int,gmsl_raw[3*(idx-1)+1:3*idx,7:end]')
        @. gmsl[scn] = gmsl[scn][:,[2,1,3]]
    end
    
    # Tide gauge stations
    idx_tg = @. ((data_raw[:,2] > 0) & (data_raw[:,2] < 10000))
    tg_raw = data_raw[idx_tg,:]
    tg_store = Array{Dict}(undef,convert(Int,sum(idx_tg)/18))

    tg_store = Dict()
    tg_store["time"] = traw

    tg_store["name"] = Array{String}(undef,convert(Int,sum(idx_tg)/18))
    tg_store["id"] = Array{Int32}(undef,convert(Int,sum(idx_tg)/18))
    tg_store["θ"] = Array{Float32}(undef,convert(Int,sum(idx_tg)/18))
    tg_store["ϕ"] = Array{Float32}(undef,convert(Int,sum(idx_tg)/18))
    for (idx,scn) ∈ enumerate(scenarios)
        tg_store[scn] = zeros(Int32,convert(Int,sum(idx_tg)/18),length(traw),3);
    end

    for tg ∈ 1:convert(Int,sum(idx_tg)/18)
        tg_idx = [1+(18*(tg-1)):18*tg...]
        tg_raw_lcl = @views tg_raw[tg_idx,:]
        tg_store["name"][tg] = tg_raw_lcl[1,1]
        tg_store["id"][tg] = tg_raw_lcl[1,2]
        tg_store["θ"][tg] = tg_raw_lcl[1,3]
        tg_store["ϕ"][tg] = tg_raw_lcl[1,4]>0 ? tg_raw_lcl[1,4] : tg_raw_lcl[1,4]+360
        for (idx,scn) ∈ enumerate(scenarios)
            tg_store[scn][tg,:,:] = convert.(Int,tg_raw_lcl[3*(idx-1)+1:3*idx,7:end]')
            @. tg_store[scn][tg,:,:] = tg_store[scn][tg,:,:][:,[2,1,3]]
        end
    end

    ## Grid cells
    idx_grid = @. (data_raw[:,2] > 10000)
    grid_raw = @views data_raw[idx_grid,:]
    grid_store = Dict()
    grid_store["time"] = traw
    grid_store["θ"] = convert.(Int32,[-90:90...])
    grid_store["ϕ"] = convert.(Int32,[0:359...])
    for (idx,scn) ∈ enumerate(scenarios)
        grid_store[scn] = zeros(Int,360,181,length(traw),3);
    end

    for pnt ∈ 1:convert(Int,sum(idx_grid)/18)
        pnt_idx = [1+(18*(pnt-1)):18*pnt...]
        grid_raw_lcl = @views grid_raw[pnt_idx,:]
        idx_θ = argmin(@. abs(grid_store["θ"]-grid_raw_lcl[1,3]))
        ϕ = grid_raw_lcl[1,4]>0 ? grid_raw_lcl[1,4] : grid_raw_lcl[1,4]+360
        idx_ϕ = argmin(@. abs(grid_store["ϕ"]-ϕ))
        for (idx,scn) ∈ enumerate(scenarios)
            @. grid_store[scn][idx_ϕ,idx_θ,:,:] = grid_raw_lcl[3*(idx-1)+1:3*idx,7:end]'
            @. grid_store[scn][idx_ϕ,idx_θ,:,:] = grid_store[scn][idx_ϕ,idx_θ,:,[2,1,3]]
        end
    end

    ### Save data
    fh = Dataset(fn_out_gmsl,"c")
    defDim(fh,"time", length(gmsl["time"]))
    defDim(fh,"percentiles",3)
    defVar(fh,"time",convert.(Int32,gmsl["time"]),("time",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)
    for (idx,scn) ∈ enumerate(scenarios)
        defVar(fh,scn,convert.(Int32,gmsl[scn]),("time","percentiles",),deflatelevel=5)
    end
    close(fh)

    ## Grid
    fh = Dataset(fn_out_grid,"c")
    defDim(fh,"time", length(grid_store["time"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"lon", length(grid_store["ϕ"]))
    defDim(fh,"lat", length(grid_store["θ"]))
    defVar(fh,"lon",grid_store["ϕ"],("lon",),deflatelevel=5)
    defVar(fh,"lat",grid_store["θ"],("lat",),deflatelevel=5)
    defVar(fh,"time",convert.(Int32,grid_store["time"]),("time",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)
    for (idx,scn) ∈ enumerate(scenarios)
        defVar(fh,scn,convert.(Int32,grid_store[scn]),("lon","lat","time","percentiles",),deflatelevel=5)
    end
    close(fh)

    ## TG
    fh = Dataset(fn_out_tg,"c")
    defDim(fh,"time", length(tg_store["time"]))
    defDim(fh,"percentiles",3)
    defDim(fh,"id", length(tg_store["id"]))

    defVar(fh,"id",tg_store["id"],("id",),deflatelevel=5)
    defVar(fh,"name",tg_store["name"],("id",),deflatelevel=5)
    defVar(fh,"lat",tg_store["θ"],("id",),deflatelevel=5)
    defVar(fh,"lon",tg_store["ϕ"],("id",),deflatelevel=5)
    defVar(fh,"time",convert.(Int32,tg_store["time"]),("time",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)
    for (idx,scn) ∈ enumerate(scenarios)
        defVar(fh,scn,convert.(Int32,tg_store[scn]),("id","time","percentiles",),deflatelevel=5)
    end
    close(fh)
end
