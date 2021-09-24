module Masks
using NetCDF
using NCDatasets
function CreateMask(settings)
    # ----------------------------------------------------------------------------
    # Make a 0.5° mask of US coastal sea level to combine with CDS/CMEMS altimetry
    # 8 basins:
    # 1. EC East Coast
    # 2. SE southeast coast
    # 3. GCE Gulf Coast east
    # 4. GCW Gulf coast west
    # 5. SWC Southwest coast
    # 6. NWC Northwest coast
    # 7. PAC Pacific Islands
    # 8. CAR Carribean Islands
    # ----------------------------------------------------------------------------
    fn_bathy = homedir()*"/Data/GEBCO/grid_05.nc"
    fn_slm   = homedir()*"/Data/GRACE/JPL_mascon/mask.nc"

    # boundaries
    EC = [282.5 292;35.2  43.0]
    SE = [278.5 286; 25.0  35.2]
    GCE = [ 270.5  278.5; 25.0  31.0]
    GCW = [261.0  270.5; 25.0  31.0]
    SWC = [234.0  242.0; 31.0  42.0]
    NWC = [234.0  242.0; 42.0  48.0] 
    PAC = [143.0  230.0; -20.0  32.0] 
    CAR = [280.0  300.0; 14.0  25.0] 

    # EC- Lat: [35.2 43]; Lon: [282.5 292]
    # SE- Lat: [25 35.2]; Lon: [278.5 286]
    # GC- Lat: [25 31]; Lon: [261 278.5]
    # WC- Lat: [31 47.5]; Lon: [234, 242]

    ϕ = reshape(ncread(fn_slm,"lon"),(:,1))
    θ = reshape(ncread(fn_slm,"lat"),(1,:))
    slm = convert.(Bool,1 .- ncread(fn_slm,"land"))
    depth = ncread(fn_bathy,"z")

    EC_mask = @. (ϕ>=EC[1,1]) & (ϕ<=EC[1,2]) & (θ>=EC[2,1]) & (θ<=EC[2,2]) & slm & (depth > -500)
    SE_mask = @. (ϕ>=SE[1,1]) & (ϕ<=SE[1,2]) & (θ>=SE[2,1]) & (θ<=SE[2,2]) & slm & (depth > -700)
    GCE_mask = @. (ϕ>=GCE[1,1]) & (ϕ<=GCE[1,2]) & (θ>=GCE[2,1]) & (θ<=GCE[2,2]) & slm & (depth > -2000)
    GCW_mask = @. (ϕ>=GCW[1,1]) & (ϕ<=GCW[1,2]) & (θ>=GCW[2,1]) & (θ<=GCW[2,2]) & slm & (depth > -2000)
    SWC_mask = @. (ϕ>=SWC[1,1]) & (ϕ<=SWC[1,2]) & (θ>=SWC[2,1]) & (θ<=SWC[2,2]) & slm & (depth > -2500)
    NWC_mask = @. (ϕ>=NWC[1,1]) & (ϕ<=NWC[1,2]) & (θ>=NWC[2,1]) & (θ<=NWC[2,2]) & slm & (depth > -2500)
    CAR_mask = @. (ϕ>=CAR[1,1]) & (ϕ<=CAR[1,2]) & (θ>=CAR[2,1]) & (θ<=CAR[2,2]) & slm & (depth > -2500)
    PAC_mask = @. (ϕ>=PAC[1,1]) & (ϕ<=PAC[1,2]) & (θ>=PAC[2,1]) & (θ<=PAC[2,2]) & slm 
    USA_mask = @. EC_mask | SE_mask | GCE_mask | GCW_mask | SWC_mask | NWC_mask

    # Mask without land removed (for tide-gauge selection)
    EC_unc = @. (ϕ>=EC[1,1]) & (ϕ<=EC[1,2]) & (θ>=EC[2,1]) & (θ<=EC[2,2])
    SE_unc = @. (ϕ>=SE[1,1]) & (ϕ<=SE[1,2]) & (θ>=SE[2,1]) & (θ<=SE[2,2]) 
    GCE_unc = @. (ϕ>=GCE[1,1]) & (ϕ<=GCE[1,2]) & (θ>=GCE[2,1]) & (θ<=GCE[2,2]) 
    GCW_unc = @. (ϕ>=GCW[1,1]) & (ϕ<=GCW[1,2]) & (θ>=GCW[2,1]) & (θ<=GCW[2,2]) 
    SWC_unc = @. (ϕ>=SWC[1,1]) & (ϕ<=SWC[1,2]) & (θ>=SWC[2,1]) & (θ<=SWC[2,2]) 
    NWC_unc = @. (ϕ>=NWC[1,1]) & (ϕ<=NWC[1,2]) & (θ>=NWC[2,1]) & (θ<=NWC[2,2]) 
    CAR_unc = @. (ϕ>=CAR[1,1]) & (ϕ<=CAR[1,2]) & (θ>=CAR[2,1]) & (θ<=CAR[2,2]) 
    PAC_unc = @. (ϕ>=PAC[1,1]) & (ϕ<=PAC[1,2]) & (θ>=PAC[2,1]) & (θ<=PAC[2,2])
    USA_unc = @. EC_unc | SE_unc | GCE_unc | GCW_unc | SWC_unc | NWC_unc

    # Save
    fh = Dataset(settings["fn_region_mask"],"c")
    defDim(fh,"lon", length(ϕ[:]))
    defDim(fh,"lat", length(θ[:]))
    defVar(fh,"lon",ϕ[:],("lon",),deflatelevel=5)
    defVar(fh,"lat",θ[:],("lat",),deflatelevel=5)
    defVar(fh,"EC_mask",convert.(UInt8,EC_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"SE_mask",convert.(UInt8,SE_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"GCE_mask",convert.(UInt8,GCE_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"GCW_mask",convert.(UInt8,GCW_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"SWC_mask",convert.(UInt8,SWC_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"NWC_mask",convert.(UInt8,NWC_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"CAR_mask",convert.(UInt8,CAR_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"PAC_mask",convert.(UInt8,PAC_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"USA_mask",convert.(UInt8,USA_mask),("lon","lat",),deflatelevel=5)

    defVar(fh,"EC_unc",convert.(UInt8,EC_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"SE_unc",convert.(UInt8,SE_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"GCE_unc",convert.(UInt8,GCE_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"GCW_unc",convert.(UInt8,GCW_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"SWC_unc",convert.(UInt8,SWC_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"NWC_unc",convert.(UInt8,NWC_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"CAR_unc",convert.(UInt8,CAR_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"PAC_unc",convert.(UInt8,PAC_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"USA_unc",convert.(UInt8,USA_unc),("lon","lat",),deflatelevel=5)
    close(fh)
    return nothing
end

function ReadMask(settings)
    mask = Dict()
    mask["ϕ"] = ncread(settings["fn_region_mask"],"lon")
    mask["θ"] = ncread(settings["fn_region_mask"],"lat")
    for region ∈ settings["regions"]
        mask[region] = convert.(Bool,ncread(settings["fn_region_mask"],region*"_mask"))
        mask[region*"_unc"] = convert.(Bool,ncread(settings["fn_region_mask"],region*"_unc"))
    end
    return mask
end

end