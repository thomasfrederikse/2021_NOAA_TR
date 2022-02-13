module Masks

using NetCDF
using NCDatasets
using Interpolations

function CreateMask(settings)
    println("  Creating region masks...")

    # ----------------------------------------------------------------------------
    # Make a 0.5° mask of US coastal sea level to combine with CDS/CMEMS altimetry
    # 8 basins:
    # 1. NE Northeast Coast
    # 2. SE southeast coast
    # 3. EGOM Eastern Gulf of Mexico
    # 4. WGOM Gulf coast west
    # 5. SW Southwest coast
    # 6. NW Northwest coast
    # 7. PAC Pacific Islands (Hawaii)
    # 8. CAR Carribean Islands
    # 9. NAL North Alaska 
    # 10. SAL South Alaska
    # ----------------------------------------------------------------------------
    fn_slm   = homedir()*"/Data/GRACE/JPL_mascon/mask.nc"

    # boundaries
    NE = [282.5 292;35.2  43.0]
    SE = [278.5 286; 25.0  35.2]
    EGOM = [ 270.5  278.5; 25.0  31.0]
    WGOM = [261.0  270.5; 25.0  31.0]
    SW = [234.0  242.0; 31.0  42.0]
    NW = [234.0  242.0; 42.0  48.0] 
    PAC = [180.0  210.0; 18.0  30.0] 
    CAR = [280.0  300.0; 14.0  25.0] 
    NAL = [187.0  219.0; 51.0  72.0] 
    SAL1 = [172.0  213.0; 47.0  60.0] 
    SAL2 = [200.0  231.0; 54.75  62.5] 

    # NE- Lat: [35.2 43]; Lon: [282.5 292]
    # SE- Lat: [25 35.2]; Lon: [278.5 286]
    # GC- Lat: [25 31]; Lon: [261 278.5]
    # WC- Lat: [31 47.5]; Lon: [234, 242]

    ϕ = reshape(ncread(fn_slm,"lon"),(:,1))
    θ = reshape(ncread(fn_slm,"lat"),(1,:))
    slm = convert.(Bool,1 .- ncread(fn_slm,"land"))
    depth = ncread(settings["fn_bathymetry"],"z")
    basin = ncread(settings["fn_basin_codes"],"basin")

    NE_mask = @. (ϕ>=NE[1,1]) & (ϕ<=NE[1,2]) & (θ>=NE[2,1]) & (θ<=NE[2,2]) & slm & (depth > -500)
    SE_mask = @. (ϕ>=SE[1,1]) & (ϕ<=SE[1,2]) & (θ>=SE[2,1]) & (θ<=SE[2,2]) & slm & (depth > -700)
    EGOM_mask = @. (ϕ>=EGOM[1,1]) & (ϕ<=EGOM[1,2]) & (θ>=EGOM[2,1]) & (θ<=EGOM[2,2]) & slm & (depth > -2000)
    WGOM_mask = @. (ϕ>=WGOM[1,1]) & (ϕ<=WGOM[1,2]) & (θ>=WGOM[2,1]) & (θ<=WGOM[2,2]) & slm & (depth > -2000)
    SW_mask = @. (ϕ>=SW[1,1]) & (ϕ<=SW[1,2]) & (θ>=SW[2,1]) & (θ<=SW[2,2]) & slm & (depth > -2500)
    NW_mask = @. (ϕ>=NW[1,1]) & (ϕ<=NW[1,2]) & (θ>=NW[2,1]) & (θ<=NW[2,2]) & slm & (depth > -2500)
    CAR_mask = @. (ϕ>=CAR[1,1]) & (ϕ<=CAR[1,2]) & (θ>=CAR[2,1]) & (θ<=CAR[2,2]) & slm & (depth > -2500)
    PAC_mask = @. (ϕ>=PAC[1,1]) & (ϕ<=PAC[1,2]) & (θ>=PAC[2,1]) & (θ<=PAC[2,2]) & slm & (depth > -5000)
    NAL_mask = @. (ϕ>=NAL[1,1]) & (ϕ<=NAL[1,2]) & (θ>=NAL[2,1]) & (θ<=NAL[2,2]) & slm & (depth > -500) & (basin !=1) & (basin !=0)

    SAL1_mask = @. (ϕ>=SAL1[1,1]) & (ϕ<=SAL1[1,2]) & (θ>=SAL1[2,1]) & (θ<=SAL1[2,2]) & slm & (depth > -500) & (basin ==1) 
    SAL2_mask = @. (ϕ>=SAL2[1,1]) & (ϕ<=SAL2[1,2]) & (θ>=SAL2[2,1]) & (θ<=SAL2[2,2]) & slm & (depth > -500) & (basin ==1) 
    SAL_mask = @. SAL1_mask | SAL2_mask
    USA_mask = @. NE_mask | SE_mask | EGOM_mask | WGOM_mask | SW_mask | NW_mask

    # Mask without land removed (for tide-gauge selection)
    NE_unc = @. (ϕ>=NE[1,1]) & (ϕ<=NE[1,2]) & (θ>=NE[2,1]) & (θ<=NE[2,2])
    SE_unc = @. (ϕ>=SE[1,1]) & (ϕ<=SE[1,2]) & (θ>=SE[2,1]) & (θ<=SE[2,2]) 
    EGOM_unc = @. (ϕ>=EGOM[1,1]) & (ϕ<=EGOM[1,2]) & (θ>=EGOM[2,1]) & (θ<=EGOM[2,2]) 
    WGOM_unc = @. (ϕ>=WGOM[1,1]) & (ϕ<=WGOM[1,2]) & (θ>=WGOM[2,1]) & (θ<=WGOM[2,2]) 
    SW_unc = @. (ϕ>=SW[1,1]) & (ϕ<=SW[1,2]) & (θ>=SW[2,1]) & (θ<=SW[2,2]) 
    NW_unc = @. (ϕ>=NW[1,1]) & (ϕ<=NW[1,2]) & (θ>=NW[2,1]) & (θ<=NW[2,2]) 
    CAR_unc = @. (ϕ>=CAR[1,1]) & (ϕ<=CAR[1,2]) & (θ>=CAR[2,1]) & (θ<=CAR[2,2]) 
    PAC_unc = @. (ϕ>=PAC[1,1]) & (ϕ<=PAC[1,2]) & (θ>=PAC[2,1]) & (θ<=PAC[2,2])
    NAL_unc = @. (ϕ>=NAL[1,1]) & (ϕ<=NAL[1,2]) & (θ>=NAL[2,1]) & (θ<=NAL[2,2])
    SAL1_unc = @. (ϕ>=SAL1[1,1]) & (ϕ<=SAL1[1,2]) & (θ>=SAL1[2,1]) & (θ<=SAL1[2,2])
    SAL2_unc = @. (ϕ>=SAL2[1,1]) & (ϕ<=SAL2[1,2]) & (θ>=SAL2[2,1]) & (θ<=SAL2[2,2])
    SAL_unc = @. SAL1_unc | SAL2_unc

    USA_unc = @. NE_unc | SE_unc | EGOM_unc | WGOM_unc | SW_unc | NW_unc

    # Make numerical mask for plot
    mask_num = zeros(Int16,size(NE_mask))
    mask_num[NE_mask] .= 1
    mask_num[SE_mask] .= 2
    mask_num[EGOM_mask] .= 3
    mask_num[WGOM_mask] .= 4
    mask_num[SW_mask] .= 5
    mask_num[NW_mask] .= 6
    mask_num[PAC_mask] .= 7
    mask_num[CAR_mask] .= 8
    mask_num[NAL_mask] .= 9
    mask_num[SAL_mask] .= 10

    # Save
    println("   Saving...")

    fh = Dataset(settings["fn_region_mask"],"c")
    defDim(fh,"lon", length(ϕ[:]))
    defDim(fh,"lat", length(θ[:]))
    defVar(fh,"lon",ϕ[:],("lon",),deflatelevel=5)
    defVar(fh,"lat",θ[:],("lat",),deflatelevel=5)
    defVar(fh,"NE_mask",convert.(UInt8,NE_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"SE_mask",convert.(UInt8,SE_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"EGOM_mask",convert.(UInt8,EGOM_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"WGOM_mask",convert.(UInt8,WGOM_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"SW_mask",convert.(UInt8,SW_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"NW_mask",convert.(UInt8,NW_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"PAC_mask",convert.(UInt8,PAC_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"CAR_mask",convert.(UInt8,CAR_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"NAL_mask",convert.(UInt8,NAL_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"SAL_mask",convert.(UInt8,SAL_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"USA_mask",convert.(UInt8,USA_mask),("lon","lat",),deflatelevel=5)
    defVar(fh,"mask_num",mask_num,("lon","lat",),deflatelevel=5)

    defVar(fh,"NE_unc",convert.(UInt8,NE_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"SE_unc",convert.(UInt8,SE_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"EGOM_unc",convert.(UInt8,EGOM_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"WGOM_unc",convert.(UInt8,WGOM_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"SW_unc",convert.(UInt8,SW_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"NW_unc",convert.(UInt8,NW_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"CAR_unc",convert.(UInt8,CAR_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"PAC_unc",convert.(UInt8,PAC_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"NAL_unc",convert.(UInt8,NAL_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"SAL_unc",convert.(UInt8,SAL_unc),("lon","lat",),deflatelevel=5)
    defVar(fh,"USA_unc",convert.(UInt8,USA_unc),("lon","lat",),deflatelevel=5)
    close(fh)
    println("  Creating region done")

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

function RegridBasinCodes(settings)
    # --------------------------------------------------------------
    # Regrid the basin code map from Eric Leuliette/NOAA
    # Data provided by the NOAA Laboratory for Satellite Altimetry."
    # --------------------------------------------------------------
    fn_in = homedir()*"/Data/Basins/basin_codes.nc"
    # Read
    ncinfo(fn_in)
    θ_in = ncread(fn_in,"lat")
    ϕ_in = ncread(fn_in,"lon")
    B_in = 1.0f0 .* ncread(fn_in,"basin")

    # Interpolate
    ϕ = [0.25f0:0.5f0:359.75f0...]
    θ = [-89.75f0:0.5f0:89.75f0...]
    B = interpolate((ϕ_in,θ_in), B_in, Gridded(Constant()))(ϕ,θ)

    # Write
    fh = Dataset(settings["fn_basin_codes"],"c")
    defDim(fh,"lon", length(ϕ))
    defDim(fh,"lat", length(θ))
    defVar(fh,"lon",ϕ,("lon",),deflatelevel=5)
    defVar(fh,"lat",θ,("lat",),deflatelevel=5)
    defVar(fh,"basin",convert.(UInt8,B),("lon","lat",),deflatelevel=5)
    fh.attrib["title"] = ["BASIN5 - Global 5-minute Basin Grid regridded to 0.5 degree"]
    fh.attrib["data_source"] = ["Data provided by the NOAA Laboratory for Satellite Altimetry, regridded afterwards"]
    fh.attrib["author"] = ["Original data: Eric Leuliette NOAA, regridded by Thomas Frederikse NASA JPL/Caltech"]
    close(fh)
    return nothing
end

end