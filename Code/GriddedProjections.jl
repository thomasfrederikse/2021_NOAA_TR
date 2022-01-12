# ------------------------------------------------------------------------------
# Read converted NCA5 files and add 15 mm baseline correction to account for the
# 2000 - 2005 difference between both
# ------------------------------------------------------------------------------
module GriddedProjections
using NetCDF
using NCDatasets

function RunGriddedProjections(settings)
    println("\nGridded scenarios...")
    # Locations and time
    fn = settings["dir_NCA5"]*"NCA5_Low_grid.nc"
    lon = ncread(fn,"lon")
    lat = ncread(fn,"lat")
    years_NCA5 = ncread(fn,"years",start=[1],count=[14])

    # Create output NetCDF file
    fh = Dataset(settings["fn_proj_gri"],"c")
    defDim(fh,"years", length(years_NCA5))
    defDim(fh,"percentiles",3)
    defDim(fh,"lon", length(lon))
    defDim(fh,"lat", length(lat))

    defVar(fh,"lon",lon,("lon",),deflatelevel=5)
    defVar(fh,"lat",lat,("lat",),deflatelevel=5)
    defVar(fh,"years",years_NCA5,("years",),deflatelevel=5)
    defVar(fh,"percentiles",convert.(Int32,[17,50,83]),("percentiles",),deflatelevel=5)

    # Write scenarios
    for scn in settings["NCA5_scenarios"]
        fn_in = settings["dir_NCA5"]*"NCA5_"*scn*"_grid.nc"
        for prc in settings["processes"]
            NCA5_prc = ncread(fn_in,prc,start=[1,1,1,1],count=[-1,-1,14,-1]);
            if (prc == "total") | (prc == "verticallandmotion")
                NCA5_prc .+= 15
            end
            defVar(fh,"rsl_"*prc*"_"*scn,NCA5_prc,("lon","lat","years","percentiles",),deflatelevel=5)
        end
    end
    close(fh)
    println("Gridded scenarios done\n")
end

end