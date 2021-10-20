module Hector
# ------------------------------------------------------------
# Hector tools
# File with all sorts of handy tools to use Hector from Julia
# Thomas ♡♡♡ Hector
# ------------------------------------------------------------
using DelimitedFiles
using Statistics
using JSON
function EstTrend(time,height;accel=false,model="White",AR=0,MA=0,SA=false,SSA=false,monthly=false)
    time = convert.(Float32,time)
    height = convert.(Float32,height)
    # Set time in MJD
    if monthly
        mjd = 30.5 * [1:length(time)...]
    else
        mjd = @. 365.25 * time - 678942
    end

    # Filenames
    dir_hector = homedir()*"/Code/Hector/hector_1.9_source/"
    fn_config = dir_hector*"estimatetrend.ctl"
    fn_input  = dir_hector*"input.mom"
    fn_output  = dir_hector*"estimatetrend.json"

    # Write input file
    if monthly
        headerline = "# sampling period "*string(30.5)*"\n"
    else
        headerline = "# sampling period "*string(minimum(diff(mjd)))*"\n"
    end
    f = open(fn_input, "w")
    write(f, headerline)
    writedlm(f, [mjd height])
    close(f)

    # Write configuration file
    config_list = String[]
    push!(config_list,"DataFile input.mom\n")
    push!(config_list,"DataDirectory ./\n")
    push!(config_list,"OutputFile trend.out\n")
    push!(config_list,"interpolate no\n")
    push!(config_list,"firstdifference no\n")
    push!(config_list,"estimateoffsets no\n")
    push!(config_list,"PhysicalUnit m\n")
    push!(config_list,accel ? "DegreePolynomial 2\n" : "DegreePolynomial 1\n")
    push!(config_list,SA ? "seasonalsignal yes\n" : "seasonalsignal no\n")
    push!(config_list,SSA ? "halfseasonalsignal yes\n" : "halfseasonalsignal no\n")
    if model == "GGM"
        push!(config_list,"NoiseModels GGM White\n")
    elseif model=="Powerlaw"
        push!(config_list,"NoiseModels Powerlaw White\n")
    elseif model=="ARMA"
        push!(config_list,"NoiseModels ARMA White\n")
    elseif model=="Matern"
        push!(config_list,"NoiseModels Matern White\n")
    elseif model=="White"
        push!(config_list,"NoiseModels White\n")
    else
        throw(ArgumentError(model* " is not a valid noise model"))
    end
    push!(config_list,"AR_p "*string(AR)*"\n");
    push!(config_list,"MA_q "*string(MA)*"\n");
    push!(config_list,"JSON yes\n");

    f = open(fn_config, "w");
    [write(f, line) for line in config_list];
    close(f);
    pdir = pwd();
    cd(dir_hector);
    run(pipeline(`./estimatetrend`; stdout=devnull, stderr=devnull));
    cd(pdir);
    output = JSON.parse(join(readlines(fn_output)));
    # rm(fn_config);
    # rm(fn_input);
    # rm(fn_output);
    return output;
end

function GenerateNoise(time,height,n_sims,n_tsteps;accel=false,model="White",AR=0,MA=0,SA=false,SSA=false,monthly=false)
    # Generate noise from time series and noise model. 
    time = convert.(Float32,time)
    height = convert.(Float32,height)

    # Set time and sampling period
    monthly ? mjd = 30.5 * [1:length(time)...] : mjd = @. 365.25 * time - 678942
    monthly ? sampling_period = 30.5 : sampling_period = minimum(diff(mjd))

    # First step: get noise parameters:
    trend_info = EstTrend(time,height;accel,model,AR,MA,SA,SSA,monthly)

    dir_hector = homedir()*"/Code/Hector/hector_1.9_source/"
    fn_config = dir_hector*"simulatenoise.ctl"

    config_list = String[]
    push!(config_list,"SimulationDir "*dir_hector*"\n")
    push!(config_list,"PhysicalUnit m\n")
    push!(config_list,"TimeNoiseStart 1\n")
    push!(config_list,"GGM_1mphi 0.5\n")
    push!(config_list,"SimulationLabel Noise\n")
    push!(config_list,"NumberOfSimulations 1\n")
    push!(config_list,"NumberOfPoints "*string(n_sims*n_tsteps)*"\n")
    push!(config_list,"SamplingPeriod "*string(sampling_period)*"\n")
    if model == "GGM"
        push!(config_list,"NoiseModels GGM White\n")
    elseif model=="Powerlaw"
        push!(config_list,"NoiseModels Powerlaw White\n")
    # elseif model=="ARMA"
    #     push!(config_list,"NoiseModels ARMA White\n")
    # elseif model=="Matern"
    #     push!(config_list,"NoiseModels Matern White\n")
    # elseif model=="White"
    #     push!(config_list,"NoiseModels White\n")
    else
        throw(ArgumentError(model* " is not a implemented noise model for noise creation "))
    end
    push!(config_list,"AR_p "*string(AR)*"\n");
    push!(config_list,"MA_q "*string(MA)*"\n");
    f = open(fn_config, "w");
    [write(f, line) for line in config_list];
    close(f);
    pdir = pwd();
    cd(dir_hector);
    if model == "GGM"
        val1=trend_info["driving_noise"]
        val2=trend_info["NoiseModel"]["GGM"]["fraction"]
        val3=trend_info["NoiseModel"]["GGM"]["fraction"]
        val4=trend_info["NoiseModel"]["GGM"]["d"]
        run(pipeline(`printf %".3f\n" $val1 $val2 $val3 $val4 `,`./simulatenoise`));

    elseif model=="Powerlaw"
        val1=trend_info["driving_noise"]
        val2=trend_info["NoiseModel"]["Powerlaw"]["fraction"]
        val3=trend_info["NoiseModel"]["White"]["fraction"]
        val4=trend_info["NoiseModel"]["Powerlaw"]["d"]
        run(pipeline(`printf %".3f\n" $val1 $val2 $val3 $val4 `,`./simulatenoise`));
    end
    noise = readdlm("Noise_0.mom",skipstart=1)[:,2]
    rm(fn_config);
    rm("Noise_0.mom");
    cd(pdir);
    noise = reshape(noise,(n_tsteps,n_sims))
    return noise
end
end