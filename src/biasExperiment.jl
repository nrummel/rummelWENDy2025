## Load everything we need 
includet("util.jl")
using ArgParse, ProgressMeter
##
save_dir = joinpath(DATA_DIR, "biasExperimentSubResults")
if !isdir(save_dir)
    mkpath(save_dir)
end  
function _getGoodwinExample()
    wstar = [3.4884, 0.0969, 2.15, 10, 0.0969, 0.0581, 0.0969, 0.0775]
    function _f!(du,u,p,t)
        du[1] = wstar[3] / (wstar[3] + u[3]^p[1]) - wstar[2] * u[1]
        du[2] = p[2]*u[1]- wstar[6]*u[2]
        du[3] = wstar[7]*u[2]-wstar[8]*u[3]
        nothing
    end
    trng = (0.0, 80.0)            
    initCond = [0.3617, 0.9137, 1.3934]
    prng = [
        (5.0, 30.0), # w4 hill coeffiecient
        (0.0,0.2),  # w5
    ]
    params = WENDyParameters( 
        radiusMinTime  = 0.5,
        radiusMaxTime  = 10.0,
        Kmax           = 500
    )
    pstar = [wstar[4], wstar[6]] 

    return SimulatedWENDyData(
        "goodwin", 
        _f!,
        initCond,
        trng,
        pstar,
        prng,
        params;
        linearInParameters=Val(false), # NONLINEAR!
        noiseDist=Val(LogNormal),
        dt=0.1
    );
end

function _goodwinExperiment(nr::AbstractFloat, mc::Int)
    ex = _getGoodwinExample()
    p0 = [
        33.08,  # w4
        0.1287,  # w5
    ]
    w1Min, w1Max = 5, 20
    w2Min, w2Max = 0.04, 0.1
    dw1=0.01
    dw2=0.001
    ww1 = w1Min:dw1:w1Max
    ww2 = w2Min:dw2:w2Max 
    @info "Running goodwin experiment monte carlo $mc and noise $nr"
    dt = @elapsed begin
        _saveFile = joinpath(save_dir,"goodwinCostFunData_nr.$(nr)_mc.$(mc).jld2")
        simParams = SimulationParameters(
            noiseRatio=nr,
            seed=mc,
            timeSubsampleRate=1, 
            corruptU0=false
        )
        algoNames =[:fslsq, :trustRegion]
        @info "  solving wendyProb..."
        wendyProb, w0, results = runExample(ex, simParams, algoNames; w0=p0, ll=Info)
        wnll = zeros(length(ww2), length(ww1))
        pbar = Progress(prod(size(wnll)))
        for (i,w1) in enumerate(ww1), (j,w2) in enumerate(ww2)
            wnll[j,i]  = try 
                wendyProb.wnll.f([w1,w2])
            catch 
                1e6
            end
            next!(pbar, desc="wnll...")
        end
        @info "  Saving data to $_saveFile..."
        JLD2.save(
            _saveFile, 
            Dict(
                "what"=>results.trustRegion.what, 
                "wnll"=>wnll
            )
        )
    end
    days = floor(dt/86400)
    dt -= days * 86400
    hr = floor(dt/3600)
    dt -= hr*3600
    minutes = floor(dt/60)
    dt -= minutes*60
    @info "  This took $days days, $hr hr, $minutes min, $(@sprintf "%.4g" dt) s"
end

function _goodwinSummary()
    p0 = [
        33.08,  # w4
        0.1287,  # w5
    ]
    w1Min, w1Max = 5, 20
    w2Min, w2Max = 0.04, 0.1
    dw1=0.01
    dw2=0.001
    ww1 = w1Min:dw1:w1Max
    ww2 = w2Min:dw2:w2Max
    noiseRatios = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25]
    MCs = 1:100
    for nr in noiseRatios
        @info "nr = $nr"
        what = zeros(2, length(MCs))
        wnll = zeros(length(ww2), length(ww1), length(MCs))
        for (k,mc) in enumerate(MCs)
            _saveFile = joinpath(DATA_DIR,"goodwinCostFunData_nr.$(nr)_mc.$(mc).jld2")
            _retDic = JLD2.load(_saveFile)
            what[:,k] .= _retDic["what"]
            wnll[:,:,k] .= _retDic["wnll"]
        end
        saveFile = joinpath(DATA_DIR, "goodwinCostFunData.nr.$nr.summary.jld2")
        @info "  Saving to $saveFile"
        JLD2.save(
            saveFile, 
            Dict(
                "what"=>what,
                "wnll"=>wnll,
                "ww1"=>ww1,
                "ww2"=>ww2,
            )
        )
    end
end

###########################################################
################ Calling as a script ######################
###########################################################
if abspath(PROGRAM_FILE) == @__FILE__
    @info "Parsing Command Line Arguments"
    s = ArgParseSettings()
    @add_arg_table! s begin 
        "--noise-ratio"
            help = "noise ratios **for the data** (default 0.01)"
            arg_type = Float64
            default = 0.01
        "--monte-carlo"
            help = "number of sample at each parameterization"
            arg_type = Int
            default = 1
        "--summary"
            help = "Summarize data "
            arg_type = Int
            default = 1
    end
    args = parse_args(s)
    println("Parsed args:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end
    if args["summary"]
        _goodwinSummary()
    else 
        nr, mc = args["noise-ratio"],  args["monte-carlo"]
        _goodwinExperiment(nr, mc)
    end
end