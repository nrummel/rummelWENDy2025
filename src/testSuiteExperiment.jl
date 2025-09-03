includet("util.jl")
using ProgressMeter: @showprogress
#
specs = (
    hindmarshRose = (
        noiseRatios = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2],
        subSampleRates = [1,2,4],
        ex=HINDMARSH_ROSE,
        J = length(HINDMARSH_WSTAR)
    ),
   goodwin = (
        noiseRatios =[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25],
        subSampleRates = [1,2,4],
        ex=GOODWIN,
        J = length(GOODWIN_WSTAR)
    ),
    sir = (
        noiseRatios =[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2],
        subSampleRates = [1,2,4],
        ex=SIR,
        J = length(SIR_WSTAR)
    )
)

exName = :hindmarshRose 
noiseRatios = specs[exName].noiseRatios
U0s = [true]
subSampleRates = specs[exName].subSampleRates
ex = specs[exName].ex
J = specs[exName].J
monteCarlos = 100
metrics = [:what, :cov, :cl2, :mean_fsl2, :final_fsl2, :fsnll, :wnll, :dt, :iters ] 
algoNames = [:init,:truth,:oels,:wls,:wendy_irls,:wendy_mle, :hybrid]
## Preallocate memory
_I, II, III = length(noiseRatios), length(U0s), length(subSampleRates)
results = NamedTuple(
    metric=>NamedTuple(
        algo=>if metric == :cov || metric == :what || metric == :w0
            Array{Any,4}(missing, _I, II, III, monteCarlos)
        else 
            NaN*ones(_I, II, III, monteCarlos)
        end
        for algo in algoNames
    )
    for metric in metrics
)
@showprogress for (i, nr) in enumerate(noiseRatios), (ii, corruptU0) in enumerate(U0s), (iii, sub) in enumerate(subSampleRates), mc in 1:monteCarlos
    # Specify the simulation params
    simParams = SimulationParameters(
        noiseRatio=nr,
        seed=mc,
        timeSubsampleRate=sub, 
        corruptU0=corruptU0
    )
    # Run all algorithms on this instance of the noise and data
    global res
    _, _, res = runExample(ex, simParams)
    for algo in algoNames, metric in metrics
        if metric == :what                 
            what_mc = getproperty(res[algo], metric)[1:J]
            results[metric][algo][i,ii,iii,mc] = what_mc
        else 
            results[metric][algo][i,ii,iii,mc] = getproperty(res[algo], metric)
        end
    end
end
saveFile = joinpath(DATA_DIR, "$(ex.name).summary.new.jld2")
@info "Saving to $saveFile"
JLD2.save( 
    saveFile,
    Dict(
        "results"=>results,
        "metrics"=>metrics, 
        "algoNames"=>algoNames,
        "subSampleRates"=>subSampleRates,
        "noiseRatios"=>noiseRatios,
        "U0s"=>U0s,
    )
)

##
saveFile = joinpath(DATA_DIR, "$(ex.name).summary.new.jld2")
