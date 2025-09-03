ex = SIR
fn = joinpath(DATA_DIR, "$(ex.name).summary.jld2")
resDic =  JLD2.load(fn)
results = resDic["results"]

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
noiseRatios = specs[Symbol(ex.name)].noiseRatios
subSampleRates = specs[Symbol(ex.name)].subSampleRates
U0s = [true]
_I = length(noiseRatios)
III = 3
II = 1 
monteCarlos = 100
_fn = joinpath(DATA_DIR, "$(ex.name).summary.new.jld2")
metrics = [:what, :cov, :cl2, :mean_fsl2, :final_fsl2, :fsnll, :wnll, :dt, :iters ] 
algoNames = [:init,:truth,:oels,:wls,:wendy_irls,:wendy_mle, :hybrid]
_results = NamedTuple(
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
nameSwitch = (
    init=:init,
    truth=:truth,
    oels=:fslsq,
    wls=:wlsq,
    wendy_irls=:irwls,
    wendy_mle=:trustRegion,
    hybrid=:hybrid_trustRegion_fslsq,
)

for metric in metrics,algo in algoNames
    _results[metric][algo] .= results[metric][nameSwitch[algo]]
end
JLD2.save(_fn, 
  Dict(
        "results"=>_results,
        "metrics"=>metrics, 
        "algoNames"=>algoNames,
        "subSampleRates"=>subSampleRates,
        "noiseRatios"=>noiseRatios,
        "U0s"=>U0s,
    )
)

##

##
algoNames = Symbol[:fslsq,:trustRegion,:init,:truth,]
monteCarlos =30
newResults = [ NamedTuple(
    algoName=>AlgoRes()
    for algoName in algoNames
) for _ in endTimes, _ in noiseRatios, _ in 1:monteCarlos]
for (i,_) in enumerate(endTimes), (j,_) in enumerate(noiseRatios), (k,_) in enumerate(1:monteCarlos)
    for algo in algoNames
        newResults[i,j,k][algo].what = getproperty(results[i,j,k], algo).what
        newResults[i,j,k][algo].cov = getproperty(results[i,j,k], algo).cov
        newResults[i,j,k][algo].wits = getproperty(results[i,j,k], algo).wits
        newResults[i,j,k][algo].cl2 = getproperty(results[i,j,k], algo).cl2
        newResults[i,j,k][algo].mean_fsl2 = getproperty(results[i,j,k], algo).mean_fsl2
        newResults[i,j,k][algo].final_fsl2 = getproperty(results[i,j,k], algo).final_fsl2
        newResults[i,j,k][algo].fsnll = getproperty(results[i,j,k], algo).fsnll
        newResults[i,j,k][algo].wl2 = getproperty(results[i,j,k], algo).wl2
        newResults[i,j,k][algo].wnll = getproperty(results[i,j,k], algo).wnll
        newResults[i,j,k][algo].dt = getproperty(results[i,j,k], algo).dt
        newResults[i,j,k][algo].iters = getproperty(results[i,j,k], algo).iters
    end
end
JLD2.save(
    joinpath(DATA_DIR, "lorenzExperiment.new.jld2"),
    Dict(
        "noiseRatios"=>noiseRatios,
        "endTimes"=>endTimes,
        "results"=>newResults
    )
)