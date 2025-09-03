## Load everything we need 
@info "Loading Dependencies"
include("util.jl")
## 
@info "Setting up Lorenz Experiment"
# Build a data set that we will truncate in time rather than subsampling
noiseRatios = [0.1] # [0, 0.01, 0.05, 0.1, 0.2] 
timeSubsampleRate, monteCarlos = 1, 30
dt = 0.01
fullEndTime = 30.0
endTimes = vcat(collect(1:1:10), collect(12:2:14), collect(15:5:fullEndTime))
# These end time caused errors which I believe is due to some numerical issue
endTimes[findmin(endTimes .- 6 )[2]] += 0.1
endTimes[findmin(endTimes .- 8 )[2]] += 0.1
# 
params = WENDyParameters(;
    optimMaxiters  = 200,
    optimTimelimit = 200.0,
    radiusMinTime  = 2*dt,
    radiusMaxTime  = 5.,
    Kᵣ             = 200,
    Kmax           = 500,
)
ex = SimulatedWENDyData(
    "lorenz", 
    LORENZ_f!,
    LORENZ_INIT_COND,
    (0.0,fullEndTime),
    LORENZ_WSTAR,
    [    
        (5.0,15.0),
        (23.0,33.0),
        (2.0,3.0)
    ],
    params;
    dt=dt, 
    Mp1=nothing,
    linearInParameters=Val(false),
    noiseDist=Val(Normal),
    forceOdeSolve=true,
    ll=Warn
);
algoNames = Symbol[:init,:truth,:fslsq,:trustRegion]
# preallocate
results = [ NamedTuple(
    algoName=>AlgoRes()
    for algoName in algoNames
) for _ in endTimes, _ in noiseRatios, _ in 1:monteCarlos]
# wendyProbs = Array{WENDyProblem,3}(undef, length(endTimes), length(noiseRatios), monteCarlos)
# loop through everything
@info "Starting Loop"
pbar = Progress(length(endTimes)*length(noiseRatios)*monteCarlos)
J = length(LORENZ_WSTAR)
for mc in 1:monteCarlos
    Random.seed!(mc)
    w0 = zeros(J)
    for _j in 1:J 
        d = Uniform(ex.wRng[_j]...)
        w0[_j] = rand(d)
    end
    # loop through noise ratios
    for (j, nr) in enumerate(noiseRatios)
        # simulate the data so that everyone has the same noise
        simParams = SimulationParameters(
            seed=mc, 
            timeSubsampleRate=timeSubsampleRate,
            noiseRatio=nr,
            corruptU0=true
        ) 
        simulate!(ex, simParams, ll=Warn)
        # extract data from their structures
        tt, U, Ustar, f! = ex.tt[], ex.U[], ex.Ustar[], ex.f!
        # Run all the algorithms and collect results s
        for (i,T) in enumerate(endTimes)
            # @info "========================================================="
            # @info "T=$(@sprintf "%.4g" T)"
            _Mp1 = findlast(tt .<= T) 
            @assert !isnothing(_Mp1) "We are computing this index improperly"
            wendyProb = WENDyProblem(tt[1:_Mp1], U[1:_Mp1,:], f!, J, Val(false), Val(Normal), ex.params, ll=Warn)
            # Loop through algos and run this example
            for algoName in algoNames 
                # @info "  $algoName"
                alg_dt = @elapsed begin   
                    what, iters, wits = if algoName == :truth 
                        ex.wTrue, 0, reshape(ex.wTrue,J,1)
                    elseif algoName == :init
                        w0, 0, reshape(w0,J,1)
                    else
                        try 
                            solve(wendyProb, w0, ex.params; alg=algoName, return_wits=true)
                        catch
                            (NaN * ones(J), -1, zeros(J,0)) 
                        end
                    end
                end
                @views computeMetrics!(
                    results[i,j,mc][algoName], wendyProb, 
                    f!, tt[1:_Mp1], Ustar[1:_Mp1,:], U[1:_Mp1,:], ex.sigTrue[], ex.wTrue, ex.params, 
                    alg_dt, what, iters, wits; ll=Warn
                ) 
            end
            next!(
                pbar, 
                desc="nr=$(@sprintf "%.4g" nr*100)%, T=$(@sprintf "%.4g" T), mc = $mc / $monteCarlos",
                showvalues=vcat([
                    ("$algo $metric", "$(@sprintf "%.4g" getproperty(results[i,j,mc][algo], metric))")
                for metric in [:dt, :cl2, :mean_fsl2], algo in [:init, :fslsq, :trustRegion]][:], [("K", wendyProb.K)])    
            )
        end
    end
end
nothing
##
savefile = joinpath(DATA_DIR, "lorenzExperiment.jld2")
@info "Saving to $savefile"
JLD2.save(
    savefile,
    Dict(
        "noiseRatios"=>noiseRatios,
        "endTimes"=>endTimes,
        "results"=>results
    )
)