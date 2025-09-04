includet("../src/util.jl")
includet("plotBase.jl")
## 
function plotSweep(;
    ex = HINDMARSH_ROSE,
    nominalVal=NaN,
    center_fun = mean,
    yaxis_type = "log", # "linear"
    yaxis_range = nothing, # (0,100),
    maxNoise = 0.3,
    algoNames = [:init,:truth,:oels,:wls,:wendy_irls,:wendy_mle, :hybrid],
    avgAcrossSubs=true,
)
    fn = joinpath(DATA_DIR, "$(ex.name).summary.jld2")
    Mfull = length(ex.tt_full)-1
    @info """Setting up summary plot for $(ex.name) 
        file $fn
    """
    odeName = String(split(splitpath(fn)[end], ".")[1])
    displayName = ex2Disp[odeName]
    resDic =  JLD2.load(fn)
    results = resDic["results"] # zeros(nr,ic,sub,mc,J)
    maxNoise
    noiseRatios = resDic["noiseRatios"]
    filter!(x->x <= maxNoise, noiseRatios)
    noiseRatios
    maxNoise = minimum([maxNoise, maximum(noiseRatios)])
    U0s = resDic["U0s"]
    subSampleRates = resDic["subSampleRates"]
    @info "Plotting $displayName "
    ii = 1
    III = length(subSampleRates)
    _I = length(noiseRatios) # size(results[:cl2][:init], 1) # number of noiseRatios
    yy = NamedTuple(
        algo=>(cl2=zeros(_I,III), cl2_std=zeros(_I,III), fsl2=zeros(_I,III), fsl2_std=zeros(_I,III), failRate=zeros(_I,III)) 
    for algo in algoNames)
    for (_, algo) in enumerate(algoNames), (iii, _) in enumerate(subSampleRates)
        # algo, metric, sub, 1:length(noiseRatios), ii, iii
        Yfsl2 = results[:mean_fsl2][algo][1:length(noiseRatios),ii,iii,:]
        Ycl2 = results[:cl2][algo][1:length(noiseRatios),ii,iii,:]
        failThresh = 25
        isFail(x) = isnan(x) || x >= failThresh
        isGood(x) = !isnan(x) && x < failThresh 
        for  (nrIx,_) in enumerate(noiseRatios)
            this_yfsl2 = Yfsl2[nrIx, :]
            this_ycl2 = Ycl2[nrIx, :]
            MCs = length(Yfsl2[nrIx, :])
            failIX = findall(isFail.(this_yfsl2) .|| isFail.(this_ycl2))
            goodIX = findall(isGood.(this_yfsl2) .&& isGood.(this_ycl2))
            @assert !isempty(goodIX)
            yy[algo].cl2[nrIx,iii] = center_fun(this_ycl2[goodIX])
            yy[algo].cl2_std[nrIx,iii] = std(this_ycl2[goodIX])
            yy[algo].fsl2[nrIx,iii] = center_fun(this_yfsl2[goodIX])
            yy[algo].fsl2_std[nrIx,iii] = std(this_yfsl2[goodIX])
            yy[algo].failRate[nrIx,iii] = length(failIX) / MCs
        end
    end 
    
    ## Plotting
    p = make_subplots(rows=2, cols=2)
    linespec = ["solid", "dash", "dot"]
    for (k,algo) in enumerate(algoNames), (iii, sub) in enumerate(subSampleRates)
        if algo == :truth 
            continue 
        end
        if avgAcrossSubs 
            if iii > 1
                continue 
            end
            cl2 = mean(yy[algo].cl2, dims=2)[:]

            cl2_std =  mean(yy[algo].cl2_std, dims=2)[:]
            fsl2 = mean(yy[algo].fsl2, dims=2)[:]
            fsl2_std =  mean(yy[algo].fsl2_std, dims=2)[:]
            failRate = mean(yy[algo].failRate, dims=2)[:]
            mean_failRate = mean(yy[algo].failRate[:])*100
            name= algo2Disp[algo]
            legendgrouptitle=nothing
            dash =nothing
        else 
            cl2 = yy[algo].cl2[:,iii]
            cl2_std = yy[algo].cl2_std[:,iii]
            fsl2 = yy[algo].fsl2[:,iii]
            fsl2_std = yy[algo].fsl2_std[:,iii]
            failRate = yy[algo].failRate[:,iii]
            mean_failRate = mean(yy[algo].failRate[:,iii])*100
            name=algo == :init ? "p₀" : "M = $(@sprintf "%d" (Mfull/sub))"
            legendgrouptitle=attr(
                legendgroup=algo,
                text= (algo ==:init ) ? nothing : algo2Disp[algo],
                font_size=15
            )
            dash = linespec[iii]
        end
        if algo == :init && iii > 1 
            continue 
        end
        add_trace!(
            p,
            scatter(
                x=noiseRatios,
                y=cl2,
                text=[
                    "$(@sprintf "%.3g" yi) +/- $(@sprintf "%.3g" stdyi), failRate = $(@sprintf "%.3g" (failRatei*100))%" 
                for (yi, stdyi, failRatei) in zip(cl2, cl2_std, failRate)],
                hoverinfo="text",
                mode="lines",
                line=attr(
                    color= algo2Color[algo],
                    dash=dash,
                    width=4,
                ),
                name=name,
                legendgrouptitle=legendgrouptitle,
                legendgroup=algo,
                showlegend=true,
            ),
            row=1, 
            col=1
        )
        add_trace!(
            p,
            scatter(
                x=noiseRatios,
                y=fsl2,
                text=[
                    "$(@sprintf "%.3g" yi) +/- $(@sprintf "%.3g" stdyi), failRate = $(@sprintf "%.3g" (failRatei*100))%" 
                for (yi, stdyi, failRatei) in zip(fsl2, fsl2_std, failRate)],
                hoverinfo="text",
                mode="lines",
                line=attr(
                    color= algo2Color[algo],
                    dash=dash,
                    width=4,
                ),
                name=algo == :init ? "p₀" : algo2Disp[algo],
                legendgroup=algo,
                showlegend=false,
            ),
            row=1, 
            col=2
        )
        if algo == :truth || algo == :init
            continue
        end
        add_trace!(
            p,
            bar(
                x=[iii],
                y=[mean_failRate+10],
                text=["$(@sprintf "%.0f" mean_failRate[1])%"],
                textposition="outside",
                marker_color=algo2Color[algo],
                legendgroup=algo,
                showlegend=false,
            ),
            row=2, 
            col=1
        )
    end
    ## set y tick vals
    if yaxis_type == "log" && isnothing(yaxis_range)
        # Obtain minimum and maximum values for 
        lcl2,ucl2 = (Inf,-Inf)
        lfsl2,ufsl2 = (Inf,-Inf)
        for algo in algoNames
            if algo == :truth 
                continue
            end
            _ycl2 = filter(!isnan, yy[algo].cl2[:])
            lcl2 = min(lcl2,floor(minimum(log10.(_ycl2))))
            ucl2 = max(ucl2,ceil(maximum(log10.(_ycl2))))
            _yfsl2 = filter(!isnan,yy[algo].fsl2[:])
            lfsl2 = min(lfsl2,floor(minimum(log10.(_yfsl2))))
            ufsl2 = max(ufsl2,ceil(maximum(log10.(_yfsl2))))
        end
        ucl2 += 0.5
        ufsl2 += 0.5
        cl2_yaxis_tickvals = [10.0^e for e in lcl2:ucl2]
        fsl2_yaxis_tickvals = [10.0^e for e in lfsl2:ufsl2]
        cl2_yaxis_ticktext = []
        fsl2_yaxis_ticktext = []
        for e in lcl2:ucl2 
            pp = 10.0^e*100 
            if pp >= 1
               pp = Int(round(10.0^e*100))
            end
            push!(cl2_yaxis_ticktext, "$(pp)%")
        end
        for e in lfsl2:ufsl2 
            pp = 10.0^e*100 
            if pp >= 1
               pp = Int(round(10.0^e*100))
            end
            push!(fsl2_yaxis_ticktext, "$(pp)%")
        end
        cl2_yaxis_range = (lcl2,ucl2)
        fsl2_yaxis_range = (lfsl2,ufsl2)
    else
        cl2_yaxis_tickvals = nothing
        fsl2_yaxis_tickvals =  nothing
        cl2_yaxis_range = nothing
        fsl2_yaxis_range =  nothing
    end
    realAlgos = filter(algo-> algo != :truth && algo != :init, algoNames)
    maxFailPercent = avgAcrossSubs ? maximum(mean(yy[algo].failRate[:]) for algo in realAlgos) : maximum(maximum(mean(yy[algo].failRate, dims=1)[:]) for algo in realAlgos)
    maxFailPercent = Int(round(maxFailPercent*100 + 5))
    while mod(maxFailPercent,5) != 0
        maxFailPercent +=1
    end
    failure_ytickvals = 
    if maxFailPercent > 50
        range(0, maxFailPercent+10, step=10)
    elseif maxFailPercent > 20
        range(0, maxFailPercent+10, step=5)
    else 
        0:maxFailPercent+10
    end
    ## Layout plot 
    relayout!(p,
        margin_t=0,
        showlegend=true,
        xaxis=attr(
            title_text= avgAcrossSubs ? "Noise Ratio" : "",
            title_font_size=15,
            tickfont_size=15,
            tickvals=0.0:0.05:maxNoise,
            ticktext=["$(Int(round(nr*100)))%" for nr in 0.0:.05:maxNoise],
            tickangle=45,
            domain=avgAcrossSubs ? [0, 0.45] : [0, 0.8],
        ),
        yaxis=attr( 
            type=yaxis_type,
            range=cl2_yaxis_range,
            tickvals=cl2_yaxis_tickvals,
            ticktext=cl2_yaxis_ticktext,
            title_text="Coefficient Error",
            title_font_size=15,
            tickfont_size=15,
            domain= avgAcrossSubs ? [0.55, 0.9] : [0.625, 0.9],
        ),
        xaxis2=attr(
            title_text= avgAcrossSubs ? "Noise Ratio" : "",
            title_font_size=15,
            tickfont_size=15,
            tickvals=0:0.05:maxNoise,
            ticktext=["$(Int(round(nr*100)))%" for nr in 0:.05:maxNoise],
            tickangle=45,
            domain=avgAcrossSubs ? [0.55, 1] : [0, 0.8],
        ),
        yaxis2=attr( 
            type=yaxis_type,
            range=fsl2_yaxis_range,
            side = avgAcrossSubs ? "right" : "left",
            tickvals=fsl2_yaxis_tickvals,
            ticktext=fsl2_yaxis_ticktext,
            title_text="Forward Simulation Error",
            title_font_size=15,
            tickfont_size=15,
            domain=avgAcrossSubs ? [0.55, 0.9] : [0.325, 0.575],
        ),
        xaxis3=attr(
            title_text="",
            tickvals=[1],
            ticktext=[""],
            title_font_size=15,
            tickfont_size=15,
            domain=avgAcrossSubs ? [0, 0.45] : [0, 0.8],
        ),
        yaxis3=attr(
            range=(5,maxFailPercent+10),
            tickvals=failure_ytickvals,
            ticktext=[p < 10 ? "" : "$(p-10)%" for p in failure_ytickvals],
            title_text="Failure Rate",
            title_font_size=15,
            tickfont_size=15,
            domain=avgAcrossSubs ? [0, 0.35] : [0, 0.275]
        ),
        uniformtext_minsize=12,
        uniformtext_mode="show",
        title=attr(
            text="$(displayName)",
            font_size=20,
            xanchor="center", 
            x=0.5,
            y=0.95
        ),
        legend=attr(
            x= avgAcrossSubs ? 0.775 : 0.99,
            y= avgAcrossSubs ? 0 : 0.5,
            xanchor="center",
            yanchor="center",
            font=(
                family="sans-serif",
                size=15,
                color="#000"
            ),
            bgcolor="#E2E2E2",
            bordercolor= "#636363",
            borderwidth= 2,
            title=attr(
                text="Legend",
                font_size=15,
                side="top center", 
            )
        ),
        hovermode="x unified",
    )
    p
end
##
ex2maxNoise = Dict(
    HINDMARSH_ROSE.name=>0.2,
    GOODWIN.name=>0.25,
    SIR.name=>0.25
)
save_dir = FIG_DIR
if !isdir(save_dir)
    mkpath(save_dir)
end 
for ex in [HINDMARSH_ROSE, GOODWIN, SIR]
    p = plotSweep(ex=ex,  maxNoise= ex2maxNoise[ex.name], avgAcrossSubs=true) 
    display(p)
    plot_file = joinpath(FIG_DIR, "$(ex.name)_accuracy_avg.pdf")
    @info "Saving to $plot_file"
    savefig(p,plot_file, height=600, width=600)
    p = plotSweep(ex=ex,  maxNoise= ex2maxNoise[ex.name], avgAcrossSubs=false) 
    display(p)
    plot_file = joinpath(FIG_DIR, "$(ex.name)_accuracy.pdf")
    @info "Saving to $plot_file"
    savefig(p,plot_file, height=1000, width=600)
end
nothing
