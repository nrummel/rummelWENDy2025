## 
include("../src/util.jl")
include("plotBase.jl")
using Tullio

## Set the params
function plotEstimatorQuality(;
    ex = HINDMARSH_ROSE,
    corruptU0 = true,
)
    ## Unpack the saved data 
    fn = joinpath(DATA_DIR, "$(ex.name).summary.jld2")
    Mfull = length(ex.tt_full)-1
    @info """Setting up UQ plot for $(ex.name) 
        file $fn
    """
    odeName = String(split(splitpath(fn)[end], ".")[1])
    displayName = ex2Disp[odeName]
    resDic =  JLD2.load(fn)
    results = resDic["results"] 
    noiseRatios = copy(resDic["noiseRatios"])
    U0s = resDic["U0s"]
    subSampleRates = resDic["subSampleRates"]
    ## Process data
    @info "Preprocessing the data"
    dashopts = ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]
    wstar = ex.wTrue
    _I = length(noiseRatios) 
    III = length(subSampleRates)
    J = length(wstar)
    numRows = 5
    numCols = Int(floor(J / numRows))
    MC = size(results.what.wendy_mle,4)
    ii = findfirst(corruptU0 .== U0s) 
    what = Array{Union{Missing,Float64},4}(missing, _I, III, MC, J);
    cov = Array{Union{Missing,Float64},4}(missing, _I, III, MC, J) ;# diagonal of covariance
    failThresh = 25 
    isGood(x) = !isnan(x) && x < failThresh 
    goodIX = Matrix(undef, _I, III)
    for i in 1:_I, iii in 1:III
        this_yfsl2 = results.mean_fsl2.wendy_mle[i,ii,iii, :]
        this_ycl2 = results.cl2.wendy_mle[i,ii,iii, :]
        goodIX[i,iii] = findall(isGood.(this_yfsl2) .&& isGood.(this_ycl2))
        for mc in goodIX[i,iii]
            what[i,iii,mc,:] .= results.what.wendy_mle[i,ii,iii,mc]
            cov[i,iii,mc,:] .= try 
                diag(results.cov.wendy_mle[i,ii,iii,mc])
            catch 
                missing 
            end
            # Remove negative covariances?
            ixBad = findall(cov[i,iii,mc,:] .< 0)
            cov[i,iii,mc,ixBad] .= 0 
        end
    end

    what_bar = Array{Union{Missing,Float64}}(missing,_I,III,J) 
    cov_bar = similar(what_bar)
    for i in 1:_I, iii in 1:III, j in 1:J 
        what_bar[i,iii,j] = mean(skipmissing(what[i,iii,goodIX[i,iii],j]))
        cov_bar[i,iii,j] = mean(skipmissing(cov[i,iii,goodIX[i,iii],j]))
        if any( skipmissing(cov[i,iii,goodIX[i,iii],j]) .< 0)
            ixBad = findall(cov[i,iii,:,j] .< 0)
            @show cov[i,iii,ixBad,j]
            @assert false
        end
    end
    ## compute bias, variance, mse, and coverage
    _bias = similar(what)
    _variance = similar(what)
    _mse = similar(what)
    _coverage = similar(what)
    @tullio _bias[i,iii,mc,j]     = what[i,iii,mc,j] - wstar[j]
    @tullio _variance[i,iii,mc,j] = (what[i,iii,mc,j] - what_bar[i,iii,j])^2
    @tullio _mse[i,iii,mc,j]      = _bias[i,iii,mc,j]^2 + _variance[i,iii,mc,j]
    @tullio _coverage[i,iii,mc,j] = abs(_bias[i,iii,mc,j]) < (2*sqrt(cov[i,iii,mc,j]))
    # take the mean accross the third dimension
    bias     = similar(cov_bar)
    variance = similar(cov_bar)
    mse      = similar(cov_bar)
    coverage = similar(cov_bar)
    for i in 1:_I, iii in 1:III, j in 1:J 
        bias[i,iii,j] = mean(skipmissing(_bias[i,iii,goodIX[i,iii],j]))
        variance[i,iii,j] = mean(skipmissing(_variance[i,iii,goodIX[i,iii],j]))
        mse[i,iii,j] = mean(skipmissing(_mse[i,iii,goodIX[i,iii],j]))
        coverage[i,iii,j] = mean(skipmissing(_coverage[i,iii,goodIX[i,iii],j]))
    end
   
    ## Compute relative bias     
    relativeBias2 = similar(bias)
    relativeVar = similar(bias)
    relativeMSE = similar(bias)
    for j in 1:J 
        relativeBias2[:,:,j] .= bias[:,:,j].^2 ./ abs(wstar[j])^2 * 100
        relativeVar[:,:,j] .= variance[:,:,j] ./ abs(wstar[j])^2 * 100
        relativeMSE[:,:,j] .= mse[:,:,j] ./ abs(wstar[j])^2 * 100
    end
    coverage*=100 # to get coverage as a percentage

    open(joinpath(FIG_DIR, "$(ex.name)_UQ_summary.tex"), "w") do f
        s = """Across all parameters, the median estimated relative squared bias, variance, and mean squared error are $(@sprintf "%.3g" median(relativeBias2[:]))\\%, $(@sprintf "%.3g" median(relativeVar[:]))\\%, $(@sprintf "%.3g" median(relativeMSE[:]))\\% respectively. The median coverage across all parameters is $(@sprintf "%.3g" median(coverage[:]))\\%."""
        write(f, s)
    end

    p = make_subplots(rows=J, cols=2, shared_xaxes=true)
    dashopts = ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]
    for j in 1:J
        rowIx = j
        for (iii, sub) in enumerate(subSampleRates)
            dash  = dashopts[mod(iii-1, length(dashopts))+1]
            for (k,(arr, name)) in enumerate(zip([relativeBias2, relativeVar, relativeMSE], ["Bias²", "Variance", "MSE"]))
                add_trace!(
                    p, 
                    scatter(
                        x=noiseRatios,
                        y=arr[:,iii,j],
                        name="M = $(Int(round(1024/sub)))",
                        showlegend=j==1,
                        legendgrouptitle=attr(text=name),
                        legendgroup="$(name)_$j",  
                        marker_color=PLOTLYJS_COLORS[k],
                        mode="lines+markers",
                        line_dash=dash
                    ), 
                    row=j, 
                    col=1
                )
            end
            add_trace!(
                p, 
                scatter(
                    x=noiseRatios,
                    y=coverage[:,iii,j],
                    name="M = $(Int(round(1024/sub)))",
                    showlegend=j==1,
                    legendgrouptitle=attr(text="Coverage"),
                    legendgroup="Bias_$j",  
                    marker_color=PLOTLYJS_COLORS[4],
                    mode="lines+markers",
                    line_dash=dash
                ), 
                row=j, 
                col=2
            )
        end
        add_trace!(
                p, 
                scatter(
                    x=[0, noiseRatios[end]],
                    y=[95,95],
                    showlegend=false,
                    marker_color=CU_BOULDER_COLORS[3],
                    mode="lines",
                    line_dash="dash"
                ), 
                row=j, 
                col=2
            )
        ix=((j-1)*2)+1
        ix2=((j-1)*2)+2
        relayout!(
            p;
            (; zip(
                [
                    j == 1 ? :xaxis : Symbol("xaxis$ix"), 
                    j == 1 ? :yaxis : Symbol("yaxis$ix"),
                    Symbol("xaxis$ix2"),
                    Symbol("yaxis$ix2")
                ], 
                [
                    j > J - numCols ? attr(
                        tickangle=45,
                        tickmode="array",
                        tickvals=filter(nr->mod(Int(round(nr*100)),5)==0, noiseRatios),
                        ticktext=["$(Int(round(nr*100)))%" for nr in filter(nr->mod(Int(round(nr*100)),5)==0, noiseRatios)],
                        tickfont_size=15,
                        # range=[0,_i*III+1],
                        showgrid=false,
                    ) : attr(
                        tickmode="array",
                        tickvals=[],
                        ticktext=[],
                        # range=[0,_i*III+1]
                    ),
                    attr(
                        title_text="\$p_{$j}\$",
                        tickformat=".2g",
                        tickfont_size=15,
                        ticksuffix="%",
                        exponentformat="E"
                    ),
                    j > J - numCols ? attr(
                        tickangle=45,
                        tickmode="array",
                        tickfont_size=15,
                        tickvals=filter(nr->mod(Int(round(nr*100)),5)==0, noiseRatios),
                        ticktext=["$(Int(round(nr*100)))%" for nr in filter(nr->mod(Int(round(nr*100)),5)==0,noiseRatios)],
                        # range=[0,_i*III+1],
                        showgrid=false,
                    ) : attr(
                        tickmode="array",
                        tickvals=[],
                        ticktext=[],
                        # range=[0,_i*III+1]
                    ),
                    attr(
                        # title_text="\$p_{$j}\$",
                        tickformat=".0f",
                        tickfont_size=15,
                        ticksuffix="%",
                        exponentformat="E"
                    ),
                ]
            )...)...
        )
    end
    relayout!(
        p, 
        title=attr(
            text="$(displayName)",#<br>$(metric2Str[metric])",
            font_size=30,
            xanchor="center", 
            x=0.5,
            y=.985
        ),
        hovermode="x unified",
        legend=attr(
            # x=.925,
            y=0.5,
            yanchor="center",
            font=(
                family="sans-serif",
                size=15,
                color="#000"
            ),
            bgcolor="#E2E2E2",
            bordercolor= "#636363",
            borderwidth= 2,
        ),
    )
    return p
end

ex2HeightWidth = Dict(
    HINDMARSH_ROSE.name=>(800,600),
    GOODWIN.name=>(800,600),
    SIR.name=>(800,600),
)
if !isdir(FIG_DIR)
    @info "Save Director DNE $FIG_DIR"
    @info "  mkdir -p $FIG_DIR"
    mkpath(FIG_DIR)
end  
for ex in [HINDMARSH_ROSE, GOODWIN, SIR] 
    height, width = ex2HeightWidth[ex.name]
    p = plotEstimatorQuality(
        ex=ex, 
    ) 
    display(p)
    plot_file = joinpath(FIG_DIR, "$(ex.name)_estimateQuality.pdf")
    @info "Saving plot to $plot_file"
    savefig(p, plot_file, height=height, width=width)
end