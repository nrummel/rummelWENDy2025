## 
include("../src/util.jl")
include("plotBase.jl")
using Tullio
## Set the params
function plotConfidenceIntervals(;
    ex = HINDMARSH_ROSE,
    numRows = 5,
    corruptU0 = true,
    maxNoise = 0.10
)
    ## Unpack the saved data 
    fn = joinpath(DATA_DIR, "$(ex.name).summary.jld2")
    @info """Setting up UQ plot for $(ex.name) 
        file $fn
    """
    odeName = String(split(splitpath(fn)[end], ".")[1])
    resDic =  JLD2.load(fn)
    results = resDic["results"] 
    noiseRatios = copy(resDic["noiseRatios"])
    U0s = resDic["U0s"]
    subSampleRates = resDic["subSampleRates"]
    ## Process data
    @info "Preprocessing the data"
    wstar = ex.wTrue
    _I = length(noiseRatios) 
    III = length(subSampleRates)
    J = length(wstar)
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
    @info "Plotting UQ"
    numCols = Int(floor(J / numRows))
    p = make_subplots(rows=numRows, cols=numCols, shared_xaxes=true)
    _nrs = filter(x->x <= maxNoise, noiseRatios)
    _i = length(_nrs) 
    for j in 1:J
        wj = []
        cj = []
        label = []
        wj = mean(what_bar[1:_i,:,j], dims=2)[:]
        cj = mean(cov_bar[1:_i,:,j],dims=2)[:]
        label = ["$(@sprintf "%.2g" nr*100)%" for nr in _nrs]
        _xx = 1:_i
        rowIx = mod(j-1, numRows) + 1
        colIx = Int(floor((j-1)/numRows)) + 1
        add_trace!(
            p, 
            scatter(
                x=_xx,
                y=wj,
                error_y=attr(
                    type="data",
                    array=2.0.*sqrt.(cj),
                    visible=true,
                    thickness=5
                ),
                name=LaTeXString("p_$(j)"),
                showlegend=false, 
                marker_color=CU_BOULDER_COLORS[2],
                mode="markers"
            ), 
            row=rowIx, 
            col=colIx
        )
        add_trace!(
            p, 
            scatter(
                x=[0, _i*2],
                y=ex.wTrue[j]*ones(2),
                name=LaTeXString("w^*_$(j)"),
                showlegend=false, 
                mode="lines",
                line_dash="dash",
                line_width=5,
                line_color=CU_BOULDER_COLORS[4]
            ), 
            row=rowIx, 
            col=colIx
        )
        relayout!(
            p;
            (; zip(
                [
                    j == 1 ? :xaxis : Symbol("xaxis$j"), 
                    j == 1 ? :yaxis : Symbol("yaxis$j")
                ], 
                [
                    j > J - numCols ? attr(
                        tickangle=45,
                        tickmode="array",
                        tickvals=_xx,
                        ticktext=label,
                        range=[0,_i+1],
                        showgrid=false,
                        tickfont_size=25
                    ) : attr(
                        tickmode="array",
                        tickvals=[],
                        ticktext=[],
                        range=[0,_i+1],
                        tickfont_size=25
                    ),
                    attr(
                        title_text="\$p_{$j}\$",
                        title_font_size=25,
                        # showexponent="all",
                        exponentformat="E",
                    )
                ]
            )...)...
        )
    end
    relayout!(p,
         title=attr(
            text="$(ex2Disp[odeName])",#<br>$(metric2Str[metric])",
            font_size=30,
            xanchor="center", 
            yanchor="bottom", 
            x=0.5,
            y=0.92
        ),
    )
    return p
end

ex2maxNoise = Dict(
    HINDMARSH_ROSE.name=>0.1,
    GOODWIN.name=>0.1,
    SIR.name=>0.1
)

ex2NumRows = Dict(
    HINDMARSH_ROSE.name=>5,
    GOODWIN.name=>4,
    SIR.name=>5,
)
ex2HeightWidth = Dict(
    HINDMARSH_ROSE.name=>(600,1000),
    GOODWIN.name=>(600,1000),
    SIR.name=>(600,1000),
)
if !isdir(FIG_DIR)
    @info "Save Director DNE $FIG_DIR"
    @info "  mkdir -p $FIG_DIR"
    mkpath(FIG_DIR)
end  
for ex in [HINDMARSH_ROSE, GOODWIN, SIR] 
    p = plotConfidenceIntervals(
        ex=ex, 
        numRows=ex2NumRows[ex.name], 
        maxNoise=ex2maxNoise[ex.name]
    ) 
    display(p)
    plot_file = joinpath(FIG_DIR, "$(ex.name)_confidenceIntervals.pdf")
    @info "Saving plot to $plot_file"
    height, width = ex2HeightWidth[ex.name]
    savefig(p, plot_file, height=height, width=width)
end
