@info "Loading Dependencies"
includet("../src/util.jl")
includet("plotBase.jl")
## Load results from file 
@info "Loading results from file"
savefile    = joinpath(DATA_DIR, "lorenzExperiment.jld2")
noiseRatios = JLD2.load(savefile)["noiseRatios"]
endTimes    = JLD2.load(savefile)["endTimes"]
results     = JLD2.load(savefile)["results"]
nothing
##
@info "Plotting results"
metric         = :cl2
error_y        = false
skipZeroNoise  = true
N, numNRm, MCs = size(results)
trs = AbstractTrace[
    scatter(
        x=endTimes[3:end], 
        y=[mean(getproperty(res.init,metric) for res in results[n,1,:]) for n in 3:N],
        name="p₀ Solution",
        legendgroup="p₀ Solution",
        mode="lines",
        line=attr(
            color=CU_BOULDER_COLORS[3],
            dash="dash"
        ),
    )
]
dashopts = ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]
for (j, nr) in enumerate(noiseRatios) 
    if nr == 0 && skipZeroNoise
        continue 
    end
    global trs,tt, y_fslsq,y_wendy
    tt = []
    y_fslsq = []
    y_wendy = []
    for n in 3:N
        tn, yn_fs, yn_w = (
                endTimes[n],
                mean(min(getproperty(res.fslsq,metric),getproperty(res.init,metric)) for res in results[n,j,:]),
                mean(min(getproperty(res.trustRegion,metric),getproperty(res.init,metric)) for res in results[n,j,:])
            )
        push!(tt, tn)
        push!(y_fslsq, yn_fs)
        push!(y_wendy, yn_w)
    end
    trs = vcat(
        trs, 
        [   
            scatter(
                x=tt, 
                y=y_fslsq,
                error_y= attr(
                    type="data",
                    array=[std(getproperty(res.fslsq,metric) for res in results[n,j,:]) for n in 3:N],
                    visible=error_y
                ),
                name=skipZeroNoise ? algo2Disp[:oels] : "$(@sprintf "%.4g" nr*100)%",
                legendgroup="flsq",
                legendgrouptitle_text=skipZeroNoise ? nothing : algo2Disp[:oels],
                line_color=algo2Color[:oels],
                line_dash=skipZeroNoise ? "solid" : dashopts[mod(j-1, length(dashopts))+1]
            ),
            scatter(
                x=tt, 
                y=y_wendy,
                error_y= attr(
                    type="data",
                    array=[std(getproperty(res.trustRegion,metric) for res in results[n,j,:]) for n in 3:N],
                    visible=error_y
                ),
                name=skipZeroNoise ? algo2Disp[:wendy_mle] : "$(@sprintf "%.4g" nr*100)%",
                legendgroup="WENDy",
                legendgrouptitle_text=skipZeroNoise ? nothing : algo2Disp[:wendy_mle],
                line_color=algo2Color[:wendy_mle],
                line_dash=skipZeroNoise ? "solid" : dashopts[mod(j-1, length(dashopts))+1]
            ),
        ]
    )
end

yrange = metric == :cl2 ? [0, 0.16] : [0, 1]
p = plot(
    trs,
    Layout(
        margin_t=0,
        # title_text="Lorenz Oscillator",
        xaxis_title_text="Maximum Time Observed (s)",
        xaxis_title_font_size=25,
        xaxis_tickfont_size=20,
        yaxis_title_font_size=25,
        yaxis_tickfont_size=20,
        xaxis_range=(2.5,30.5),
        yaxis_title_text=metric == :cl2 ? "Relative Coefficient<br>Error" : "Relative Forward<br>Simulation Error",
        yaxis_range = yrange,
        # yaxis_type="log",
        hovermode="x unified",
        legend=attr(
            # x=.925,
            y=0.5,
            yanchor="center",
            font=(
                family="sans-serif",
                size=20,
                color="#000"
            ),
            bgcolor="#E2E2E2",
            bordercolor= CU_BOULDER_COLORS[3],
            borderwidth= 2,
        ),
    )
)
display(p)
##
if !isdir(FIG_DIR)
    mkpath(FIG_DIR)
end   
plot_file = joinpath(FIG_DIR, "lorenzExperiment.pdf")
@info "Saving to $plot_file"
savefig(p,plot_file, height=300, width=900 )
