##
includet("../src/util.jl")
includet("plotBase.jl")
##
w1Min, w1Max = 5, 20
w2Min, w2Max = 0.04, 0.1
dw1=0.01
dw2=0.001
ww1 = w1Min:dw1:w1Max
ww2 = w2Min:dw2:w2Max
wstar = [10, 0.0581]
noiseRatios = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25]
p1 =nothing
for nr in noiseRatios
    global p1
    @info "nr = $nr"
    savefile = joinpath(DATA_DIR, "goodwinCostFunData.nr.$nr.summary.jld2")
    @info "  Loading data from $savefile..."
    resDic = JLD2.load(savefile)
    wnll = resDic["wnll"]
    what = resDic["what"]
    wnll = mean(wnll, dims=3)[:,:,1]
    what = mean(what, dims=2)[:]
    min_wnll  = minimum(wnll[:])
    corrected_log(w) = log10(w - min_wnll + 1e-6 )
    p1 = plot(
        [
            contour(
                x=ww1, # horizontal axis
                y=ww2, # vertical axis
                z=corrected_log.(wnll),
                zmin=-7,
                zmax=5,
                colorbar_tickfont_size=25
            ),
            scatter(
                x=[what[1]],
                y=[what[2]],
                # text=["Iteration $i<br>z=$(@sprintf "%.2g" corrected_log(wnll_i))" for (i, wnll_i) in enumerate(wnll_its)],
                mode="markers",
                marker_size=20,
                marker_symbol="circle",
                opacity=0.8,
                marker_color="black",
            ) ,
            scatter( 
                x=[wstar[1]],
                y=[wstar[2]],
                # text="Truth<br>cost = $(@sprintf "%.2g" corrected_log(wnll_star))))",
                name="Truth",
                mode="markers",
                marker_symbol="star",
                marker_color="green",
                opacity=0.8,
                marker_size=20,
            )
        ],
        Layout(
            title=attr(
                text="Noise Ratio $(@sprintf "%.2g" nr*100)%", 
                font_size=40,
                x=0.5,
                xanchor="center",
            ),
            showlegend=false,
            xaxis=attr(
                title=attr(
                    text=L"\huge{p_4}",
                    font=attr(size=60)
                ),
                tickfont=attr(
                    size=25
                )
            ),
            yaxis=attr(
                title=attr(
                    text=L"\huge{p_5}",
                    font_size=60
                ),
                tickfont=attr(
                    size=25
                )
            )
        )
    )
    display(p1)
    savefig(
        p1, 
        joinpath(FIG_DIR, "goodwinExperiment_wnll_costSurface_nr.$nr.pdf"), width=550, height=500
    )
end
