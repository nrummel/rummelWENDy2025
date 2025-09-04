includet("plotBase.jl")
function fit_line(Ms, times, alpha=1)
    N = length(times)
    ix =3
    log_mstar = log10(Ms[ix])
    log_tstar = log10(times[ix])
    _f(m) = 10^(alpha*(log10.(m) - log_mstar) + log_tstar)
    alpha, _f.(Ms) 
end 
##
@info "M Experiment "
savefile = joinpath(DATA_DIR, "computationalCostExperiment_M.jld2")
resDic = JLD2.load(savefile)
Ms,l_wnll_times,l_grad_times,l_hess_times,nl_wnll_times,nl_grad_times,nl_hess_times = (resDic["Ms"],resDic["l_wnll_times"],resDic["l_grad_times"],resDic["l_hess_times"],resDic["nl_wnll_times"],resDic["nl_grad_times"],resDic["nl_hess_times"])
α_l_wnll, y_l_wnll = fit_line(Ms, l_wnll_times)
α_l_grad, y_l_grad = fit_line(Ms, l_grad_times)
α_l_hess, y_l_hess = fit_line(Ms, l_hess_times)
α_nl_wnll, y_nl_wnll = fit_line(Ms, nl_wnll_times)
α_nl_grad, y_nl_grad = fit_line(Ms, nl_grad_times)
α_nl_hess, y_nl_hess = fit_line(Ms, nl_hess_times)
p = plot(
    [
        scatter(
            x=Ms, 
            y=l_hess_times, 
            mode="markers", 
            marker_symbol="diamond", 
            marker_color=PLOTLYJS_COLORS[3], 
            legendgroup="hess", 
            name="Linear", 
            legendgrouptitle=attr(
                text=L"\nabla^2 \ell",
                font_size=20
            )
        ),
        scatter(
            x=Ms, 
            y=y_l_hess, 
            marker_color=PLOTLYJS_COLORS[3], 
            legendgroup="hess", 
            name="",showlegend=false, 
            line_dash="dash", 
            mode="lines"
        ),
        scatter(
            x=Ms, 
            y=nl_hess_times, 
            mode="markers", marker_symbol="circle", 
            marker_color=PLOTLYJS_COLORS[3], 
            legendgroup="hess", 
            name="Nonlinear"
        ),
        scatter(
            x=Ms, 
            y=y_nl_hess,
            marker_color=PLOTLYJS_COLORS[3], 
            legendgroup="hess", 
            name="",showlegend=false, 
            line_dash="dot", 
            mode="lines"
        ),
        scatter(
            x=Ms, 
            y=l_grad_times, 
            mode="markers", 
            marker_symbol="diamond", 
            marker_color=PLOTLYJS_COLORS[2], 
            legendgroup="grad", 
            name="Linear", 
            legendgrouptitle=attr(
                text=L"\nabla \ell",
                font_size=20
            )
        ),
        scatter(
            x=Ms, 
            y=y_l_grad, 
            marker_color=PLOTLYJS_COLORS[2], 
            legendgroup="grad", 
            name="",showlegend=false, 
            line_dash="dash", 
            mode="lines"
        ),
        scatter(
            x=Ms, 
            y=l_wnll_times, 
            mode="markers", marker_symbol="diamond", marker_color=PLOTLYJS_COLORS[1], 
            legendgroup="wnll", 
            name="Linear", 
            legendgrouptitle=attr(
                text=L"\ell",
                font_size=20
            )
        ),
        scatter(
            x=Ms, 
            y=y_l_wnll, 
            marker_color=PLOTLYJS_COLORS[1], 
            legendgroup="wnll", 
            name="",showlegend=false,  
            line_dash="dash", 
            mode="lines"
        ),
        scatter(
            x=Ms, 
            y=nl_wnll_times, 
            mode="markers", 
            marker_symbol="circle", 
            marker_color=PLOTLYJS_COLORS[1], 
            legendgroup="wnll", 
            name="Nonlinear"
        ),
        scatter(
            x=Ms, 
            y=y_nl_wnll,
            marker_color=PLOTLYJS_COLORS[1], 
            legendgroup="wnll", 
            name="",showlegend=false, 
            line_dash="dot", 
            mode="lines"
        ),
        scatter(
            x=Ms, 
            y=nl_grad_times, 
            mode="markers", 
            marker_symbol="circle", 
            marker_color=PLOTLYJS_COLORS[2], 
            legendgroup="grad", 
            name="Nonlinear"
        ),
        scatter(
            x=Ms, 
            y=y_nl_grad,
            marker_color=PLOTLYJS_COLORS[2], 
            legendgroup="grad", 
            name="",showlegend=false,
            line_dash="dot", mode="lines"
        ),
    ],
    Layout( 
        # title="Computational Cost Comparison<br>Hindmarsh Rose",
        margin_t=0,
        xaxis=attr(title_text="M", type="log",title_font_size=20,tickfont_size=20),
        yaxis=attr(title_text="time (s)", type="log", title_font_size=20,tickfont_size=20),
        legend=attr(
            # x=.925,
            y=0.5,
            yanchor="center",
            font=(
                family="sans-serif",
                size=18,
                color="#000"
            ),
            bgcolor="#E2E2E2",
            bordercolor= "#636363",
            borderwidth= 2,
        ),
    )
)
display(p)
savefig(p,joinpath(FIG_DIR, "computationalCost_Mp1.pdf"), width=550, height=400)
##
@info "K Experiment "
savefile = joinpath(DATA_DIR, "computationalCostExperiment_K.jld2")
resDic = JLD2.load(savefile)
Ks,l_wnll_times,l_grad_times,l_hess_times,nl_wnll_times,nl_grad_times,nl_hess_times = (resDic["Ks"],resDic["l_wnll_times"],resDic["l_grad_times"],resDic["l_hess_times"],resDic["nl_wnll_times"],resDic["nl_grad_times"],resDic["nl_hess_times"])
α_l_wnll, y_l_wnll = fit_line(Ks, l_wnll_times, 3)
α_l_grad, y_l_grad = fit_line(Ks, l_grad_times, 3)
α_l_hess, y_l_hess = fit_line(Ks, l_hess_times, 3)
α_nl_wnll, y_nl_wnll = fit_line(Ks, nl_wnll_times, 3)
α_nl_grad, y_nl_grad = fit_line(Ks, nl_grad_times, 3)
α_nl_hess, y_nl_hess = fit_line(Ks, nl_hess_times, 3)
p = plot(
    [
        scatter(x=Ks, y=l_hess_times, mode="markers", marker_symbol="diamond", marker_color=PLOTLYJS_COLORS[3], legendgroup="hess", name="Linear", legendgrouptitle=attr(
            text=L"\nabla^2 \ell",
            font_size=20)
        ),
        scatter(x=Ks, y=y_l_hess, marker_color=PLOTLYJS_COLORS[3], legendgroup="hess", name="",showlegend=false, line_dash="dash", mode="lines"),
        scatter(x=Ks, y=nl_hess_times, mode="markers", marker_symbol="circle", marker_color=PLOTLYJS_COLORS[3], legendgroup="hess", name="Nonlinear"),
        scatter(x=Ks, y=y_nl_hess,marker_color=PLOTLYJS_COLORS[3], legendgroup="hess", name="",showlegend=false, line_dash="dot", mode="lines"),
        scatter(x=Ks, y=l_grad_times, mode="markers", marker_symbol="diamond", marker_color=PLOTLYJS_COLORS[2], legendgroup="grad", name="Linear", legendgrouptitle=attr(
            text=L"\nabla \ell",
            font_size=20)
        ),
        scatter(x=Ks, y=y_l_grad, marker_color=PLOTLYJS_COLORS[2], legendgroup="grad", name="",showlegend=false, line_dash="dash", mode="lines"),
        scatter(x=Ks, y=nl_grad_times, mode="markers", marker_symbol="circle", marker_color=PLOTLYJS_COLORS[2], legendgroup="grad", name="Nonlinear"),
        scatter(x=Ks, y=y_nl_grad,marker_color=PLOTLYJS_COLORS[2], legendgroup="grad", name="",showlegend=false, line_dash="dot", mode="lines"),
        scatter(x=Ks, y=l_wnll_times, mode="markers", marker_symbol="diamond", marker_color=PLOTLYJS_COLORS[1], legendgroup="wnll", name="Linear", legendgrouptitle=attr(
            text=L"\ell",
            font_size=20)
        ),
        scatter(x=Ks, y=y_l_wnll, marker_color=PLOTLYJS_COLORS[1], legendgroup="wnll", name="",showlegend=false,  line_dash="dash", mode="lines"),
        scatter(x=Ks, y=nl_wnll_times, mode="markers", marker_symbol="circle", marker_color=PLOTLYJS_COLORS[1], legendgroup="wnll", name="Nonlinear"),
        scatter(x=Ks, y=y_nl_wnll,marker_color=PLOTLYJS_COLORS[1], legendgroup="wnll", name="",showlegend=false, line_dash="dot", mode="lines"),
    ],
    Layout( 
        # title="Computational Cost Comparison for \$K\$",
        margin_t=0,
        xaxis=attr(title_text="K", type="log",title_font_size=20,tickfont_size=20),
        yaxis=attr(title_text="time (s)", type="log", title_font_size=20,tickfont_size=20),
        legend=attr(
            # x=.925,
            y=0.5,
            yanchor="center",
            font=(
                family="sans-serif",
                size=18,
                color="#000"
            ),
            bgcolor="#E2E2E2",
            bordercolor= "#636363",
            borderwidth= 2,
        ),
    )
)
display(p)
savefig(p,joinpath(FIG_DIR, "computationalCost_K.pdf"), width=550, height=400)
