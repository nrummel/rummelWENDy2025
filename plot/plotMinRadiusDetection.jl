includet("../src/util.jl")
includet("plotBase.jl")
using WENDy: _minRadius
##
# hyperparameters
sub = 2
radiusMinTime = 0.01
radiusMaxTime = 1
numRadii = 100
testFunSubRate = 2.0 ## timeSubSample Rate
testFunReductionScale = 2
Kᵣ = nothing
# data 
f! = LOGISTIC_GROWTH.f!
pstar = LOGISTIC_GROWTH.wTrue
tt = LOGISTIC_GROWTH.tt_full[1:sub:end]
U_exact = LOGISTIC_GROWTH.U_full[1:sub:end,:]
Mp1, D = size(LOGISTIC_GROWTH.U_full)
dt = mean(diff(tt))
function f(u,t)
    du = similar(u)
    f!(du, u, pstar, t)
    du
end
ff_exact = reduce(vcat, f(U_exact[m,:], t) for (m,t) in enumerate(tt))
dt = tt[2] - tt[1]
nothing
## Mininum Radius detection 
@info "Minimum Radius Detection"
soltrs = AbstractTrace[]
trs = AbstractTrace[]
shapes = []
noiseRatios = [0, 0.01, 0.1]
uu = Vector(undef, length(noiseRatios))
for (i,σ) in enumerate(noiseRatios)
    uu[i],_,_ = generateNoise(U_exact, σ, 1, true, Val(Normal)) 
    # for d in 1:D
    #     push!(soltrs, scatter(x=tt, y=uu[i][:,d],  name="σ=$(@sprintf "%.2g" σ), d=$d"))
    # end
    radiusMin = Int(max(ceil(radiusMinTime/dt), 2))
    radiusMax = Int(floor(radiusMaxTime/dt))
    ix, radii, errs = _minRadius(
        uu[i], dt, 
        radiusMin, numRadii, radiusMax, testFunSubRate, nothing; 
        debug=true
    )
    rstar = radii[ix]
    push!(trs, 
        scatter(
            x=dt*radii, 
            y=errs, 
            name="σ = $(@sprintf "%.2g" σ)", 
            marker_color=PLOTLYJS_COLORS[i+2]
        )
    )
    push!(shapes,attr(
        type="line",
        x0=rstar*dt,
        y0=100,
        x1=rstar*dt,
        y1=1e-12,
        line= attr(
                dash="dash",
                color=PLOTLYJS_COLORS[i+2],
                width= 3
            )
        ),
    )
end 
p_minRad = plot(
    trs, 
    Layout(
        title="Minimum Radius Selection", 
        yaxis_type="log", 
        yaxis_range=[-9, 1.5],
        xaxis_type="log", 
        shapes=shapes, 
        xaxis_title="Radius (s)", 
        yaxis_title=L"\text{log10}\Big(\hat{\mathcal{F}}\big[\lfloor \tfrac{M}{2}\rfloor\big]\Big)", 
        include_mathjax=true
    )
)
savefig(p_minRad, joinpath(FIG_DIR, "minRadiusDetection.pdf"))
display(p_minRad)
