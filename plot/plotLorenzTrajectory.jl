## Load everything we need 
includet("../src/util.jl")
includet("plotBase.jl")
##
params = WENDyParameters(;
    optimMaxiters  = 200,
    optimTimelimit = 200.0,
    radiusMinTime  = 0.02,
    radiusMaxTime  = 5.,
    Kᵣ             = 10,
    Kmax           = 500,
)
ex = SimulatedWENDyData(
    "lorenz", 
    LORENZ_f!,
    LORENZ_INIT_COND,
    (0.0,30.0),
    LORENZ_WSTAR,
    [    
        (5.0,15.0),
        (23.0,33.0),
        (2.0,3.0)
    ],
    params;
    dt = 0.01, 
    Mp1=nothing,
    linearInParameters=Val(false),
    noiseDist=Val(Normal),
    forceOdeSolve=true,
    ll=Warn
);

simParams = SimulationParameters(
    noiseRatio=0.1,
    seed=2,
    timeSubsampleRate=1, 
    corruptU0=true
)
simulate!(ex, simParams)

p = plot(
    [
        scatter3d(
            x=ex.U[][:,1],
            y=ex.U[][:,2],
            z=ex.U[][:,3],
            name=L"\mathbf{u}", 
            mode="markers" ,
            marker_color=CU_BOULDER_COLORS[4], 
            marker_size=5,
            marker_opacity=0.5,
            marker_line=attr(width=1,color="black"),
            showlegend=false
        ),
        scatter3d(
            x=ex.Ustar[][:,1],
            y=ex.Ustar[][:,2],
            z=ex.Ustar[][:,3],
            name=L"\mathbf{u}^*", 
            mode="lines",
            showlegend=false ,
            line_color=CU_BOULDER_COLORS[1], 
        )
    ], 
    Layout(
        scene = attr(
            xaxis_title=attr(
                text="X",
                font_size=25,
            ),
            xaxis_tickfont_size=15,
            yaxis_title=attr(
                text="Y",
                font_size=25,
            ),
            yaxis_tickfont_size=15,
            zaxis_title=attr(
                text="Z",
                font_size=25,
            ),
            zaxis_tickfont_size=15,
            camera=attr(eye= attr(x=-2, y=2, z = 1))
        )
    )
)
##
display(p)
plot_file = joinpath(FIG_DIR, "lorenz_solution3d.pdf")
@info "saving to $plot_file"
savefig(p, plot_file, height=600, width=600)