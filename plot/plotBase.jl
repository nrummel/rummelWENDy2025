using JLD2, Latexify, ColorSchemes, Colors, LaTeXStrings # for plotting and saving captions to tex
using PlotlyJS, Statistics, Printf
using PlotlyKaleido: restart as restart_kaleido
restart_kaleido(plotly_version = "2.35.2", mathjax = true) 
using PlotlyJS: savefig, scatter, Layout, AbstractTrace, attr
import PlotlyJS.plot as plotjs
Latexify.set_default(fmt = "%.4g")
##
FIG_DIR = joinpath(@__DIR__, "../../NonLinearWENDyPaper/figs")
DATA_DIR = joinpath(@__DIR__, "../data")
PLOTLYJS_COLORS = [
    colorant"#1f77b4",  # muted blue
    colorant"#ff7f0e",  # safety orange
    colorant"#2ca02c",  # cooked asparagus green
    colorant"#d62728",  # brick red
    colorant"#9467bd",  # muted purple
    colorant"#8c564b",  # chestnut brown
    colorant"#e377c2",  # raspberry yogurt pink
    colorant"#7f7f7f",  # middle gray
    colorant"#bcbd22",  # curry yellow-green
    colorant"#17becf"   # blue-teal
]
CU_BOULDER_COLORS = [
    colorant"#000000", # black 
    colorant"#CFB87C", # gold 
    colorant"#565A5C", # dark grey 
    colorant"#A2A4A3", # light grey 
]
algo2Color = (
    init=CU_BOULDER_COLORS[1],
    truth=CU_BOULDER_COLORS[3],
    oels=PLOTLYJS_COLORS[1],
    wls=PLOTLYJS_COLORS[2],
    wendy_irls=PLOTLYJS_COLORS[3],
    wendy_mle=CU_BOULDER_COLORS[2],
    hybrid=PLOTLYJS_COLORS[6],
)
algo2Disp = (
    init="w₀ Solution",
    truth="w* Solution",
    oels="OE-LS",#"Forward Solve Nonlinear Least Squares",
    wls="WLS", 
    wendy_irls="WENDy-IRLS",
    wendy_mle="WENDy-MLE",
    hybrid="Hybrid",
    # hybrid_wlsq_trustRegion="Hybrid WLS-WENDy"
)
ex2Disp = Dict(
    "sir"=>"SIR-TDI",
    "goodwin"=>"Goodwin",
    "robertson"=>"Robertson",
    "hindmarshRose"=>"Hindmarsh Rose", 
    "logisticGrowth"=>"Logistic Growth", 
    "multimodal"=>"Goodwin 2D", 
    "lorenz"=>"Lorenz System", 
);