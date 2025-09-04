includet("util.jl")
using BenchmarkTools
##
@info "M Experiment "
savefile = joinpath(DATA_DIR, "computationalCostExperiment_M.jld2")
Ms = [Int(2^e) for e in 8:13]
l_wnll_times = zeros(length(Ms))
l_grad_times = zeros(length(Ms))
l_hess_times = zeros(length(Ms))
nl_wnll_times = zeros(length(Ms))
nl_grad_times = zeros(length(Ms))
nl_hess_times = zeros(length(Ms))
K = 100
for (i,M) in enumerate(Ms)
    @info "$i, $M"
    params = WENDyParameters(
        radiiParams       = [4,6,8],
        maxTestFunCondNum = Inf,
        minTestFunInfoNum = 0,
        Kmax              = K
    )
    global ex_l = SimulatedWENDyData(
        "hindmarshRose", 
        HINDMARSH_f!,
        HINDMARSH_INIT_COND,
        HINDMARSH_TRNG,
        HINDMARSH_WSTAR,
        HINDMARSH_WRNG,
        params;
        linearInParameters=Val(true),
        noiseDist=Val(Normal),
        Mp1=M+1
    );
    global ex_nl = SimulatedWENDyData(
        "hindmarshRose", 
        HINDMARSH_f!,
        HINDMARSH_INIT_COND,
        HINDMARSH_TRNG,
        HINDMARSH_WSTAR,
        HINDMARSH_WRNG,
        params;
        linearInParameters=Val(false),
        noiseDist=Val(Normal),
        Mp1=M+1
    );
    simParams = SimulationParameters(
        noiseRatio=0.1,
        seed=1,
        timeSubsampleRate=1, 
        corruptU0=true
    )
    @info " Building WENDyProblems..."
    @time "linear problem" global l_wendyProb, _, _ = runExample(ex_l, simParams, [:init],ll=Warn)
    @time "nonlinear problem" global nl_wendyProb, _, _ = runExample(ex_nl, simParams, [:init],ll=Warn)
    @info "  Checking for size discrepensies"
    @assert l_wendyProb.K == nl_wendyProb.K == K
    @assert l_wendyProb.Mp1 == nl_wendyProb.Mp1 == M+1
    @assert l_wendyProb.J == nl_wendyProb.J == 10
    @assert l_wendyProb.D == nl_wendyProb.D == 3 
    global l_wnll    = l_wendyProb.wnll.f
    global l_∇wnll!  = l_wendyProb.wnll.∇f!
    global l_Hwnll!  = l_wendyProb.wnll.Hf!
    global nl_wnll   = nl_wendyProb.wnll.f
    global nl_∇wnll! = nl_wendyProb.wnll.∇f!
    global nl_Hwnll! = nl_wendyProb.wnll.Hf!
    ##
    global w = ex.wTrue
    global J = l_wendyProb.J
    global g = similar(w)
    global H = zeros(J,J)
    @info "  Calling once to get precompilation out of the way..."
    l_wnll(w)
    l_wnll(w)
    l_∇wnll!(g, w)
    nl_wnll(w)
    nl_wnll(w)
    nl_∇wnll!(g, w)
    ## 
    @info "  Calling for the record"
    @info "     linear wnll"
    l_wnll_times[i] = @belapsed l_wnll(w) #@belapsed
    @info "      $(@sprintf "%.4g" l_wnll_times[i] ) s"
    @info "     nonlinear wnll"
    nl_wnll_times[i] = @belapsed nl_wnll(w) #@belapsed
    @info "      $(@sprintf "%.4g" nl_wnll_times[i] ) s"
    @info "     linear ∇wnll!"
    l_grad_times[i] = @belapsed l_∇wnll!(g, w) #@belapsed
    @info "      $(@sprintf "%.4g" l_grad_times[i] ) s"
    @info "     nonlinear ∇wnll!"
    nl_grad_times[i] = @belapsed nl_∇wnll!(g, w) #@belapsed
    @info "      $(@sprintf "%.4g" nl_grad_times[i] ) s"
    @info "     linear Hwnll!"
    l_hess_times[i] = @belapsed l_Hwnll!(H, w) #@belapsed
    @info "      $(@sprintf "%.4g" l_hess_times[i] ) s"
    @info "     nonlinear Hwnll!"
    nl_hess_times[i] = @belapsed nl_Hwnll!(H, w) #@belapsed
    @info "      $(@sprintf "%.4g" nl_hess_times[i] ) s"
end
JLD2.save(
    savefile, 
    Dict( 
        "Ms"=>Ms,
        "l_wnll_times"=>l_wnll_times,
        "l_grad_times"=>l_grad_times,
        "l_hess_times"=>l_hess_times,
        "nl_wnll_times"=>nl_wnll_times,
        "nl_grad_times"=>nl_grad_times,
        "nl_hess_times"=>nl_hess_times,
    )
)



##
@info "K Experiment "
savefile = joinpath(DATA_DIR, "computationalCostExperiment_K.jld2")
Ks = [2 ^e for e in 6:10]
l_wnll_times = zeros(length(Ks))
l_grad_times = zeros(length(Ks))
l_hess_times = zeros(length(Ks))
nl_wnll_times = zeros(length(Ks))
nl_grad_times = zeros(length(Ks))
nl_hess_times = zeros(length(Ks))
M = 2000
for (i,K) in enumerate(Ks)
    @info "i=$i/$(length(Ks)), K=$K, M=$M"
    params = WENDyParameters(
        radiiParams       = [4,6,8],
        maxTestFunCondNum = Inf,
        minTestFunInfoNum = 0,
        Kmax              = K
    )
    global ex_l = SimulatedWENDyData(
        "hindmarshRose", 
        HINDMARSH_f!,
        HINDMARSH_INIT_COND,
        HINDMARSH_TRNG,
        HINDMARSH_WSTAR,
        HINDMARSH_WRNG,
        params;
        linearInParameters=Val(true),
        noiseDist=Val(Normal),
        Mp1=M+1
    );
    global ex_nl = SimulatedWENDyData(
        "hindmarshRose", 
        HINDMARSH_f!,
        HINDMARSH_INIT_COND,
        HINDMARSH_TRNG,
        HINDMARSH_WSTAR,
        HINDMARSH_WRNG,
        params;
        linearInParameters=Val(false),
        noiseDist=Val(Normal),
        Mp1=M+1
    );
    simParams = SimulationParameters(
        noiseRatio=0.1,
        seed=1,
        timeSubsampleRate=1, 
        corruptU0=true
    )
    @info " Building WENDyProblems..."
    @time "linear problem" global l_wendyProb, _, _ = runExample(ex_l, simParams, [:init],ll=Warn)
    @time "nonlinear problem" global nl_wendyProb, _, _ = runExample(ex_nl, simParams, [:init],ll=Warn)
    @info "  Checking for size discrepensies"
    @assert l_wendyProb.K == nl_wendyProb.K == K
    @assert l_wendyProb.Mp1 == nl_wendyProb.Mp1 == M+1
    @assert l_wendyProb.J == nl_wendyProb.J == 10
    @assert l_wendyProb.D == nl_wendyProb.D == 3 
    @info "  Building Cost Function" 
    global l_wnll    = l_wendyProb.wnll.f
    global l_∇wnll!  = l_wendyProb.wnll.∇f!
    global l_Hwnll!  = l_wendyProb.wnll.Hf!
    global nl_wnll   = nl_wendyProb.wnll.f
    global nl_∇wnll! = nl_wendyProb.wnll.∇f!
    global nl_Hwnll! = nl_wendyProb.wnll.Hf!
    ##
    global w = ex.wTrue
    global J = l_wendyProb.J
    global g = similar(w)
    global H = zeros(J,J)
    @info "  Calling once to get precompilation out of the way..."
    l_wnll(w)
    nl_wnll(w)
    l_∇wnll!(g, w)
    nl_∇wnll!(g, w)
    l_Hwnll!(H, w)
    nl_Hwnll!(H, w)
    ## 
    @info "  Calling for the record"
    @info "     linear wnll"
    l_wnll_times[i] = @belapsed l_wnll(w)
    @info "      $(@sprintf "%.4g" l_wnll_times[i] ) s"
    @info "     nonlinear wnll"
    nl_wnll_times[i] = @belapsed nl_wnll(w)
    @info "      $(@sprintf "%.4g" nl_wnll_times[i] ) s"
    @info "     linear ∇wnll!"
    l_grad_times[i] = @belapsed l_∇wnll!(g, w)
    @info "      $(@sprintf "%.4g" l_grad_times[i] ) s"
    @info "     nonlinear ∇wnll!"
    nl_grad_times[i] = @belapsed nl_∇wnll!(g, w)
    @info "      $(@sprintf "%.4g" nl_grad_times[i] ) s"
    @info "     linear Hwnll!"
    l_hess_times[i] = @belapsed l_Hwnll!(H, w)
    @info "      $(@sprintf "%.4g" l_hess_times[i] ) s"
    @info "     nonlinear Hwnll!"
    nl_hess_times[i] = @belapsed nl_Hwnll!(H, w)
    @info "      $(@sprintf "%.4g" nl_hess_times[i] ) s"
end
JLD2.save(
    savefile, 
    Dict( 
        "Ks"=>Ks,
        "l_wnll_times"=>l_wnll_times,
        "l_grad_times"=>l_grad_times,
        "l_hess_times"=>l_hess_times,
        "nl_wnll_times"=>nl_wnll_times,
        "nl_grad_times"=>nl_grad_times,
        "nl_hess_times"=>nl_hess_times,
    )
)
