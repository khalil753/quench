using HCubature, Plots, SpecialFunctions, ProgressBars
include("params.jl")
include("tools/SwitchingFuncs.jl")
include("tools/2_D&W.jl")

struct _W <: Function 
    ε::Float64
    L::Float64
end

function (W::_W)(t, t′)
    iε = im*W.ε
    L = W.L
    -1/(4*π^2*((t - t′- iε)^2 - L^2))
end

θ(t) = t>=0 ? 1.0 : 0.0
χ = χs["gauss"]

Ls = LinRange(1/10, 5, 20)
Ms = []
for L in Ls    
    ε = 1e-2
    # W = _W(ε, L)
    W = _Ws["flat"]
    m(t, t′) = χ(t)*χ(t′) * exp(-im*(t+t′)) * (W([t-im*ε,0],[t′, L])*θ(t-t′) + W([t′, L], [t-im*ε,0])*θ(t′-t)) 
    # m(t, t′) = χ(t)*χ(t′) * exp(-im*(t+t′)) * (W(t,t′)*θ(t-t′) + W(t′,t)*θ(t′-t)) 
    m(ts) = m(ts...)
    M_num = hcubature(m, [-5,-5], [5,5], maxevals=70000, rtol=int_tol)[1]
    push!(Ms, M_num)
end

M(L) = im/(4*√π*L)*exp(-1 - L^2/4)*(erf(im*L/2) - 1)

p = plot(Ls, [(abs ∘ M).(Ls), abs.(Ms)], labels=["theoretical" "numerical"], ylims=[-1e-6, max(abs.(Ms)...)*3/2])
display(p)
