# include("misc.jl")
# include("1_Trajectories.jl")
include("2_D&W.jl")

get_Xs(W::DistributionWithTrajectories) = W.X, W.X′
function get_Xs(Wττ′) 
  try 
    Wττ′.Rf.inner.f.X, Wττ′.Rf.inner.f.X′
  catch e
    if isa(e, ErrorException)
      Wττ′.f.X, Wττ′.f.X′
    end
  end
end

get_γ(X::InertialTrajectory) = τ -> 1
get_γ(X::QuenchTrajectory)   = τ -> 1/X.χ0

function get_l(W, λ, Ω, χ)  
  γA, γB = get_γ.(get_Xs(W))
  l(τs) =  λ^2*χ(τs[1])*χ(τs[2]) * W(τs[1], τs[2])*exp(-im*Ω*(τs[1] - τs[2])) * γA(τs[1])*γB(τs[2])
end

function get_m(D, λ, Ω, χ) 
  γA, γB = get_γ.(get_Xs(D))
  m(τs) = -λ^2*χ(τs[1])*χ(τs[2]) * D(τs[1], τs[2])*exp( im*Ω*(τs[1] + τs[2])) * γA(τs[1])*γB(τs[2])
end 

function get_ls(Ws, λ, Ω, χ) 
  Dict("AA" => get_l(Ws["AA"], λ, Ω, χ),
       "BB" => get_l(Ws["BB"], λ, Ω, χ),
       "AB" => get_l(Ws["AB"], λ, Ω, χ))
end

function get_pole_distance(l_or_m, dist_func)
  local W
  try 
    W = l_or_m.W
  catch e
    if isa(e, ErrorException) 
      W = l_or_m.D 
    else 
      throw(e)
    end
  end
  X, X′ = get_Xs(W)
  pole_distance(τs) = dist_func(X(τs[1]), X′(τs[2]))
end

function complexify_l_or_m(l_or_m, deform_func, dist_func, ε)
  pole_distance = get_pole_distance(l_or_m, dist_func)
  complexify(l_or_m, deform_func, pole_distance, ε)   
end
complexify_l_or_m(l_or_m, ε) = complexify(l_or_m, ε)   

function complexify_ls(ls, deform_func, dist_func, ε)
  Dict([(k,complexify_l_or_m(l, deform_func, dist_func, ε)) for (k,l) in ls])
end
complexify_ls(ls, ε) = Dict([(k,complexify_l_or_m(l, ε)) for (k,l) in ls])

