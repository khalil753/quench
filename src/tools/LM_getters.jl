include("misc.jl")

get_l(W, λ, Ω, χ) = τs ->  λ^2*χ(τs[1])*χ(τs[2]) * W(τs[1], τs[2])*exp(-im*Ω*(τs[1] - τs[2]))
get_m(D, λ, Ω, χ) = τs -> -λ^2*χ(τs[1])*χ(τs[2]) * D(τs[1], τs[2])*exp( im*Ω*(τs[1] + τs[2]))
function get_ls(Ws, λ, Ω, χ) 
  Dict("AA" => get_l(Ws["AA"], λ, Ω, χ),
       "BB" => get_l(Ws["BB"], λ, Ω, χ),
       "AB" => get_l(Ws["AB"], λ, Ω, χ))
end

get_Xs(W::DistributionWithTrajectories) = W.X, W.X′
get_Xs(Wττ′) = Wττ′.Rf.inner.f.X, Wττ′.Rf.inner.f.X′

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

