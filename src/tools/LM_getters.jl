include("misc.jl")

get_Xs(W::DistributionWithTrajectories) = W.X, W.X′
get_Xs(Wττ′) = Wττ′.Rf.inner.f.X, Wττ′.Rf.inner.f.X′

function get_m(Dττ′, λ, Ω, χ, dist_func,deform_func, ε_contour) 
  m(τs) = -λ^2*χ(τs[1])*χ(τs[2]) * Dττ′(τs[1], τs[2])*exp( im*Ω*(τs[1] + τs[2]))
  X, X′ = get_Xs(Dττ′)
  pole_distance(τs) = dist_func(X(τs[1]), X′(τs[2]))
  return complexify(m, deform_func, pole_distance, ε_contour)
end

function get_l(Wττ′, λ, Ω, χ, dist_func, deform_func, ε_contour) 
  l(τs) = λ^2*χ(τs[1])*χ(τs[2]) * Wττ′(τs[1], τs[2])*exp(-im*Ω*(τs[1] - τs[2]))
  X, X′ = get_Xs(Wττ′)
  pole_distance(τs) = dist_func(X(τs[1]), X′(τs[2]))
  return complexify(l, deform_func, pole_distance, ε_contour)
end

function get_ls(Ws, λ, Ω, χ, dist_func, deform_func, ε_contour) 
  Dict("AA" => real ∘ get_l(Ws["AA"], λ, Ω, χ, dist_func, deform_func, ε_contour),
       "BB" => real ∘ get_l(Ws["BB"], λ, Ω, χ, dist_func, deform_func, ε_contour),
       "AB" =>        get_l(Ws["AB"], λ, Ω, χ, dist_func, deform_func, ε_contour))
end