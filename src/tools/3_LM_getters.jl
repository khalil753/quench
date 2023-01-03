# include("misc.jl")
# include("1_Trajectories.jl")
include("2_D&W.jl")
include("ComplexifyFuncs.jl")

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

function get_l(W, λ, Ω, χs, ε_contour)
  χD , χD′ = isa(χs, Vector) ? χs : (χs, χs)
  γD′, γD  = get_γ.(get_Xs(W))
  W = complexify(W, ε_contour)
  l(τs) =  λ^2*χD(τs[1])*χD′(τs[2]) * W(τs[2], τs[1])*exp(im*Ω*(τs[1] - τs[2])) * γD(τs[1])*γD′(τs[2])
end

function get_m(D, λ, Ω, χs, ε_contour)
  χA, χB = isa(χs, Vector) ? χs : (χs, χs)
  γA, γB = get_γ.(get_Xs(D))
  D = complexify(D, ε_contour)
  m(τs) = -λ^2*χA(τs[1])*χB(τs[2]) * D(τs[1], τs[2])*exp(im*Ω*(τs[1] + τs[2])) * γA(τs[1])*γB(τs[2])
end 

function get_ls(Ws, λ, Ω, χs, ε_contour) 
  χA, χB = isa(χs, Vector) ? χs : (χs, χs)
  Dict("AA" => get_l(Ws["AA"], λ, Ω, [χA, χA], ε_contour),
       "BB" => get_l(Ws["BB"], λ, Ω, [χB, χB], ε_contour),
       "AB" => get_l(Ws["AB"], λ, Ω, [χA, χB], ε_contour))
end

# function get_pole_distance(l_or_m, dist_func)
#   local W
#   try 
#     W = l_or_m.W
#   catch e
#     if isa(e, ErrorException) 
#       W = l_or_m.D 
#     else 
#       throw(e)
#     end
#   end
#   X, X′ = get_Xs(W)
#   pole_distance(τs) = dist_func(X(τs[1]), X′(τs[2]))
# end

# function complexify_l_or_m(l_or_m, deform_func, dist_func, ε)
#   """ I had to define this function to extract the pole distance of out l/m to then pass this to complexify"""
#   pole_distance = get_pole_distance(l_or_m, dist_func)
#   complexify(l_or_m, deform_func, pole_distance, ε)   
# end
# complexify_l_or_m(l_or_m, ε) = complexify(l_or_m, ε)   

# function complexify_ls(ls, deform_func, dist_func, ε)
#   Dict([(k,complexify_l_or_m(l, deform_func, dist_func, ε)) for (k,l) in ls])
# end
# complexify_ls(ls, ε) = Dict([(k,complexify_l_or_m(l, ε)) for (k,l) in ls])

