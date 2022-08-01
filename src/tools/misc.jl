# Abbreviations
Fl, Vec, C = Float64, Vector, ComplexF64

is_before_quench(X) = X[1] <  0
is_after_quench(X)  = X[1] >= 0

function get_crossed_derivative(f::Function, ε::Fl) ::Function
  fxy(x::Fl, y::Fl)::C = (f(x+ε, y+ε) - f(x+ε, y-ε) - f(x-ε, y+ε) + f(x-ε, y-ε))/(4*ε^2)
end

map_dict(f::Function, d::Dict)::Dict = Dict([(k, f(v)) for (k, v) in d])

get_l(W, λ, Ω, χ) = τs ->  λ^2*χ(τs[1])*χ(τs[2]) * W(τs[1], τs[2])*exp(-im*Ω*(τs[1] - τs[2]))
get_m(D, λ, Ω, χ) = τs -> -λ^2*χ(τs[1])*χ(τs[2]) * D(τs[1], τs[2])*exp( im*Ω*(τs[1] + τs[2]))

function get_ls(Ws, λ, Ω, χ) 
  Dict("AA" => real ∘ get_l(Ws["AA"], λ, Ω, χ),
       "BB" => real ∘ get_l(Ws["BB"], λ, Ω, χ),
       "AB" =>        get_l(Ws["AB"], λ, Ω, χ))
end

function create_distributions(_Ws, _Ds, wightman_funct_name,
                              χ0A, χ0B, b, 
                              dist_ε, numeric_derivative_ε)
  _W, _D = _Ws[wightman_funct_name], _Ds[wightman_funct_name]
  Ws = Dict()
  Ws["AB"] = DistributionWithTrajectories(_W, χ0A, χ0B, b)
  Ws["AA"] = DistributionWithTrajectories(_W, χ0A, χ0A, b)
  Ws["BB"] = DistributionWithTrajectories(_W, χ0B, χ0B, b)

  _gcd(f::Function) = get_crossed_derivative(f, numeric_derivative_ε)
  Wττ′s = map_dict(_gcd, Ws)
  Dττ′  = _gcd(DistributionWithTrajectories(_D, χ0A, χ0B, b))
  return  Wττ′s,  Dττ′
end

function deform_trajectories(X::QuenchTrajectory, X′::QuenchTrajectory, deform_func::Function, ε_contour::Float64)
  function Z(τ, τ′)
    Δη, Δy = X(τ) - X′(τ′)
    if abs(Δη - Δy) < ε_deformation 
      iΔτ = deform_func(τ, τ′)
      return τ 
    else
      iΔτ = im*ε_deformation*deform_func(τ, τ′)/ε_deformation
      return τ + iΔτ
    end
  end
end