# Abbreviations
Fl, Vec, C = AbstractFloat, Vector, Complex

function get_crossed_derivative(f::Function, ε::Fl) ::Function
  fxy(x::Fl, y::Fl)::C = (f(x+ε, y+ε) - f(x+ε, y-ε) - f(x-ε, y+ε) + f(x-ε, y-ε))/(4*ε^2)
end

map_dict(f::Function, d::Dict)::Dict = Dict([(k, f(v)) for (k, v) in d])

get_l(W, λ, Ω, χ) = τs ->  λ*χ(τs[1])*χ(τs[2]) * W(τs[1], τs[2])*exp(-im*Ω*(τs[1] - τs[2]))
get_m(D, λ, Ω, χ) = τs -> -λ*χ(τs[1])*χ(τs[2]) * D(τs[1], τs[2])*exp( im*Ω*(τs[1] + τs[2]))
get_ls(Ws, λ, Ω, χ) = Dict([(k, get_l(W, λ, Ω, χ)) for (k,W) in Ws])