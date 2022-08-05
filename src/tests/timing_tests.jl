using BenchmarkTools, ForwardDiff

include("../quench.jl")

Wττ′s, Dττ′ = create_distributions(_Ws[wightman_funct_name], _Ds[wightman_funct_name], 
                                   χ0A, χ0B, b)

χ(τ) = χs[switching_func_name](τ/σ)
df = deform_funcs[deform_func_name]
dist_func = distance_funcs["quench"]
m, lAA = get_m(Dττ′, λ, Ω, χ), get_l(Wττ′s["AA"], λ, Ω, χ)
lAA = complexify_l_or_m(lAA, df, dist_func, ε_contour)

function _time()
  # X, X′ = Vec{C}([1,2]), Vec{C}([1, -1])
  # @btime W₀($(C(1)), $(C(2)))
  # @btime _W_quench($X, $X′)
  # _D = _Ds["quench"]
  # @btime _D($X, $X′)

  # XA, XB = QuenchTrajectory(χ0A, b), QuenchTrajectory(χ0B, b)
  # xa, xb = Vec{Fl}(XA(1)), Vec{Fl}(XB(2))
  # W = DistributionWithTrajectories(_Ws["quench"] , XA, XB)
  # D = DistributionWithTrajectories(_Ds["quench"] , XA, XB)

  # @btime ($XA(1), $XB(1))
  # @btime _W_quench($xa, $xb)
  # @btime $W(1,2);# @btime $D(1,2)
  # @btime $Dττ′(1,2)
  @btime $lAA([0.2, 0.5])
end

_time()
