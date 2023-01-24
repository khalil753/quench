using HCubature, Plots, QuantumInformation, BenchmarkTools
include("../params.jl")
include("../tools/misc.jl")
# include("../tools/1_Trajectories.jl")
# include("../tools/2_D&W.jl")
include("../tools/3_LM_getters.jl")
include("../tools/SwitchingFuncs.jl")
include("../tools/DeformFuncs.jl")

Wττ′s, Dττ′, χ, df, dist_func = initialize_stuff()

integrate(f) = hcubature(f, [-3σ, -3σ], [3σ, 3σ], maxevals=50000, rtol=rtol)[1]
Cplify(l_or_m) = complexify_l_or_m(l_or_m, df, dist_func, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument
# Cplify(l_or_m) = complexify_l_or_m(l_or_m, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument
rotate(f) = τs -> f([τs[1]+τs[2], τs[2]-τs[1]])

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

  m, ls = get_m(Dττ′, λ, Ω, χ), get_ls(Wττ′s, λ, Ω, χ)
  println("Eval time of m and ls before any contour deformation")
  @btime $m([1,2])
  @btime $ls["AA"]([1,2])
  @btime $ls["AB"]([1,2])
  @btime $ls["BB"]([1,2])

  Cplify(l_or_m) = complexify_l_or_m(l_or_m, df, dist_func, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument
  m, ls = Cplify(m)           , map_dict(Cplify, ls)
  println("\nEval time of m and ls after contour deformation")
  @btime $m([1,2])
  @btime $ls["AA"]([1,2])
  @btime $ls["AB"]([1,2])
  @btime $ls["BB"]([1,2])
  @btime $m([1,2])

  Cplify(l_or_m) = complexify_l_or_m(l_or_m, ε_contour)
  m, ls = Cplify(m)           , map_dict(Cplify, ls)
  println("Eval time of m and ls after iε prescription")
  @btime $m([1,2])
  @btime $ls["AA"]([1,2])
  @btime $ls["AB"]([1,2])
  @btime $ls["BB"]([1,2])
  @btime $m([1,2])
  
  # X_flat = InertialTrajectory(0.0, 0.0, 0.0)
  # pole_distance_flat(τs) = distance_funcs["flat"](X_flat(τs[1]), X_flat(τs[2]))
  # X_quench = QuenchTrajectory(χ0A, b)
  # pole_distance_quench(τs) = distance_funcs["quench"](X_quench(τs[1]), X_quench(τs[2]))
  
  # @btime get_Δτ([1,2], $df, $pole_distance_flat  , $ε_contour)
  # @btime get_Δτ([1,2], $df, $pole_distance_quench, $ε_contour)

  # _get_Δτ_flat(τs) = get_Δτ(τs, df, pole_distance_flat, ε_contour)
  # get_∇Δτ_flat(τs) = ForwardDiff.gradient(_get_Δτ_flat, τs)
  # _get_Δτ_quench(τs) = get_Δτ(τs, df, pole_distance_quench, ε_contour)
  # get_∇Δτ_quench(τs) = ForwardDiff.gradient(_get_Δτ_quench, τs)

  # @btime $get_∇Δτ_flat(  [1,2])
  # @btime $get_∇Δτ_quench([1,2])



  # m, ls = get_m(Dττ′, λ, Ω, χ), get_ls(Wττ′s, λ, Ω, χ)
  # m, ls = Cplify(m)           , map_dict(Cplify, ls)
  # # @btime integrate($m)        
  # test = rotate(ls["AA"])
  # # @btime integrate($ls["AA"])        
  # @btime integrate($test)        
  # # @btime integrate($ls["AB"])        
  # # @btime map_dict($integrate, $ls)
end

function plot_l_and_m()
  X_flatA, X_flatB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(10.0, 0.0, 0.0)
  W_flat, D_flat = DistributionWithTrajectories(_Ws["flat"], X_flatA, X_flatA), DistributionWithTrajectories(_Ds["flat"], X_flatA, X_flatB)
  l_flat, m_flat = get_l(W_flat, λ, Ω, χ), get_m(D_flat, λ, Ω, χ)
  l_flat, m_flat = complexify_l_or_m(l_flat, df, distance_funcs["flat"], ε_contour), complexify_l_or_m(m_flat, df, distance_funcs["flat"], ε_contour)

  # @btime l_flat([1,1])
  # @btime integrate($l_flat)
  x = LinRange(-2,2,400)
  y = LinRange(-0.1,0.1,400)
  display(surface(x, y, (x,y) -> abs(l_flat([(x+y)/2,(x-y)/2])), camera=(85,25)))
  # x = LinRange(-1,1,400)
  # y = LinRange(-1,1,400)
  # display(surface(x, y, (x,y)->abs(m_flat([x,y])))#, camera=(85,25)))
end

# plot_l_and_m()
_time()
