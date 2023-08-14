include("tools/misc.jl")
include("params.jl")

χ2y(χ) = asinh(χ/b)
χ0B = χ0Bs[1]
y0A, y0B = χ2y.([χ0A, χ0B])

initial_τs, final_τs = initialize_integration_ranges(ηcA_or_τcA, ηcB_or_τcB, χ0A, χ0B, b, using_ηs, space_time, switching_func_name)

sqrtA, sqrtB = sqrt(χ0A^2+b^2), sqrt(χ0B^2+b^2)
η0A, ηfB = initial_τs[1]/sqrtA, final_τs[2]/sqrtB
η0B, ηfA = initial_τs[2]/sqrtB, final_τs[1]/sqrtA

if space_time=="quench"
  function are_spacelike(X,X′)
    (η, y), (η′, y′) = X, X′
    return abs(y - y′) > abs(η - η′)
  end
  X0A, X0B, XfA, XfB = [η0A, y0A], [η0B, y0B], [ηfA, y0A], [ηfB, y0B]

  if are_spacelike(X0A, XfB) println("Detector A's initial position and detector B's final   one ARE     spacelike separated") 
  else                       println("Detector A's initial position and detector B's final   one are NOT spacelike separated") end
  if are_spacelike(X0B, XfA) println("Detector A's final   position and detector B's initial one ARE     spacelike separated") 
  else                       println("Detector A's final   position and detector B's initial one are NOT spacelike separated") end
else
  are_spacelike(X, Y) = (X[1] - Y[1])^2 - sum((X[2:end] - Y[2:end]).^2) < 0

  X0A = [χ0A*sinh(η0A), χ0A*cosh(η0A)]
  XfA = [χ0A*sinh(ηfA), χ0A*cosh(ηfA)]
  X0B = [χ0B*sinh(η0B), χ0B*cosh(η0B)]
  XfB = [χ0B*sinh(ηfB), χ0B*cosh(ηfB)]

  if are_spacelike(X0A, XfB) println("Detector A's initial position and detector B's final   one ARE spacelike separated") 
  else                       println("Detector A's initial position and detector B's final   one are NOT spacelike separated") end
  if are_spacelike(X0B, XfA) println("Detector A's final   position and detector B's initial one ARE spacelike separated") 
  else                       println("Detector A's final   position and detector B's initial one are NOT spacelike separated") end
end


