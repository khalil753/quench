using HCubature, QuadGK

function memoized_integrate(l_or_m, Ms_Lss, initial_τs, final_τs, maxevals, rtol)
    if l_or_m in keys(Ms_Lss) return Ms_Lss[l_or_m] end
    Ms_Lss[l_or_m] = hcubature(l_or_m, initial_τs, final_τs, maxevals=maxevals, rtol=rtol)[1]
    # Ms_Lss[l_or_m] = quad_integrate(l_or_m, initial_τs, final_τs, rtol)
    return Ms_Lss[l_or_m]
  end
  
  function quad_integrate(f, t0s, tfs, rtol)
    f1(x) = quadgk(y -> f([x, y]), t0s[2], tfs[2], rtol=rtol)[1]
    return quadgk(f1, t0s[1], tfs[1], rtol=rtol)[1]
  end
  
struct MemoizedIntegrator{T} <: Function
    Ms_Lss    :: Dict
    initial_τs:: Vector{T}
    final_τs  :: Vector{T} 
    maxevals  :: Int
    rtol   :: AbstractFloat
end
MemoizedIntegrator(initial_τs, final_τs, maxevals, rtol) = MemoizedIntegrator(Dict(), initial_τs, final_τs, maxevals, rtol)
(MI::MemoizedIntegrator)(m_or_l) = memoized_integrate(m_or_l, MI.Ms_Lss, MI.initial_τs, MI.final_τs, MI.maxevals, MI.rtol)