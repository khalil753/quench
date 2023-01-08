using HCubature

struct MemoizedIntegrator{T} <: Function
    Ms_Lss    :: Dict
    initial_τs:: Vector{T}
    final_τs  :: Vector{T} 
    maxevals  :: Int
    int_tol   :: AbstractFloat
end
MemoizedIntegrator(initial_τs, final_τs, maxevals, int_tol) = MemoizedIntegrator(Dict(), initial_τs, final_τs, maxevals, int_tol)
(MI::MemoizedIntegrator)(m_or_l) = memoized_integrate(m_or_l, MI.Ms_Lss, MI.initial_τs, MI.final_τs, MI.maxevals, MI.int_tol)