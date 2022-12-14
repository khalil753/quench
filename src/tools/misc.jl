# Abbreviations
using ForwardDiff

Fl, Vec, C = Float64, Vector, ComplexF64
VecCFl = Union{Vector{C}, Vector{Fl}}

struct NotImplmented <:Exception
  str::String
  function NotImplmented(str)
    println(str)    
  end
end

is_before_quench(X) = real(X[1]) <  0
is_after_quench(X)  = real(X[1]) >= 0

map_dict(f::Function, d::Dict)::Dict = Dict([(k, f(v)) for (k, v) in d])

function function_names()
  func_names = names(Main)
  filter!(name -> isa(eval(name),Function), func_names)
  return func_names
end

function make_Ws_dict()
  func_names = function_names()
  _Ws = Dict()
  for name in func_names
    name_str = String(name)
    if length(name_str) < 3 continue end
    if name_str[1:3] != "_W_" continue end
    spacetime = split(name_str, "_")[3]
    _Ws[spacetime] = eval(name)
  end
  return _Ws
end

get_crossed_derivative(f, ε) = (x,y) -> 1/4ε^2*(f(x + ε, y + ε) -  f(x - ε, y + ε) - f(x + ε, y - ε) + f(x - ε, y - ε))

function initialize_trajs(space_time, χ0A, χ0B, b)
  if     space_time=="quench"  XA, XB = QuenchTrajectory(χ0A, b)   , QuenchTrajectory(χ0B, b)
  elseif space_time=="flat"    XA, XB = InertialTrajectory(0, 0, 0), InertialTrajectory(χ0B, 0, 0) 
  elseif space_time=="rindler" XA, XB = AcceleratedTrajectory(χ0A,0,0,0), AcceleratedTrajectory(χ0B,0,0,0) end
  return XA, XB
end

function initialize_Ws(W_func, XA, XB)
  Ws = Dict()
  Ws["AB"] = DistributionWithTrajectories(W_func, XB, XA)
  Ws["AA"] = DistributionWithTrajectories(W_func, XA, XA)
  Ws["BB"] = DistributionWithTrajectories(W_func, XB, XB)
  return Ws
end

initialize_distributions(W_func::Function, D_func::Function, XA, XB) = (initialize_Ws(W_func, XA, XB), DwT(D_func, XA, XB))

function add_crossed_derivatives(Ws, D, ε_numeric_derivative)
  _gcd(f) = get_crossed_derivative(f, ε_numeric_derivative)
  return map_dict(_gcd, Ws), _gcd(D)
end

function initialize_distributions(space_time::String, χ0A::Fl, χ0B::Fl, b::Fl, with_derivative_coupling, ε_numeric_derivative)
  XA, XB = initialize_trajs(space_time, χ0A, χ0B, b)
  Ws, D  = initialize_distributions(_Ws[space_time], _Ds[space_time], XA, XB)
  if with_derivative_coupling Ws, D = add_crossed_derivatives(Ws, D, ε_numeric_derivative) end
  return Ws, D
end

function memoized_integrate(l_or_m, Ms_Lss, initial_τs, final_τs, maxevals, int_tol)
  if l_or_m in keys(Ms_Lss) return Ms_Lss[l_or_m] end
  Ms_Lss[l_or_m] = hcubature(l_or_m, initial_τs, final_τs, maxevals=maxevals, rtol=int_tol)[1]
  return Ms_Lss[l_or_m]
end

function get_ρ(M, Ls) 
  #                         00        01              10        11
  ρ = [1 - Ls["AA"] - Ls["BB"]         0               0   conj(M); #00
                             0  Ls["BB"]   conj(Ls["AB"])        0; #01
                             0  Ls["AB"]        Ls["AA"]         0; #10
                             M         0               0         0] #11

  if imag(ρ[2,2]) >   1e-5 println("the imaginary part of PB is quite big so there could be a numerical problem: PB = $(ρ[2,2])") end  
  if imag(ρ[3,3]) >   1e-5 println("the imaginary part of PA is quite big so there could be a numerical problem: PB = $(ρ[3,3])") end  
  if real(ρ[2,2]) <= -1e-3 println("PB is quite negative so there could be a numerical problem: PB = $(ρ[2,2])") end
  if real(ρ[3,3]) <= -1e-3 println("PA is quite negative so there could be a numerical problem: PA = $(ρ[3,3])") end
  if real(ρ[2,2]) >   1    println("PB is too big so there could be a numerical problem: PB = $(ρ[2,2])") end
  if real(ρ[3,3]) >   1    println("PA is too big so there could be a numerical problem: PA = $(ρ[3,3])") end
  return ρ
end

function negativity(ρ) 
  ρ22, ρ33 = real(ρ[2,2]), real(ρ[3,3])
  ρ22, ρ33 = max(ρ22, 0), max(ρ33, 0)
  max(0, sqrt(abs2(ρ[1,4]) + ((ρ22 - ρ33)/2)^2) - (ρ22 + ρ33)/2)
end

function concurrence(ρ)
  ρ22, ρ33 = real(ρ[2,2]), real(ρ[3,3])
  ρ22, ρ33 = max(ρ22, 0), max(ρ33, 0)
  2*max(0, abs(ρ[1,4]) - sqrt(ρ22*ρ33))
end

function make_img(χ0Bs, Ωs, Cs)
  p = contourf(χ0Bs, Ωs, Cs, linewidth=-0.0, xlabel="χB₀", ylabel="Ω", title="Concurrence")
  display(p)
  img_name = replace("$(now())", ":"=>"_")
  savefig(p, "plots/new_plots/$(img_name).png")
  return img_name
end

function store_in_df(params, img_name, run_duration)
  new_row = insertcols!(params, 1, "Image_Name"   => [img_name],
                                   "Run_Duration" => [Int(round(run_duration/60))])
  println(new_row)
  append = "df.csv" in readdir("plots")
  CSV.write("plots/df.csv", new_row, append=append)
end