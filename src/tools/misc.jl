# Abbreviations
using ForwardDiff, CairoMakie

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
  filter!(name -> name != :run_duration, func_names)
  filter!(name -> isa(eval(name),Function), func_names)
  return func_names
end

function make_Ws_dict()
  func_names = function_names()
  _Ws = Dict()
  for name in func_names
    name_str = String(name)
    if length(name_str) < 3 || name_str=="run_duration" continue end
    if name_str[1:3] != "_W_" continue end
    spacetime = join(split(name_str, "_")[3:end], "_")
    _Ws[spacetime] = eval(name)
  end
  return _Ws
end

get_crossed_derivative(f, ε) = (x,y) -> 1/4ε^2*(f(x + ε, y + ε) -  f(x - ε, y + ε) - f(x + ε, y - ε) + f(x - ε, y - ε))

function check_boundaries_above_quench(initial_τs)
  if space_time != "quench" return end
  if initial_τs[1] == σ*1e-2 println("Detector A's lower integration boundary is touching the quench boundary\n") end
  if initial_τs[2] == σ*1e-2 println("Detector B's lower integration boundary is touching the quench boundary\n") end
end

function check_boundaries_under_lightcone(final_τs, χ0A, χ0B, b, sqrtA, sqrtB)
  ηfA, ηfB = final_τs[1]/sqrtA, final_τs[2]/sqrtB
  y0A, y0B = asinh(χ0A/b)     , asinh(χ0B/b)
  if ηfA >= abs(y0A) println("Detector A's upper integration boundary is touching the light cone: ηfA = $ηfA >= $y0A = y0A \n") end
  if ηfB >= abs(y0B) println("Detector B's upper integration boundary is touching the light cone: ηfB = $ηfB >= $y0B = y0B \n") end
end

function initialize_integration_ranges(ηcA_or_τcA, ηcB_or_τcB, χ0A, χ0B, b, using_ηs, space_time, switching_func_name)
  sqrtA, sqrtB = sqrt(χ0A^2 + b^2), sqrt(χ0B^2 + b^2)
  if using_ηs τcA, τcB = ηcA_or_τcA*sqrtA, ηcB_or_τcB*sqrtB
  else        τcA, τcB = ηcA_or_τcA      , ηcB_or_τcB end
  Δτ = switching_func_name == "cos4" ? 0.5σ : 5σ
  if   space_time=="quench"
    initial_τs, final_τs = [max(σ*1e-2, τcA - Δτ), max(σ*1e-2, τcB - Δτ)],
                           [            τcA + Δτ ,             τcB + Δτ ]
    check_boundaries_above_quench(initial_τs)                                                    
  else initial_τs, final_τs = [τcA - Δτ, τcB - Δτ], [τcA + Δτ, τcB + Δτ] end  
  
  check_boundaries_under_lightcone(final_τs, χ0A, χ0B, b, sqrtA, sqrtB)
  return initial_τs, final_τs
end

initialize_integration_ranges(ηcA_or_τcA, ηcB_or_τcB, χ0A, χ0B, b, using_ηs) = initialize_integration_ranges(ηcA_or_τcA, ηcB_or_τcB, χ0A, χ0B, b, using_ηs, space_time, switching_func_name)

function initialize_trajs(space_time, χ0A, χ0B, b)
  if     space_time=="quench"   XA, XB = QuenchTrajectory(χ0A, b)             , QuenchTrajectory(χ0B, b)
  elseif space_time=="flat_2d"  XA, XB = InertialTrajectory(χ0A, 0, 0)        , InertialTrajectory(χ0B, 0, 0) 
  elseif space_time=="rindler"  XA, XB = AcceleratedTrajectory2D(χ0A)         , AcceleratedTrajectory2D(χ0B) 
  elseif space_time=="rindler2" XA, XB = AcceleratedTrajectory(χ0A, χ0A, 0, 0), AcceleratedTrajectory(χ0B, χ0B, 0, 0) 
  else   println("I don't know what trajectory to analize given the spacetime under study: space_time = $space_time") end
  return XA, XB
end

function initialize_Ws(W_func, XA, XB, with_derivative_coupling, ε_numeric_derivative)
  Ws = Dict()
  Ws["AB"] = DistributionWithTrajectories(W_func, XB, XA)
  Ws["AA"] = DistributionWithTrajectories(W_func, XA, XA)
  Ws["BB"] = DistributionWithTrajectories(W_func, XB, XB)
  if with_derivative_coupling 
    _gcd(f) = get_crossed_derivative(f, ε_numeric_derivative)
    Ws = map_dict(_gcd, Ws)
  end
  return Ws
end

function initialize_D(W, XA, XB, with_derivative_coupling, ε_numeric_derivative)
  WAB = DwT(W, XA, XB)
  WBA = DwT(W, XB, XA)
  if with_derivative_coupling 
    _gcd(f) = get_crossed_derivative(f, ε_numeric_derivative)
    WAB, WBA = map(_gcd, [WAB, WBA])
  end
  function D(τA, τB)
    if     XA(τA)[1]  > XB(τB)[1] return WAB(τA, τB) 
    elseif XA(τA)[1] <= XB(τB)[1] return WBA(τB, τA) end
  end
  return D
end

function initialize_distributions(W::Function, XA, XB, with_derivative_coupling, ε_numeric_derivative) 
   return (initialize_Ws(W, XA, XB, with_derivative_coupling, ε_numeric_derivative), 
           initialize_D( W, XA, XB, with_derivative_coupling, ε_numeric_derivative))
end

initialize_distributions(W, XA, XB) = initialize_distributions(W, XA, XB, false, 0.0) 

function integrate_ls(ls, initial_τs, final_τs, integrate_func)
  τiA, τiB = initial_τs
  τfA, τfB = final_τs
  Ls = Dict()
  Ls["AA"] = integrate_func(ls["AA"], [τiA, τiA], [τfA, τfA])
  Ls["BB"] = integrate_func(ls["BB"], [τiB, τiB], [τfB, τfB])
  if "AB" in keys(ls) Ls["AB"] = integrate_func(ls["AB"], [τiA, τiB], [τfA, τfB]) end
  return Ls
end

function integrate_m_and_ls(m, ls, initial_τs, final_τs, integrate_func) 
  return integrate_func(m, initial_τs, final_τs), integrate_ls(ls, initial_τs, final_τs, integrate_func)
end

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

function get_ρ(M, Ls)
   #                         00        01              10        11
  ρ = [1 - Ls["AA"] - Ls["BB"]                  0                        0   conj(M); #00
                             0            Ls["BB"]  conj(get(Ls, "AB", NaN))       0; #01
                             0  get(Ls, "AB", NaN)                 Ls["AA"]        0; #10
                             M                  0                        0         0] #11

  if imag(ρ[2,2]) >   1e-5 println("the imaginary part of PB is quite big so there could be a numerical problem: PB = $(ρ[2,2])\n") end  
  if imag(ρ[3,3]) >   1e-5 println("the imaginary part of PA is quite big so there could be a numerical problem: PB = $(ρ[3,3])\n") end  
  if real(ρ[2,2]) <= -1e-3 println("PB is quite negative so there could be a numerical problem: PB = $(ρ[2,2])\n") end
  if real(ρ[3,3]) <= -1e-3 println("PA is quite negative so there could be a numerical problem: PA = $(ρ[3,3])\n") end
  if real(ρ[2,2]) >   1    println("PB is too big so there could be a numerical problem: PB = $(ρ[2,2])\n") end
  if real(ρ[3,3]) >   1    println("PA is too big so there could be a numerical problem: PA = $(ρ[3,3])\n") end
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

function make_img(χ0Bs, Ωs, Cs, path, experiment_name, save_bool)
  f, ax, hm = CairoMakie.contourf(χ0Bs, Ωs, Cs', linewidth=-0.0)
  ax.xlabel = "χB₀"; ax.ylabel = "Ω" #; ax.title = "Concurrence"
  Colorbar(f[:, end+1], hm)
  display(f)
  img_name = replace("$(experiment_name)_$(now())", ":"=>"_")
<<<<<<< HEAD
  if save_bool save("$path/$(img_name).svg", f) end
=======
  if save_bool save("$path/$(img_name).pdf", f) end
>>>>>>> 42ee51a3dd6d56ad0323f7282165f99107e7c0c5
  return img_name
end

function plot_C_vs_L(path, experiment_name, ΔLss, Ωs, Cs, save_img=true)
  img_names = []
  for Ω in Ωs
    f, ax, l = CairoMakie.lines(ΔLss[Ω], vec(Cs[Ω]/λ^2), fontsize = 12)
    CairoMakie.ylims!(ax, C_ranges[Ω])
    ax.title = "Ω = $Ω"
    ax.titlesize = 25
    display(f)
    img_name = replace("$(experiment_name)_$(now())", ":"=>"_")
    if save_img CairoMakie.save("$path/$(img_name).pdf", f, pt_per_unit = 1) end
    push!(img_names, img_name)
  end
  return img_names
end

function store_in_df(path, file_name, params, img_names, run_durations)
  for (img_name, run_duration) in zip(img_names, run_durations)
    if hasproperty(params, :Image_Name) params[:,"Image_Name"] = [img_name]
    else
      insertcols!(params, 1, "Image_Name"   => [img_name],
                             "Run_Duration" => [Int.(round.(run_duration./60))])
    end
    append = file_name in readdir(path)
    CSV.write("$path/$file_name", params, append=append)
  end
  df = CSV.read( "$path/$file_name", DataFrame)
  sort!(df, [:Image_Name])
  CSV.write("$path/$file_name", df)
end