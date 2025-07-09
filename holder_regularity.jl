"In this file, we compute and plot the Holder regularity of the probability of extinction.
We use the bound on the Lyapunov exponent lambda_u and lambda_f find thanks to the files
lower_bound.jl and upper_bound.jl."

### Files used 

include("structure.jl")
include("upper_bound.jl")
include("lower_bound.jl")

### Plot Holder regularity

function plot_bound_holder_reg(phi, T ::T, first_parameter ::Float64, size_step ::Float64, nb_step ::Int64, M_orbit ::Int64, iteration ::Int64, plt ::Bool, big_float = nothing ::Union{Nothing, Int64}) ::Tuple{Vector{Float64}, Vector{Float64}}
    "This function plot the Holder regularity of the probability of exctinction of the law of reproduction
    phi and the transformation T for parameter l = first_paramter + size_step * i for i in 0::nb_step-1. 
    To do this, we compute the periodic orbit of T of size smaller than M_orbit.
    To compute the upper bound with a big precision, we can use BigFloat instead of Float64 by assigning
    the integer corresponding to the number of bits of precision desired to the big_float variable."
    if big_float != nothing
        setprecision(big_float)
    end
    "step 1: periodic orbit"
    orbit, epsi = periodic_orbit(T, 40, M_orbit)
    "step 2: lambda_u"
    lower_lambda_u = maximum([log_derivate_T_lower(T, orbx, epsi) for orbx in orbit])
    upper_lambda_u = upper_bound_lambda_u(T, 6, iteration, big_float)
    "step 3: lambda_f"
    lower_lambda_f = zeros(Float64, nb_step)
    for p in 1:nb_step
        l = first_parameter + (p - 1) * size_step
        lower_lambda_f[p] = maximum([F_lower(phi(l), orbx, length(orbx), 100, epsi) for orbx in orbit])
    end
    upper_lambda_f = zeros(Float64, nb_step)
    for p in 1:nb_step
        l = first_parameter + (p - 1) * size_step
        upper_lambda_f[p] = min(upper_bound_lambda_f(phi(l), T, 6, iteration, 12, 0.99, 0.01, 16, big_float), 0)
    end
    "step 4: plot"
    lower_regularity = upper_lambda_f ./ (-upper_lambda_u)
    upper_regularity = lower_lambda_f ./ (-lower_lambda_u)
    if plt == true
        L = [first_parameter + i * size_step for i in 0:nb_step-1]
        p = plot(L.-size_step, upper_regularity, seriestype = :step, label="upper bound", size =(1800,1200),
        linewidth=1, legendfontsize=20, guidefontsize=22, tickfontsize=22)
        plot!(p, L, lower_regularity, color = :black, seriestype = :step, label="lower bound", linewidth=1, legendfontsize=20)
        t=Dates.format(now(), "yyyy-mm-dd_HH_MM-SS")
        savefig("$first_parameter-$size_step-$nb_step-$iteration-$t.png")
        display(p)
    end
   return upper_regularity, lower_regularity
end 