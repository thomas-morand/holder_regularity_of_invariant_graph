"The four files: structure.jl, upper_bound.jl, lower_bound.jl, and holder_regularity.jl
allow us to compute the Hölder regularity of the probability of extinction of a Galton-Watson process
in dynamical environnments in the uniformly supercritical case.
To do this, we need to estimate two Lyapunov exponents lambda_u and lambda_f. 
We estimate a lower bound of these Lyapunov exponents in the file lower_bound.jl
and a upper bound in the file upper_bound.jl.
The file holder_regularity.jl is used to link the different files 
and plot the Holder regularity of the probability of extinction of a family of Galton-Watson process
in dynamical environnments parametrised by a real parameter.
Here is the article associated with these programmes: ##############################
"

"This file contains the structures used in the code.
The first structure contains the necessary information about the transformation,
the second about the laws of reproduction.
Finally, a mutable structure allows us to manage intervals for interval arithmetic."

### Packages

using IntervalArithmetic
using Plots
using Dates
using StatsBase

### Transformations

"The transformation T is a C^1 uniformly dilating transformation of the circle seen as [0,1].
We also ask that T(0)=0 (but this is not limiting: if T(0)!=0,
simply apply the corresponding rotation to phi)."

struct T
    "The transformation T"
    value ::Function
    "The derivate of T"
    derivate ::Function
    "The maximum of the derivate of T"
    max_derivate ::Float64
    "The topological degree of T"
    degree ::Int64
end

### Examples of transformations

T1(n ::Int64, e ::Float64) = T(
    x -> n*x + e*sin(2*π*x),
    x -> n + 2*π*e*cos(2*π*x),
    n + 2*π*e,
    n
)
"n is an integer bigger than 2, and e is a float such that 
the absolute value of e is less than (n-1)/(π*2)"

### Laws of reproduction

"phi is a family of probability generating function.
The first variable x allow us to parametrize the family and take value in the circle.
So, s->phi(x,s) is the probability generating function of law mu_x"

struct phi
    "The law of reproduction"
    value ::Function
    "The logarithm of the derivate with respect to s of the law of reproduction"
    log_sderivate ::Function
end

### Examples of laws of reproduction

phi_poisson(l ::Float64) = phi(
    (x, s) -> exp((s-1) * exp(l + cos(2*x*π))),
    (x, s) -> log(exp(l + cos(2*x*π)) * phi_poisson(l).value(x, s))
)

phi_poisson_with_parameter(l ::Float64, a ::Float64)= phi(
    (x, s) -> phi_poisson(l).value(x + a, s),
    (x, s) -> phi_poisson(l).log_sderivate(x + a, s)
)

### Intervals

"The structure str_interval allows us to manage intervals in the context of interval arithmetic."

mutable struct str_interval
    "The interval"
    inter ::Union{Interval{Float64}, Interval{BigFloat}}
    "The integer N such that 2^(-N) is the length of the interval"
    length ::Int64
    "A boolean to know if the interval was obtained at the last stage"
    new ::Bool
    "The bound obtained on the interval (for the quantity we wish to bound on this interval)."
    bound ::Union{Float64, BigFloat}
end