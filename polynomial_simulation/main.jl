using CairoMakie
using LinearAlgebra
include("utils.jl")
include("one_time_dependency.jl")

a_can, b_can, c_can = 0.006457467585167605, 0.03131130708309966, 0.03979932858576989
a_bar, b_bar, c_bar = -0.006425438605566449, 0.0038910266591298437, 0.025223914926830745
a_tau, b_tau, c_tau = 0.013406834311699567, -0.044575846347551726, 0.05473263795143025

#=
0.006457467585167605 0.03131130708309966 0.03979932858576989
-0.006573285204480409 0.0038910266591298437 0.025023914926830745
0.01302505750247126 -0.044575846347551726 0.05473263795143025
=#


a, b, c, p, q = get_params(GeneralParabolaParams(a_can, b_can, c_can))
canonical = ParabolaParams(a, b, c, p, q)

a, b, c, p, q = get_params(GeneralParabolaParams(a_bar, b_bar, c_bar))
barrier = ParabolaParams(a, b, c, p, q)

a, b, c, p, q = get_params(GeneralParabolaParams(a_tau, b_tau, c_tau))
tautomerical = ParabolaParams(a, b, c, p, q)

x_c, x_t = -1.03, 1.15

y_can = canonical.a * (x_c - canonical.p)^2 + canonical.q
y_bar = barrier.a * (x_c - barrier.p)^2 + barrier.q
y_tau = tautomerical.a * (x_t - tautomerical.p)^2 + tautomerical.q        
println(abs(y_can - y_bar)) 
println(abs(y_tau - y_bar))

#
# testing
# println("barrier test: ", barrier_coefficient(canonical, barrier, x_c), "\n", a_bar)
# println("a_bar = ", a_bar)
# test_continuity(canonical, barrier, x_c)
# println(" ======")
# println("a_tau = ", a_tau)
# test_continuity(barrier, tautomerical, x_t)
#
plot_at_instance(canonical, barrier, tautomerical)

## ONE TIME DEPENDENCY

update_params(canonical, barrier, tautomerical)
make_gif_from_series(canonical, barrier, tautomerical, "parabolas.gif")


## three independent parabolas
include("three_independent.jl")
x_c = get_xc(canonical, barrier)
x_t = get_xc(barrier, tautomerical)

update_params_independent(canonical, barrier, tautomerical)
make_gif_independent_series(canonical, barrier, tautomerical, "parabolas_independent.gif")