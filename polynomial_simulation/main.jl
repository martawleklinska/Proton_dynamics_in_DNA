using CairoMakie
using LinearAlgebra
include("utils.jl")
include("dynamics.jl")
include("visualize.jl")

## =================== G-C base pair ===================

a_can, b_can, c_can = 0.006457467585167605, 0.03131130708309966, 0.03979932858576989
a_bar, b_bar, c_bar = -0.006425438605566449, 0.0038910266591298437, 0.025223914926830745
a_tau, b_tau, c_tau = 0.013406834311699567, -0.044575846347551726, 0.05473263795143025

#=
0.006457467585167605 0.03131130708309966 0.03979932858576989
-0.006573285204480409 0.0038910266591298437 0.025023914926830745
0.01302505750247126 -0.044575846347551726 0.05473263795143025
=#

canonical, barrier, tautomerical = create_params_struct(
    a_can, b_can, c_can,
    a_bar, b_bar, c_bar,
    a_tau, b_tau, c_tau
)

x_c, x_t = -1.03, 1.15

## ONLY CANONICAL TIME DEPENDENCY [X_T AND X_C FIXED]
can, bar, tau = update_params(canonical, barrier, tautomerical)
make_gif_from_series(can, bar, tau, "graphics/parabolas_gc.gif")


## CANONICAL, BARRIER AND TAUTOMERICAL TIME DEPENDENCY [X_T AND X_C MOVING]
can, bar, tau = update_params_independent(canonical, barrier, tautomerical)
make_gif_independent_series(can, bar, tau, "graphics/parabolas_independent_gc.gif")

## =================== A-T base pair ==================
a_can, b_can, c_can = 0.01548757014342916, 0.0544195348802887, 0.04817678375503837
a_bar, b_bar, c_bar = -0.010377758242695038, 0.007613726939775605, 0.0310333936186076
a_tau, b_tau, c_tau = 0.0125386460155263, -0.04565578211748337, 0.0619375921034291

#=
0.014348757014342916 0.0544195348802887 0.05317678375503837
-0.009542828706624436 0.007613726939775605 0.0310333936186076
0.012870124673305842 -0.04565578211748337 0.0619375921034291
=#

canonical, barrier, tautomerical = create_params_struct(
    a_can, b_can, c_can,
    a_bar, b_bar, c_bar,
    a_tau, b_tau, c_tau
)
x_c, x_t = -0.51, 1.21


## ONLY CANONICAL TIME DEPENDENCY [X_T AND X_C FIXED]

can, bar, tau = update_params(canonical, barrier, tautomerical)
make_gif_from_series(can, bar, tau, "graphics/parabolas_at.gif", false)

## CANONICAL, BARRIER AND TAUTOMERICAL TIME DEPENDENCY [X_T AND X_C MOVING]

can, bar, tau = update_params_independent(canonical, barrier, tautomerical, false)
make_gif_independent_series(can, bar, tau, "graphics/parabolas_independent_at.gif")