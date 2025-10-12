using CairoMakie
using LinearAlgebra
include("utils.jl")
include("dynamics.jl")
include("visualize.jl")

## =================== G-C base pair ===================

a_can, b_can, c_can = 0.006457467585167605, 0.03131130708309966, 0.03979932858576989
a_bar, b_bar, c_bar = -0.006425438605566449, 0.0038910266591298437, 0.025223914926830745
a_tau, b_tau, c_tau = 0.013406834311699567, -0.044575846347551726, 0.05473263795143025

canonical, barrier, tautomerical = create_params_struct(
    a_can, b_can, c_can,
    a_bar, b_bar, c_bar,
    a_tau, b_tau, c_tau
)

x_c, x_t = -1.03, 1.15

# ONLY CANONICAL TIME DEPENDENCY [X_T AND X_C FIXED]
# can, bar, tau = evolve_canonical(canonical, barrier, tautomerical)
# make_gif_from_series(can, bar, tau, "graphics/model/parabolas_gc.gif")

# CANONICAL, BARRIER AND TAUTOMERICAL TIME DEPENDENCY [X_T AND X_C MOVING]
L_series = readdlm("data/L_series_gc.txt") |> vec
R_series = readdlm("data/R_series_gc.txt") |> vec
can, bar, tau, x_c_series, x_t_series = evolve_all_forms(canonical, barrier, tautomerical)
make_gif_independent_series(can, bar, tau, x_c_series, x_t_series, L_series, R_series, "graphics/model/parabolas_independent_gc.gif")

## =================== A-T base pair ==================
a_can, b_can, c_can = 0.01548757014342916, 0.0544195348802887, 0.04817678375503837
a_bar, b_bar, c_bar = -0.010377758242695038, 0.007613726939775605, 0.0310333936186076
a_tau, b_tau, c_tau = 0.0125386460155263, -0.04565578211748337, 0.0619375921034291

canonical, barrier, tautomerical = create_params_struct(
    a_can, b_can, c_can,
    a_bar, b_bar, c_bar,
    a_tau, b_tau, c_tau
)
x_c, x_t = -0.51, 1.21

# ONLY CANONICAL TIME DEPENDENCY [X_T AND X_C FIXED]
can, bar, tau = evolve_canonical(canonical, barrier, tautomerical, false)
make_gif_from_series(can, bar, tau, "graphics/model/parabolas_at.gif", false)

# CANONICAL, BARRIER AND TAUTOMERICAL TIME DEPENDENCY [X_T AND X_C MOVING]
L_series = readdlm("data/L_series_at.txt") |> vec
R_series = readdlm("data/R_series_at.txt") |> vec

can, bar, tau, x_c_series, x_t_series = evolve_all_forms(canonical, barrier, tautomerical, false)
make_gif_independent_series(can, bar, tau, x_c_series, x_t_series, L_series, R_series, "graphics/model/parabolas_independent_at.gif", is_gc_base_pair=false)

## hermite approximation
include("hermite.jl")

# A-T + G-C
plot_analytical_expansion()
plot_hermite_expansion(false, 10)

# calculate the proximity of the hermite polynomial expansion and the Slocombe2022 and Godbeer2015
println(are_functions_close(false, 10))
println(are_functions_close())

## calculate the differences betwween energy eigenvalues of harmonic model and extrapolated functions 
is_at = true
x_range = is_at ? LinRange(-3.0, 3.0, 200) : LinRange(-4.0, 2.9, 200)
_, _, ene_left, ene_right = get_wavefunctions_qho(x_range, 14, 6)

include("../true_calc_at/schrodinger_eqn.jl")
"""
g-c: canonical eigenstates: 1-6, 7-11 - naprzemiennie
a-t: canonical eigenstates: 1-5, 6-11 - naprzemiennie
"""

energies, _, _ = solve_schrodinger(13)
println("A-T:energy eigenstates differences")
for i in 1:11
    if i < 6
        println("(canonical diff) n = ", i, "\t", energies[i] - ene_left[i])
    elseif i > 5 && i % 2 == 0 
        println("(canonical diff) n = ", i, "\t", energies[i] - ene_left[i])
    elseif i > 5 && i % 2 != 0
        println("(tautomerical diff) n = ", i, "\t", energies[i] - ene_right[i-5])
    end
end

## g-c
is_at = true
x_range = is_at ? LinRange(-3.0, 3.0, 200) : LinRange(-4.0, 2.9, 200)
_, _, ene_left, ene_right = get_wavefunctions_qho(x_range, 14, 6)

energies, _, _ = solve_schrodinger(12, 1000, (-4., 2.9), false)
println("G-C:energy eigenstates differences")
for i in 1:12
    if i < 7
        println("(canonical diff) n = ", i, "\t", energies[i] - ene_left[i])
    elseif i > 6 && i % 2 == 0 
        println("(canonical diff) n = ", i, "\t", energies[i] - ene_left[i])
    elseif i > 6 && i % 2 != 0
        println("(tautomerical diff) n = ", i, "\t", energies[i] - ene_right[i-5])
    end
end