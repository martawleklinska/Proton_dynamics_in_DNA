include("calculate.jl")
include("visualize.jl")

## =================== A-T base pair ==================

x_vals = LinRange(-3.5, 3, 100)
x_vals = collect(x_vals)
p_vals = LinRange(-12, 12, 100)
p_vals = collect(p_vals)

L_series = readdlm("data/L_series_at.txt") |> vec
R_series = readdlm("data/R_series_at.txt") |> vec
WDF_series_at = compute_WDF_for_time_varying_potential(x_vals, p_vals, L_series, R_series, 7, 1, title = "at")

x_vals = collect(LinRange(-3.5, 3, 100))
p_vals = collect(LinRange(-12, 12, 100))
make_wdf_gif_from_series(x_vals, p_vals, WDF_series_at, "graphics/wdf_at.gif", false)

## =================== G-C base pair ==================
x_vals = LinRange(-4.5, 2.7, 100)
x_vals = collect(x_vals)
p_vals = LinRange(-12, 12, 100)
p_vals = collect(p_vals)

L_series = readdlm("data/L_series_gc.txt") |> vec
R_series = readdlm("data/R_series_gc.txt") |> vec

x_vals = collect(LinRange(-4.5, 2.7, 100))
p_vals = collect(LinRange(-12, 12, 100))
WDF_series_gc = compute_WDF_for_time_varying_potential(x_vals, p_vals, L_series, R_series, 7, 1, title = "gc")

make_wdf_gif_from_series(x_vals, p_vals, WDF_series_gc, "graphics/wdf_gc.gif", true)
