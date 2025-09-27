include("calculate.jl")
using CairoMakie
using Plots
using LaTeXStrings
using DelimitedFiles
using JLD2

##
a_left, a_right = 0.01548757014342916, 0.0125386460155263
dL, dR = -1.968, 1.968
x_vals = range(-3.5, 3, length=1000)
p_vals = range(-12, 12, length=1000)
R = sqrt(2 * a_right / 1836) * 1836
L = sqrt(2 * a_left / 1836) * 1836
function plot_WDF_map(x_vals, p_vals, Z)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"$x$", ylabel = L"$p$",
        title = L"\mathrm{WDF\; for\; A-T\; base-pair}", ylabelsize = 30,
        xlabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20,)
    hm = CairoMakie.heatmap!(ax, x_vals, p_vals, Z, colormap = :seismic, colorrange = (-0.3, 0.3))
    CairoMakie.Colorbar(fig[1, 2], hm, label = L"\varrho(x,\; p; \; t=0)", labelsize = 30, ticklabelsize = 20)
    display(fig)
end

# Z = compute_WDF_map(x_vals, p_vals, dL, dR, L, R, n=9, k=1)
# plot_WDF_map(x_vals, p_vals, Z)

## getting gif for the time-varying potential

using JLD2
@load "data/WDF_series_at.jld2" WDF_series
WDF_at = WDF_series

@load "data/WDF_series_gc.jld2" WDF_series
WDF_gc = WDF_series


function make_wdf_gif_from_series(x_vals, p_vals, WDF_series, filename::String="wdf_at.gif", is_gc_base_pair::Bool = false)
    nframes = size(WDF_series, 3)
    data = Observable(WDF_series[:, :, 1])
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],
        xlabel = L"$x$", ylabel = L"$p$",
        title = is_gc_base_pair ? L"\mathrm{WDF\; for\; G-C\; base-pair}" : L"\mathrm{WDF\; for\; A-T\; base-pair}",
        ylabelsize = 30, xlabelsize = 30, titlesize = 30,
        xticklabelsize = 20, yticklabelsize = 20)
    hm = CairoMakie.heatmap!(ax, x_vals, p_vals, data, colormap = :seismic, colorrange = (-0.3, 0.3))
    CairoMakie.Colorbar(fig[1, 2], hm, label = L"\varrho(x,\; p; \; t)", labelsize = 30, ticklabelsize = 20)

    record(fig, filename, 1:nframes; framerate = 20) do i
        data[] = WDF_series[:, :, i]
    end
end

x_vals = collect(LinRange(-3.5, 3, 100))
p_vals = collect(LinRange(-12, 12, 100))

make_wdf_gif_from_series(x_vals, p_vals, WDF_at, "graphics/wdf_at.gif", false)

