include("calculate.jl")
using CairoMakie
using Plots
using LaTeXStrings
using DelimitedFiles
using JLD2

##
x_vals = range(-3.5, 3, length=1000)
p_vals = range(-12, 12, length=1000)

function plot_WDF_map(x_vals, p_vals, Z)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"$x$", ylabel = L"$p$",
        title = L"\mathrm{WDF\; for\; A-T\; base-pair}", ylabelsize = 30,
        xlabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20,)
    hm = CairoMakie.heatmap!(ax, x_vals, p_vals, Z, colormap = :seismic, colorrange = (-0.3, 0.3))
    CairoMakie.Colorbar(fig[1, 2], hm, label = L"\varrho(x,\; p; \; t=0)", labelsize = 30, ticklabelsize = 20)
    display(fig)
end

## getting gif for the time-varying potential


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

