using CairoMakie
include("utils.jl")
include("dynamics.jl")
include("visualize.jl")

# === GC base pair ===
gc_can = ParabolaParams(get_params(GeneralParabolaParams(0.00645, 0.0313, 0.0398))...)
gc_bar = ParabolaParams(get_params(GeneralParabolaParams(-0.00642, 0.00389, 0.0252))...)
gc_tau = ParabolaParams(get_params(GeneralParabolaParams(0.0134, -0.0445, 0.0547))...)

can_series, bar_series, tau_series = evolve_independent(gc_can, gc_bar, gc_tau; is_gc_base_pair=true)
make_gif_from_series(gc_can, gc_bar, gc_tau, can_series, bar_series, tau_series; filename="graphics/parabolas_gc.gif", is_gc_base_pair=true)

# === AT base pair ===
at_can = ParabolaParams(get_params(GeneralParabolaParams(0.0154, 0.0544, 0.0481))...)
at_bar = ParabolaParams(get_params(GeneralParabolaParams(-0.0103, 0.00761, 0.0310))...)
at_tau = ParabolaParams(get_params(GeneralParabolaParams(0.0125, -0.0456, 0.0619))...)

can_series, bar_series, tau_series = evolve_independent(at_can, at_bar, at_tau; is_gc_base_pair=false)
make_gif_from_series(at_can, at_bar, at_tau, can_series, bar_series, tau_series; filename="graphics/parabolas_at.gif", is_gc_base_pair=false)
