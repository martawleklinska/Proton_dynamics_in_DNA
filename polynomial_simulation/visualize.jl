using CairoMakie
include("utils.jl")


function plot_at_instance(canonical::ParabolaParams, barrier::ParabolaParams, tautomerical::ParabolaParams, is_gc_base_pair::Bool = true)
    if is_gc_base_pair
        xs_can = -4.5:0.01:x_c
        xs_bar = x_c:0.01:x_t
        xs_tau = x_t:0.01:2.7
    else
        xs_can = -3.7:0.01:x_c
        xs_bar = x_c:0.01:x_t
        xs_tau = x_t:0.01:4

    end

    ys_can = [canonical.a*(x-canonical.p)^2 + canonical.q for x in xs_can]

    ys_bar = [barrier.a*(x-barrier.p)^2 + barrier.q for x in xs_bar]

    ys_tau = [tautomerical.a*(x-tautomerical.p)^2 + tautomerical.q for x in xs_tau]
    fig = Figure(size=(600,400))
    ax = Axis(fig[1,1]; xlabel="x", ylabel="y", title="Parabolas at instance")
    CairoMakie.ylims!(ax, -0.002, 0.04)
    lines!(xs_can, ys_can)
    lines!(xs_bar, ys_bar)
    lines!(xs_tau, ys_tau)
    display(fig)
end



function make_gif_from_series(a_can_series, a_bar_series, a_tau_series, filename::String="parabolas.gif", is_gc_base_pair::Bool = true)
    nframes = length(a_can_series) 
    xs_can = is_gc_base_pair ? collect(-4.5:0.01:x_c) : collect(-3.7:0.01:x_c)
    xs_bar = collect(x_c:0.01:x_t)
    xs_tau = is_gc_base_pair ? collect(x_t:0.01:2.7) : collect(x_t:0.01:4.0) 
    
    fig = Figure(size=(600,400))
    ax = Axis(fig[1,1]; xlabel=L"$x$", ylabel=L"$y$", title=is_gc_base_pair ? "GC base pair" : "AT base pair")
    CairoMakie.ylims!(ax, -0.002, 0.04)

    ys_can = Observable(xs_can .* 0)
    ys_bar = Observable(xs_bar .* 0)
    ys_tau = Observable(xs_tau .* 0)

    lines!(ax, xs_can, ys_can, color=:blue)
    lines!(ax, xs_bar, ys_bar, color=:red)
    lines!(ax, xs_tau, ys_tau, color=:green)

    sc_c = Observable(Point2f(0.0, 0.0))
    sc_t = Observable(Point2f(0.0, 0.0))
    CairoMakie.scatter!(ax, sc_c, color=:black, markersize=12)
    CairoMakie.scatter!(ax, sc_t, color=:black, markersize=12)

    record(fig, filename, nframes; framerate=20) do i
        canonical.a = a_can_series[i]
        canonical.p, canonical.q = get_pq_params(
                GeneralParabolaParams(canonical.a, canonical.b, canonical.c)
            )

        barrier.a = a_bar_series[i]
        y_c = canonical.a * (x_c - canonical.p)^2 + canonical.q
        barrier.q = y_c - barrier.a * (x_c - barrier.p)^2

        tautomerical.a = a_tau_series[i]
        y_t = barrier.a * (x_t - barrier.p)^2 + barrier.q
        tautomerical.q = y_t - tautomerical.a * (x_t - tautomerical.p)^2

        ys_can[] = [canonical.a*(x-canonical.p)^2 + canonical.q for x in xs_can]
        ys_bar[] = [barrier.a*(x-barrier.p)^2 + barrier.q for x in xs_bar]
        ys_tau[] = [tautomerical.a*(x-tautomerical.p)^2 + tautomerical.q for x in xs_tau]
        
        sc_c[] = Point2f(x_c, y_c)
        sc_t[] = Point2f(x_t, y_t)
    end
end

## Three

function make_gif_independent_series(a_can_series, a_bar_series, a_tau_series, x_c_series, x_t_series, 
    L_series, R_series, filename::String="parabolas_independent.gif";
     is_gc_base_pair::Bool = true)
    energies = Observable(L_series[1])
    energies_r = Observable(R_series[1])

    nframes = length(a_can_series)
    title=is_gc_base_pair ? Observable("GC base pair") : Observable("AT base pair")

    fig = Figure(size=(600,400))
    ax = Axis(fig[1,1]; xlabel=L"$x$", ylabel=L"$y$", title=title)
    if is_gc_base_pair
        CairoMakie.ylims!(ax, -0.001, 0.03)

        # for starters - the stationary values
        xs_can = Observable(LinRange(-4.5, -1.03, 100))
        xs_bar = Observable(LinRange(-1.03, 1.15, 100))
        xs_tau = Observable(LinRange(1.15, 2.7, 100))
    else
        CairoMakie.ylims!(ax, -0.001, 0.04)
        xs_can = Observable(LinRange(-3.7, -0.51, 100))
        xs_bar = Observable(LinRange(-0.51, 1.21, 100))
        xs_tau = Observable(LinRange(1.21, 4., 100))

    end

    ys_can = Observable(zeros(length(xs_can[])))
    ys_bar = Observable(zeros(length(xs_bar[])))
    ys_tau = Observable(zeros(length(xs_tau[])))
    lines!(ax, xs_can, ys_can, color=:blue, label="canonical")
    lines!(ax, xs_bar, ys_bar, color=:red, label="barrier")
    lines!(ax, xs_tau, ys_tau, color=:green, label="tautomerical")


    sc_c = Observable(Point2f(0.0, 0.0))
    sc_t = Observable(Point2f(0.0, 0.0))
    CairoMakie.scatter!(ax, sc_c, color=:black, markersize=12)
    CairoMakie.scatter!(ax, sc_t, color=:black, markersize=12)
    CairoMakie.axislegend(ax; position = :rb, framevisible = false, textsize=15)

    # energies plot 
    canonical.a = a_can_series[1]
    tautomerical.a = a_tau_series[1]
    canonical.p, canonical.q = get_pq_params(
            GeneralParabolaParams(canonical.a, canonical.b, canonical.c)
        )
    tautomerical.p, tautomerical.q = get_pq_params(
            GeneralParabolaParams(tautomerical.a, tautomerical.b, tautomerical.c)
        )
    for i in 0:10
        y_level = lift(energies) do e
            e/1836 * (i + 0.5) + canonical.q
        end
        x = is_gc_base_pair ? [-2.7, -2.3] : [-1.95, -1.55]
        lines!(ax, x, lift(y_level) do y
            [y, y]
        end, color=:gray, linestyle=:dash)   
    end
    for i in 0:3
        y_level_r = lift(energies_r) do e
            e/1836 * (i + 0.5) + tautomerical.q
        end
        x = is_gc_base_pair ? [1.5, 1.9] : [1.65, 2.05]
        lines!(ax, x, lift(y_level_r) do y
            [y, y]
        end, color=:gray, linestyle=:dash)   
    end

    #

    record(fig, filename, 1:nframes; framerate = 20) do i
        canonical.a = a_can_series[i]
        barrier.a = a_bar_series[i]
        tautomerical.a = a_tau_series[i]

        energies[] = L_series[i]
        energies_r[] = R_series[i]
        title[] = is_gc_base_pair ? "GC base pair, t=$(i)" : "AT base pair, t=$(i)"

        # update params
        canonical.p, canonical.q = get_pq_params(
            GeneralParabolaParams(canonical.a, canonical.b, canonical.c)
        )
        barrier.p, barrier.q = get_pq_params(
            GeneralParabolaParams(barrier.a, barrier.b, barrier.c))
        tautomerical.p, tautomerical.q = get_pq_params(
            GeneralParabolaParams(tautomerical.a, tautomerical.b, tautomerical.c)
        )

        
        xc = x_c_series[i]
        xt = x_t_series[i]

        y_c = canonical.a * (xc - canonical.p)^2 + canonical.q
        barrier.q = y_c - barrier.a * (xc - barrier.p)^2
        y_t = barrier.a * (xt - barrier.p)^2 + barrier.q
        tautomerical.q = y_t - tautomerical.a * (xt - tautomerical.p)^2

        if is_gc_base_pair
            xs_can[] = LinRange(-4.5, xc, 100)
            xs_bar[] = LinRange(xc, xt, 100)
            xs_tau[] = LinRange(xt, 2.7, 100)
        else
            xs_can[] = LinRange(-3.7, xc, 100)
            xs_bar[] = LinRange(xc, xt, 100)
            xs_tau[] = LinRange(xt, 4.0, 100)
        end
        ys_can[] = [canonical.a*(x-canonical.p)^2 + canonical.q for x in xs_can[]]
        ys_bar[] = [barrier.a*(x-barrier.p)^2 + barrier.q for x in xs_bar[]]
        ys_tau[] = [tautomerical.a*(x-tautomerical.p)^2 + tautomerical.q for x in xs_tau[]]

        sc_c[] = Point2f(xc, y_c)
        sc_t[] = Point2f(xt, y_t)

    end
end
