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



function make_gif_from_series(canonical::ParabolaParams, barrier::ParabolaParams, 
    tautomerical::ParabolaParams, filename::String="parabolas.gif", is_gc_base_pair::Bool = true)
    if is_gc_base_pair
        xs_can = collect(-4.5:0.01:x_c)
        xs_bar = collect(x_c:0.01:x_t)
        xs_tau = collect(x_t:0.01:2.7)
    else
        xs_can = collect(-3.7:0.01:x_c)
        xs_bar = collect(x_c:0.01:x_t)
        xs_tau = collect(x_t:0.01:4)
    end    
    
    fig = Figure(size=(600,400))
    ax = Axis(fig[1,1]; xlabel="x", ylabel="y", title="Parabolas over time")
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

    record(fig, filename, 1:length(a_can_series); framerate=20) do i
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

function make_gif_independent_series(canonical::ParabolaParams, barrier::ParabolaParams, 
    tautomerical::ParabolaParams, filename::String="parabolas_independent.gif";
                                     nframes::Int=length(a_can_series), is_gc_base_pair::Bool = true)


    fig = Figure(size=(600,400))
    title = is_gc_base_pair ? "GC base pair" : "AT base pair"
    ax = Axis(fig[1,1]; xlabel=L"$x$", ylabel=L"$y$", title=title)
    if is_gc_base_pair
        CairoMakie.ylims!(ax, -0.001, 0.03)

        # for starters - the stationary values
        xs_can = Observable(LinRange(-4.5, 0.0, 100))
        xs_bar = Observable(LinRange(0.0, 1.0, 100))
        xs_tau = Observable(LinRange(1.0, 2.7, 100))
    else
        CairoMakie.ylims!(ax, -0.001, 0.06)
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

    record(fig, filename, 1:nframes; framerate = 20) do i
        canonical.a = a_can_series[i]
        barrier.a = a_bar_series[i]
        tautomerical.a = a_tau_series[i]

        # update params
        canonical.p, canonical.q = get_pq_params(
            GeneralParabolaParams(canonical.a, canonical.b, canonical.c)
        )
        barrier.p, barrier.q = get_pq_params(
            GeneralParabolaParams(barrier.a, barrier.b, barrier.c))
        tautomerical.p, tautomerical.q = get_pq_params(
            GeneralParabolaParams(tautomerical.a, tautomerical.b, tautomerical.c)
        )

        y_c = canonical.a * (x_c - canonical.p)^2 + canonical.q
        barrier.q = y_c - barrier.a * (x_c - barrier.p)^2
        y_t = barrier.a * (x_t - barrier.p)^2 + barrier.q
        tautomerical.q = y_t - tautomerical.a * (x_t - tautomerical.p)^2

        xc = get_xc(canonical, barrier)
        xt = get_xc(barrier, tautomerical)
        
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





# function make_gif_from_series(canonical::ParabolaParams, barrier::ParabolaParams, tau::ParabolaParams,
#                               can_series, bar_series, tau_series;
#                               filename="parabolas.gif", is_gc_base_pair=true,
#                               x_c=-1.0, x_t=1.0)

#     xs_can = is_gc_base_pair ? collect(-4.5:0.01:x_c) : collect(-3.7:0.01:x_c)
#     xs_bar = collect(x_c:0.01:x_t)
#     xs_tau = is_gc_base_pair ? collect(x_t:0.01:2.7) : collect(x_t:0.01:4.0)

#     fig = Figure(size=(600,400))
#     ax = Axis(fig[1,1]; xlabel="x", ylabel="y", title=is_gc_base_pair ? "GC base pair" : "AT base pair")
#     CairoMakie.ylims!(ax, -0.002, is_gc_base_pair ? 0.03 : 0.06)

#     ys_can = Observable(xs_can .* 0)
#     ys_bar = Observable(xs_bar .* 0)
#     ys_tau = Observable(xs_tau .* 0)

#     lines!(ax, xs_can, ys_can, color=:blue, label="canonical")
#     lines!(ax, xs_bar, ys_bar, color=:red, label="barrier")
#     lines!(ax, xs_tau, ys_tau, color=:green, label="tautomerical")

#     sc_c = Observable(Point2f(0.0, 0.0))
#     sc_t = Observable(Point2f(0.0, 0.0))
#     CairoMakie.scatter!(ax, sc_c, color=:black, markersize=12)
#     CairoMakie.scatter!(ax, sc_t, color=:black, markersize=12)

#     record(fig, filename, 1:length(can_series); framerate=20) do i
#         canonical.a = can_series[i]
#         barrier.a   = bar_series[i]
#         tau.a       = tau_series[i]

#         canonical.p, canonical.q = get_pq_params(GeneralParabolaParams(canonical.a, canonical.b, canonical.c))
#         barrier.p, barrier.q     = get_pq_params(GeneralParabolaParams(barrier.a, barrier.b, barrier.c))
#         tau.p, tau.q             = get_pq_params(GeneralParabolaParams(tau.a, tau.b, tau.c))

#         ys_can[] = [canonical.a*(x-canonical.p)^2 + canonical.q for x in xs_can]
#         ys_bar[] = [barrier.a*(x-barrier.p)^2 + barrier.q for x in xs_bar]
#         ys_tau[] = [tau.a*(x-tau.p)^2 + tau.q for x in xs_tau]
#     end
# end
