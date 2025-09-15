# using LinearAlgebra
using Plots
include("utils.jl")

function test_continuity(parabola1::ParabolaParams, parabola2::ParabolaParams, x_c::Float64)
    output = 1/x_c^2 * (parabola1.a * x_c ^ 2 + x_c * (2 * parabola2.a * parabola2.p - 2 * parabola1.a * parabola1.p) + parabola1.q - parabola2.q + parabola1.a * parabola1.p^2-parabola2.a * parabola2.p^2)
    return println("continuity test: ", output)
end


function barrier_coefficient(parabola1::ParabolaParams, parabola2::ParabolaParams, x_c::Float64)
    output = 1/x_c^2 * (parabola1.a * x_c ^ 2 + x_c * (2 * parabola2.a * parabola2.p - 2 * parabola1.a * parabola1.p) + parabola1.q - parabola2.q + parabola1.a * parabola1.p^2-parabola2.a * parabola2.p^2)
    return output
end


function tautomerical_coefficient(parabola1::ParabolaParams, parabola2::ParabolaParams, x_c::Float64)
    output = 1/x_c^2 * (parabola1.a * x_c ^ 2 + x_c * (2 * parabola2.a * parabola2.p - 2 * parabola1.a * parabola1.p) + parabola1.q - parabola2.q + parabola1.a * parabola1.p^2-parabola2.a * parabola2.p^2)
    return output
end

a_can_series = Float64[]
a_bar_series = Float64[]
a_tau_series = Float64[]

function plot_at_instance(canonical::ParabolaParams, barrier::ParabolaParams, tautomerical::ParabolaParams)
    xs_can = -4.5:0.01:x_c
    xs_bar = x_c:0.01:x_t
    xs_tau = x_t:0.01:2.7

    ys_can = [canonical.a*(x-canonical.p)^2 + canonical.q for x in xs_can]

    ys_bar = [barrier.a*(x-barrier.p)^2 + barrier.q for x in xs_bar]

    ys_tau = [tautomerical.a*(x-tautomerical.p)^2 + tautomerical.q for x in xs_tau]

    lines(xs_can, ys_can)
    lines!(xs_bar, ys_bar)
    lines!(xs_tau, ys_tau)
    display(current_figure())
end

function update_params(canonical::ParabolaParams, barrier::ParabolaParams, tautomerical::ParabolaParams)
    T = 10.
    dt = 0.1
    times = 0:dt:T
    for t in times
        canonical.a += .00005*sin(5*t) 

        canonical.p, canonical.q = get_pq_params(
            GeneralParabolaParams(canonical.a, canonical.b, canonical.c)
        )

        barrier.a = barrier_coefficient(canonical, barrier, x_c)
        y_c = canonical.a * (x_c - canonical.p)^2 + canonical.q
        barrier.q = y_c - barrier.a * (x_c - barrier.p)^2

        tautomerical.a = tautomerical_coefficient(barrier, tautomerical, x_t)
        y_t = barrier.a * (x_t - barrier.p)^2 + barrier.q
        tautomerical.q = y_t - tautomerical.a * (x_t - tautomerical.p)^2

        push!(a_can_series, canonical.a)
        push!(a_bar_series, barrier.a)
        push!(a_tau_series, tautomerical.a)

        y_can = canonical.a * (x_c - canonical.p)^2 + canonical.q
        y_bar = barrier.a * (x_c - barrier.p)^2 + barrier.q
        y_bar_2 = barrier.a * (x_t - barrier.p)^2 + barrier.q
        y_tau = tautomerical.a * (x_t - tautomerical.p)^2 + tautomerical.q   

        println(abs(y_can - y_bar)) 
        println(abs(y_tau - y_bar_2))


    end
end


function make_gif_from_series(canonical::ParabolaParams, barrier::ParabolaParams, tautomerical::ParabolaParams, filename::String="parabolas.gif")
    xs_can = collect(-4.5:0.01:x_c)
    xs_bar = collect(x_c:0.01:x_t)
    xs_tau = collect(x_t:0.01:2.7)

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

