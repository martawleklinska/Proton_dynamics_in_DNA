include("utils.jl")
## Three Independent Parabola
function get_xc(parabola1::ParabolaParams, parabola2::ParabolaParams)
    A = parabola1.a - parabola2.a
    B = 2 * (parabola2.a * parabola2.p - parabola1.a * parabola1.p)
    C = parabola1.q - parabola2.q + parabola1.a * parabola1.p^2 - parabola2.a * parabola2.p^2
    D = B^2 - 4*A*C

    if D < 0
        return nothing   # no real intersection
    end

    x1 = (-B + sqrt(D)) / (2*A)
    x2 = (-B - sqrt(D)) / (2*A)

    # pick the physically relevant one 
    return x1 < x2 ? x1 : x2
end

function update_params_independent(canonical::ParabolaParams, barrier::ParabolaParams, tautomerical::ParabolaParams)
    T = 10.0
    dt = 0.1
    times = 0:dt:T

    empty!(a_can_series)
    empty!(a_bar_series)
    empty!(a_tau_series)

    global x_c, x_t

    for t in times
        # IMPORTANT: changing the paramaters may alter the continuity and the ability of finding real x_c and x_t
        canonical.a    += .00005 * sin(5*t) 
        barrier.a      += .00006 * sin(4*t + .89)
        tautomerical.a += .00002 * sin(1.5*t + 2.0)

        canonical.p, canonical.q = get_pq_params(
            GeneralParabolaParams(canonical.a, canonical.b, canonical.c)
        )
        barrier.p, barrier.q = get_pq_params(
            GeneralParabolaParams(barrier.a, barrier.b, barrier.c)
        )
        tautomerical.p, tautomerical.q = get_pq_params(
            GeneralParabolaParams(tautomerical.a, tautomerical.b, tautomerical.c)
        )

        # recompute moving junctions (cannot be imaginary)
        x_c = get_xc(canonical, barrier)
        x_t = get_xc(barrier, tautomerical)

        if x_c !== nothing
            y_can = canonical.a * (x_c - canonical.p)^2 + canonical.q
            y_bar = barrier.a * (x_c - barrier.p)^2 + barrier.q
            println("continuity error (can-bar): ", abs(y_can - y_bar))
        else
            println("No real intersection between canonical and barrier at t = $t")
        end

        if x_t !== nothing
            y_tau = tautomerical.a * (x_t - tautomerical.p)^2 + tautomerical.q
            y_bar_2 = barrier.a * (x_t - barrier.p)^2 + barrier.q
            println("continuity error (bar-tau): ", abs(y_tau - y_bar_2))
        else
            println("No real intersection between barrier and tautomerical at t = $t")
        end
        
        push!(a_can_series, canonical.a)
        push!(a_bar_series, barrier.a)
        push!(a_tau_series, tautomerical.a)
    end
end

function make_gif_independent_series(canonical::ParabolaParams, barrier::ParabolaParams, 
    tautomerical::ParabolaParams, filename::String="parabolas_independent.gif";
                                     nframes::Int=length(a_can_series))


    fig = Figure(size=(600,400))
    ax = Axis(fig[1,1]; xlabel=L"$x$", ylabel=L"$y$", title="Independent parabolas")
    CairoMakie.ylims!(ax, -0.001, 0.03)

    # for starters - the stationary values
    xs_can = Observable(LinRange(-4.5, 0.0, 100))
    xs_bar = Observable(LinRange(0.0, 1.0, 100))
    xs_tau = Observable(LinRange(1.0, 2.7, 100))

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

        xs_can[] = LinRange(-4.5, xc, 100)
        xs_bar[] = LinRange(xc, xt, 100)
        xs_tau[] = LinRange(xt, 2.7, 100)

        ys_can[] = [canonical.a*(x-canonical.p)^2 + canonical.q for x in xs_can[]]
        ys_bar[] = [barrier.a*(x-barrier.p)^2 + barrier.q for x in xs_bar[]]
        ys_tau[] = [tautomerical.a*(x-tautomerical.p)^2 + tautomerical.q for x in xs_tau[]]

        sc_c[] = Point2f(xc, y_c)
        sc_t[] = Point2f(xt, y_t)

    end
end

