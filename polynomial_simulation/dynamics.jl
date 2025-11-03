include("utils.jl")
using DelimitedFiles

a_can_series = Float64[]
a_bar_series = Float64[]
a_tau_series = Float64[]


function evolve_canonical(canonical::ParabolaParams, barrier::ParabolaParams, tautomerical::ParabolaParams,
                          is_gc_base_pair::Bool=true)
    T = 500.
    dt = 5.0
    times = 0:dt:T
    x_c = is_gc_base_pair ? -1.03 : -0.51
    x_t = is_gc_base_pair ? 1.15 : 1.21

    can_series, bar_series, tau_series = Float64[], Float64[], Float64[]

    canonical = deepcopy(canonical)
    barrier   = deepcopy(barrier)
    tautomerical = deepcopy(tautomerical)

    for t in times
        # frequency to equal thermal excitations: Ω ~ (9.8183 ± 0.1583) 10^{-4} a.u.    
        canonical.a += 0.000021349044066346842 * sin(9.8183e-04*t) 

        canonical.p, canonical.q = get_pq_params(
            GeneralParabolaParams(canonical.a, canonical.b, canonical.c)
        )

        barrier.a = barrier_coefficient(canonical, barrier, x_c)
        y_c = canonical.a * (x_c - canonical.p)^2 + canonical.q
        barrier.q = y_c - barrier.a * (x_c - barrier.p)^2

        tautomerical.a = tautomerical_coefficient(barrier, tautomerical, x_t)
        y_t = barrier.a * (x_t - barrier.p)^2 + barrier.q
        tautomerical.q = y_t - tautomerical.a * (x_t - tautomerical.p)^2

        push!(can_series, canonical.a)
        push!(bar_series, barrier.a)
        push!(tau_series, tautomerical.a)

    end
    return can_series, bar_series, tau_series
end


## Three Independent Parabola
function get_xc(parabola1::ParabolaParams, parabola2::ParabolaParams, is_gc_base_pair::Bool=true)
    max_iters = 1000
    iter = 0
    A = parabola1.a - parabola2.a
    B = 2 * (parabola2.a * parabola2.p - parabola1.a * parabola1.p)
    C = parabola1.q - parabola2.q + parabola1.a * parabola1.p^2 - parabola2.a * parabola2.p^2
    D = B^2 - 4*A*C
    while true
        A = parabola1.a - parabola2.a
        B = 2 * (parabola2.a * parabola2.p - parabola1.a * parabola1.p)
        C = parabola1.q - parabola2.q + parabola1.a * parabola1.p^2 - parabola2.a * parabola2.p^2
        D = B^2 - 4*A*C

        if D >= 0 || iter >= max_iters
            break
        end
        println(D)
        parabola1.a -= 0.00001 
        iter +=1

    end

    if iter == max_iters # still no real intersection
        return nothing
    end
    x1 = (-B + sqrt(D)) / (2*A)
    x2 = (-B - sqrt(D)) / (2*A)

    # pick the physically relevant one 
    if is_gc_base_pair
        return x1 > x2 ? x1 : x2
    else
        return x1 > x2 ? x1 : x2
    end
end

function evolve_all_forms(canonical_::ParabolaParams, barrier_::ParabolaParams, 
    tautomerical_::ParabolaParams, is_gc_base_pair::Bool = true)
    T = 20000
    dt = 100.
    times = 0:dt:T
    a_can_series, a_bar_series, a_tau_series, x_c_series, x_t_series = Float64[], Float64[], Float64[], Float64[], Float64[]
    L_series, R_series = Float64[], Float64[]
    
    canonical = deepcopy(canonical_)
    barrier = deepcopy(barrier_)
    tautomerical = deepcopy(tautomerical_)

    A = is_gc_base_pair ? 0.000067 : 0.0021 # limitation of the A coefficient (energy levels<>thermal energy) 
    A_tau = is_gc_base_pair ? 0.0018579056611541696 : 0.0017320978687266177 # same here
    m = 1836
    omega0 = sqrt(2 * canonical.a / m)
    omega0_barrier = sqrt(-2 * barrier.a /m)
    omega0_tautomer = sqrt(2 * tautomerical.a /m)
    # !!! frequency to equal thermal excitations: Ω ~ (9.8183 ± 0.1583) 10^{-4} a.u. !!! ============
    Ω = 9.8183e-04

    for t in times
        # IMPORTANT: changing the paramaters may alter the continuity and the ability of finding real x_c and x_t
        ωt = omega0 + A * sin(Ω * t)
        canonical.a    = 0.5 * m * ωt^2
        barrier.a      = -0.5 * m * (omega0_barrier - A * sin(Ω/3 * t+100))^2
        tautomerical.a = 0.5 * m * (omega0_tautomer + A_tau/50 * sin(Ω/3 * t+100))^2

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
        x_c = get_xc(canonical, barrier, is_gc_base_pair)
        x_t = get_xc(barrier, tautomerical, is_gc_base_pair)

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
        push!(x_c_series, x_c)
        push!(x_t_series, x_t)
        L = sqrt(abs(2 * canonical.a * 1836))
        R = sqrt(abs(2 * tautomerical.a * 1836))
        push!(L_series, L)
        push!(R_series, R)
    end

    if is_gc_base_pair
        writedlm("data/L_series_gc.txt", L_series)
        writedlm("data/R_series_gc.txt", R_series)
    else
        writedlm("data/L_series_at.txt", L_series)
        writedlm("data/R_series_at.txt", R_series)
    end

    return a_can_series, a_bar_series, a_tau_series, x_c_series, x_t_series, L_series, R_series
end
