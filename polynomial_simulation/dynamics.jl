include("utils.jl")


a_can_series = Float64[]
a_bar_series = Float64[]
a_tau_series = Float64[]


function update_params(canonical::ParabolaParams, barrier::ParabolaParams, tautomerical::ParabolaParams)
    T = 10.
    dt = 0.1
    times = 0:dt:T
    canonical = deepcopy(canonical)
    barrier = deepcopy(barrier)
    tautomerical = deepcopy(tautomerical)
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
    return canonical, barrier, tautomerical
end


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

function update_params_independent(canonical::ParabolaParams, barrier::ParabolaParams, 
    tautomerical::ParabolaParams, is_gc_base_pair::Bool = true)
    T = 10.0
    dt = 0.1
    times = 0:dt:T

    empty!(a_can_series)
    empty!(a_bar_series)
    empty!(a_tau_series)

    canonical = deepcopy(canonical)
    barrier = deepcopy(barrier)
    tautomerical = deepcopy(tautomerical)

    global x_c, x_t

    for t in times
        # IMPORTANT: changing the paramaters may alter the continuity and the ability of finding real x_c and x_t
        if is_gc_base_pair
            canonical.a    += .00005 * sin(5*t) 
            barrier.a      += .00006 * sin(4*t + .89)
            tautomerical.a += .00002 * sin(1.5*t + 2.0)
        else    
            canonical.a    += .00005 * sin(5*t) 
            barrier.a      += .00006 * sin(3*t + 0.2)
            tautomerical.a += .00002 * sin(1.5*t + 2.0)
        end

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
    return canonical, barrier, tautomerical
end


# "Canonical-driven evolution (barrier and tau depend on canonical)."
# function evolve_canonical(canonical::ParabolaParams, barrier::ParabolaParams, tau::ParabolaParams;
#                           T=10.0, dt=0.1, x_c=-1.0, x_t=1.0)
#     times = 0:dt:T
#     can_series, bar_series, tau_series = Float64[], Float64[], Float64[]

#     canonical = deepcopy(canonical)
#     barrier   = deepcopy(barrier)
#     tau       = deepcopy(tau)

#     for t in times
#         canonical.a += 0.00005 * sin(5t)
#         canonical.p, canonical.q = get_pq_params(GeneralParabolaParams(canonical.a, canonical.b, canonical.c))

#         barrier.a = continuity_error(canonical, barrier, x_c) # or barrier_coefficient if you keep it
#         y_c = canonical.a*(x_c-canonical.p)^2 + canonical.q
#         barrier.q = y_c - barrier.a*(x_c - barrier.p)^2

#         tau.a = continuity_error(barrier, tau, x_t) # or tau_coefficient
#         y_t = barrier.a*(x_t-barrier.p)^2 + barrier.q
#         tau.q = y_t - tau.a*(x_t - tau.p)^2

#         push!(can_series, canonical.a)
#         push!(bar_series, barrier.a)
#         push!(tau_series, tau.a)
#     end

#     return can_series, bar_series, tau_series
# end


# "Independent evolution: each parabola wiggles on its own."
# function evolve_independent(canonical::ParabolaParams, barrier::ParabolaParams, tau::ParabolaParams;
#                             T=10.0, dt=0.1, is_gc_base_pair=true)
#     times = 0:dt:T
#     can_series, bar_series, tau_series = Float64[], Float64[], Float64[]

#     canonical = deepcopy(canonical)
#     barrier   = deepcopy(barrier)
#     tau       = deepcopy(tau)

#     for t in times
#         if is_gc_base_pair
#             canonical.a    += 0.00005 * sin(5t)
#             barrier.a      += 0.00006 * sin(4t + 0.89)
#             tau.a          += 0.00002 * sin(1.5t + 2.0)
#         else
#             canonical.a    += 0.00005 * sin(5t)
#             barrier.a      += 0.00006 * sin(3t + 0.2)
#             tau.a          += 0.00002 * sin(1.5t + 2.0)
#         end

#         canonical.p, canonical.q = get_pq_params(GeneralParabolaParams(canonical.a, canonical.b, canonical.c))
#         barrier.p, barrier.q     = get_pq_params(GeneralParabolaParams(barrier.a, barrier.b, barrier.c))
#         tau.p, tau.q             = get_pq_params(GeneralParabolaParams(tau.a, tau.b, tau.c))

#         push!(can_series, canonical.a)
#         push!(bar_series, barrier.a)
#         push!(tau_series, tau.a)
#     end

#     return can_series, bar_series, tau_series
# end
