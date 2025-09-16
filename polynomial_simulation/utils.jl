mutable struct GeneralParabolaParams
    a::Float64
    b::Float64
    c::Float64
end

mutable struct ParabolaParams
    a::Float64
    b::Float64
    c::Float64
    p::Float64
    q::Float64
end


function get_pq_params(gen_params::GeneralParabolaParams)
    a = gen_params.a
    b = gen_params.b
    c = gen_params.c
    p = -b/(2a)
    q = c - (b^2)/(4a)
    return (p, q)
end

function get_params(gen_params::GeneralParabolaParams)
    a = gen_params.a
    b = gen_params.b
    c = gen_params.c
    (p, q) = get_pq_params(gen_params)
    return (a, b, c, p, q)
end

function create_params_struct(a_can::Float64, b_can::Float64, c_can::Float64,
                              a_bar::Float64, b_bar::Float64, c_bar::Float64,
                              a_tau::Float64, b_tau::Float64, c_tau::Float64)
    a, b, c, p, q = get_params(GeneralParabolaParams(a_can, b_can, c_can))
    canonical = ParabolaParams(a, b, c, p, q)

    a, b, c, p, q = get_params(GeneralParabolaParams(a_bar, b_bar, c_bar))
    barrier = ParabolaParams(a, b, c, p, q)

    a, b, c, p, q = get_params(GeneralParabolaParams(a_tau, b_tau, c_tau))
    tautomerical = ParabolaParams(a, b, c, p, q)
    return canonical, barrier, tautomerical
    
end

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

#=
y_can = canonical.a * (x_c - canonical.p)^2 + canonical.q
y_bar = barrier.a * (x_c - barrier.p)^2 + barrier.q
y_tau = tautomerical.a * (x_t - tautomerical.p)^2 + tautomerical.q        
println(abs(y_can - y_bar)) 
println(abs(y_tau - y_bar))

#
# testing
println("barrier test: ", barrier_coefficient(canonical, barrier, x_c), a_bar)
println("a_bar = ", a_bar)
test_continuity(canonical, barrier, x_c)
println(" ======")
println("a_tau = ", a_tau)
test_continuity(barrier, tautomerical, x_t)
#
=#