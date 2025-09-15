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


function test_continuity(parabola1::ParabolaParams, parabola2::ParabolaParams, x_c::Float64)
    output = 1/x_c^2 * (parabola1.a * x_c ^ 2 + x_c * (2 * parabola2.a * parabola2.p - 2 * parabola1.a * parabola1.p) + parabola1.q - parabola2.q + parabola1.a * parabola1.p^2-parabola2.a * parabola2.p^2)
    return println("continuity test: ", output)
end
