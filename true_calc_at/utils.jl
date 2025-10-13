using LinearAlgebra
using CairoMakie
using QuadGK
using LaTeXStrings

function fourth_order(v)
    alpha = 1.963
    x = v/alpha
    a4 = 0.0207
    a3 = -0.0053
    a2 = -0.0414
    a1 = 0.0158
    a0 = 0.0312
    return a4*x.^4 + a3*x.^3 + a2*x.^2 + a1*x + a0
end

function get_coefficients_at()
    alpha = 1.963
    a4 = 0.0207
    a3 = -0.0053
    a2 = -0.0414
    a1 = 0.0158
    a0 = 0.0312
    return [a0, a1/alpha, a2/alpha^2, a3/alpha^3, a4/alpha^4]
end

function morse(V1, V2, a1, a2, q, r1, r2)
    exp1 = exp(-2 * a1 * (q - r1))
    exp2 = exp(-a1 * (q - r1))
    exp3 = exp(-2 * a2 * (r2 - q))
    exp4 = exp(-a2 * (r2 - q))
    output = V1 * (exp1 - 2 * exp2) + V2 * (exp3 - 2 * exp4) + 0.166 + 0.00019
    return output
end
function double_morse(q)
    V1 = 0.1617
    V2 = 0.082
    a1 = 0.305
    a2 = 0.755
    r1 = -2.7
    r2 = 2.1
    output = morse(V1, V2, a1, a2, q, r1, r2)
    return output
end

function eval_hermite(n::Int, x::Float64)
    if n == 0
        return one(x)
    elseif n == 1
        return 2 * x
    else
        l0 = one(x)       
        l1 = 2 * x        
        for k in 1:(n-1)
            l2 = (2 * x * l1 - 2 * k * l0)
            l0, l1 = l1, l2
        end
        return l1
    end
end

function get_integral_hermite(idx::Int, potential)
    integral, _ = quadgk(x -> exp(-x^2) * eval_hermite(idx, x) * potential(x), -10, 10, rtol=1e-8)
    return integral
end


function hermite_coefficients(nmax; is_at::Bool = true)
    potential = is_at ? (x -> fourth_order(x)) : (x -> double_morse(x))
    a_coefs = Float64[]
    for i in 0:nmax
        constan = 1/(sqrt(π) * 2^i * factorial(i))
        integral = get_integral_hermite(i, potential)
        push!(a_coefs, constan * integral)
    end
    return a_coefs
end


function hermite_approximation(x, a_coefs)
    s = zero(x)
    for (i, a) in enumerate(a_coefs)
        s += a * eval_hermite(i - 1, x)  
    end
    return s
end