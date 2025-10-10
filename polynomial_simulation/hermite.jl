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

##
function plot_hermite_expansion(is_at::Bool, nmax::Int)
    a_coefs = hermite_coefficients(nmax; is_at=is_at)

    potential = is_at ? (x -> fourth_order(x)) : (x -> double_morse(x))
    title_str = is_at ? L"\text{A-T}" : L"\text{G-C}"

    xs = is_at ? range(-3, 2.7, length=400) : range(-4.0, 2.7, length = 400)
    xlims = is_at ? (-3., 2.7) : (-4.0, 2.7)
    ylims = is_at ? (-0.005, 0.045) : (-0.002, 0.033)
    ys_true = [potential(x) for x in xs]
    ys_approx = [hermite_approximation(x, a_coefs) for x in xs]

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"$x$ (a.u.)", ylabel = L"$U(x)$ (a.u.)", 
            title = title_str, xlabelsize = 35, xticklabelsize = 25, yticklabelsize = 25,
            ylabelsize = 35, titlesize = 35, limits = (xlims, ylims))
    lines!(ax, xs, ys_true, color=:blue, label = is_at ? L"\text{potencjał}" : L"\text{potencjał}")
    lines!(ax, xs, ys_approx, color=:red, linestyle=:dash, label=L"Rozwinięcie $\mathcal{H}_n$ ($n=10$)")
    axislegend(ax, position=:rb, labelsize = 25)
    filename = is_at ? "graphics/model/hermite_expan_g-c.pdf" : "graphics/model/hermite_expan_a-t.pdf" 
    # save(filename, fig)
    fig
end

## analytical:A-T
function doublefactorial(number)
    fact = one(number)
    for m in iseven(number)+1:2:number
        fact *= m
    end
    return fact
end
function gaussian_moment(k::Int)
    if k % 2 != 0
        return 0.0
    elseif k == 0
        return sqrt(π)
    else
        return sqrt(π) * doublefactorial((k - 1)) / 2^(k/2)
    end
end

function analytic_hermite_coeffs_at(; nmax::Int=4)
    coeffs = get_coefficients_at()

    a = Float64[]
    for n in 0:nmax
        integral = 0.0
        for (k, c) in enumerate(coeffs)
            power = k - 1 
            if n == 0
                # H_0 = 1
                integral += c * gaussian_moment(power)
            elseif n == 1
                # H_1 = 2x
                integral += c * 2 * gaussian_moment(power + 1)
            elseif n == 2
                # H_2 = 4x^2 - 2
                integral += c * (4 * gaussian_moment(power + 2) - 2 * gaussian_moment(power))
            elseif n == 3
                # H_3 = 8x^3 - 12x
                integral += c * (8 * gaussian_moment(power + 3) - 12 * gaussian_moment(power + 1))
            elseif n == 4
                # H_4 = 16x^4 - 48x^2 + 12
                integral += c * (16 * gaussian_moment(power + 4) - 48 * gaussian_moment(power + 2) + 12 * gaussian_moment(power))
            end
        end
        coeff = integral / (sqrt(π) * 2^n * factorial(n))
        push!(a, coeff)
    end
    return a
end
function get_analytical_coefficients_gc(idx::Int = 10)
    V1, V2 = 0.1617, 0.082
    a1, a2 = 0.305, 0.755
    r1, r2 = -2.7, 2.1
    const_term = 0.166 + 0.00019

    betas = [-2*a1, -a1, 2*a2, a2]
    gammas = [2*a1*r1, a1*r1, -2*a2*r2, -a2*r2]
    A = [V1, -2*V1, V2, -2*V2] 

    a_hermite = Float64[]
    for n in 0:idx
        sum_val = 0.0
        for i in eachindex(betas)
            sum_val += A[i] * exp(gammas[i]) * (betas[i]^n) * exp(betas[i]^2 / 4)
        end
        push!(a_hermite, (sum_val / (2^n * factorial(n))) + const_term * (n == 0 ? 1 : 0))
    end
    return a_hermite
end

## plotting 

function plot_analytical_expansion(is_at::Bool = true)
    a_coefs = is_at ? analytic_hermite_coeffs_at() : get_analytical_coefficients_gc()

    potential = is_at ? (x -> fourth_order(x)) : (x -> double_morse(x))
    title_str = is_at ? L"\text{A-T}" : L"\text{G-C}"

    xs = is_at ? range(-3, 2.7, length=400) : range(-4.0, 2.7, length = 400)
    xlims = is_at ? (-3., 2.7) : (-4.0, 2.7)
    ylims = is_at ? (-0.005, 0.045) : (-0.002, 0.033)
    ys_true = [potential(x) for x in xs]
    ys_approx = [hermite_approximation(x, a_coefs) for x in xs]

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"$x$ (a.u.)", ylabel = L"$U(x)$ (a.u.)", 
            title = title_str, xlabelsize = 35, xticklabelsize = 25, yticklabelsize = 25,
            ylabelsize = 35, titlesize = 35, limits = (xlims, ylims))
    lines!(ax, xs, ys_true, color=:blue, linewidth = 3, label = is_at ? L"\text{potencjał}" : L"\text{potencjał}")
    lines!(ax, xs, ys_approx, color=:red, linestyle=:dash, linewidth = 5, label=L"Rozwinięcie $\mathcal{H}_n$ ($n=4$)")
    axislegend(ax, position=:rb, labelsize = 25)
    filename = is_at ? "graphics/model/hermite_expan_g-c.pdf" : "graphics/model/hermite_expan_a-t.pdf" 
    save(filename, fig) 
    # fig
end
## test optimization of max_idx
#(\int_{a}^b dx p(x) (f(x) - w(x))^2)^{1/2}
# \sum_{i=1}^N (f(x_i) - w(x_i))^2p(xi)
function are_functions_close(is_at::Bool = true, nmax = 4)
    # a_coefs = hermite_coefficients(nmax; is_at=is_at) # numerical
    a_coefs = is_at ? analytic_hermite_coeffs_at() : get_analytical_coefficients_gc(nmax)
    x_range = is_at ? LinRange(-3.0, 3.0, 100) : LinRange(-4.0, 2.7, 100)
    potential = is_at ? [fourth_order(x) for x in x_range] : [double_morse(x) for x in x_range]
    hermite = [hermite_approximation(x, a_coefs) for x in x_range] 
    proxi = 0
    for i in eachindex(x_range)
        proxi +=  (potential[i] - hermite[i]) ^ 2
    end
    return proxi
end

