using QuadGK, LinearAlgebra
using CairoMakie

const m_p = 1836.15267245  
const sqrtπ = sqrt(pi)
is_at = true
const dL = is_at ? -1.7737 : -2.42442
const dR = is_at ? 1.8963 : 1.787

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

xs = LinRange(-3.0, 3.0, 2001)
vals = [fourth_order(x) for x in xs]
i0 = argmin(vals)
x0 = xs[i0]
h = 1e-4
V0 = fourth_order(x0)
Vpp = (fourth_order(x0+h) - 2*V0 + fourth_order(x0-h)) / h^2   
a_local = Vpp / 2.0

c, α = local_hermite_coeffs(fourth_order, x0, a_local, nmax=6)
ys_approx = [local_hermite_approx(x, x0, α, c) for x in xs]
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y")
lines!(ax, xs, ys_approx)
fig