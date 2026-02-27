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
#=
0.006457467585167605 0.03131130708309966 0.03979932858576989
-0.006573285204480409 0.0038910266591298437 0.025023914926830745
0.01302505750247126 -0.044575846347551726 0.05473263795143025
=#
#=
0.014348757014342916 0.0544195348802887 0.05317678375503837
-0.009542828706624436 0.007613726939775605 0.0310333936186076
0.012870124673305842 -0.04565578211748337 0.0619375921034291
=#
##
function energy_differences()
    ene_at, _, _ = solve_schrodinger(14, 1000, (-3.5, 3.0), true)
    ene_gc, _, _ = solve_schrodinger(14, 1000, (-4.0, 2.9), false)

    ene_at_ho, _, _ = solve_schrodinger_sum_harmonic(14, 1000, (-3.5, 3.0), true)
    ene_gc_ho, _, _ = solve_schrodinger_sum_harmonic(14, 1000, (-4.0, 2.9), false)

    n_at = min(length(ene_at), length(ene_at_ho))
    n_gc = min(length(ene_gc), length(ene_gc_ho))

    ΔE_at = ene_at[1:n_at] .- ene_at_ho[1:n_at]
    ΔE_gc = ene_gc[1:n_gc] .- ene_gc_ho[1:n_gc]

    ns_at = 0:(n_at-1)
    ns_gc = 0:(n_gc-1)

    fig = Figure(resolution = (1000, 400))
    ax1 = Axis(fig[1,1], xlabel = L"$n$", ylabel = L"$\Delta E \; (10^{-15} \; \text{a.u.})$", title = L"\text{A-T}",
                xlabelsize = 28, ylabelsize = 28, titlesize = 28, xticklabelsize = 24, yticklabelsize = 24)
    ax2 = Axis(fig[1,2], xlabel = L"$n$", ylabel = L"$\Delta E \; (10^{-15} \; \text{a.u.})$", title = L"\text{G-C}",
                xlabelsize = 28, ylabelsize = 28, titlesize = 28, xticklabelsize = 24, yticklabelsize = 24)
    ax1.yticks = [0, 1.0, 2.0, 3.0]
    ax2.yticks = [0, 1, 2, 3]
    scatter!(ax1, ns_at, abs.(ΔE_at).*1e15, marker = :circle, color = :violetred4, label = "ΔE (exact - HO)")
    scatter!(ax2, ns_gc, abs.(ΔE_gc)*1e15, marker = :circle, color = :palegreen4, label = "ΔE (exact - HO)")

    save("graphics/energy_diff.pdf", fig)
    return fig
end

function get_coth_approx()
    kb = 3.167e-6
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1,1], 
        xlabel = L"$x = {2k_B T}/{(\hbar\omega)}$", 
        ylabel = L"$\varepsilon / (\hbar\omega/2)$",
        limits = ((0.0, 1.5), (-5, 3)),
    xlabelsize = 30, ylabelsize = 30, xticklabelsize = 25, yticklabelsize = 25)
    x = LinRange(0.0001, 2, 500)

    x_fill = LinRange(0, 2*300.0*kb/(4.107e-03), 500)
    lines!(ax, fill(2*300.0*kb/(4.107e-03), 500), LinRange(-10, 3, 500), color=:gray, linestyle=:dashdot, linewidth=5)
    
    approx1(x) = x 
    approx2(x) = x + 1 / (3 * x)
    approx3(x) = x + 1 / (3 * x) - 1 / (45 * x^3)
    ax.yticks = [-4, -2, 0, 2]
    text!(ax, 2*300.0*kb/(4.107e-03)+0.03, -1, text = L"T = 300\text{ K}", fontsize = 30)
    
    coth_all(x) = (exp(1/x) + exp(-1/x)) / (exp(1/x) - exp(-1/x))

    lines!(ax, x, approx1, label = L"I: x", color = :red, linewidth = 5)
    lines!(ax, x, approx2, label = L"II: x + \frac{1}{3x}\quad ", color = :blue, linewidth = 5)
    lines!(ax, x, approx3, label = L"III: x + \frac{1}{3x} - \frac{1}{45x^3}", color = :green, linewidth = 5)
    lines!(ax, x, coth_all, label = L"\mathrm{ctgh}(1/x)", color = :black, linestyle = :dash, linewidth = 5)
    
    axislegend(ax, position = :rb, labelsize = 30, framevisible = false)
    save("graphics/ctgh_x.pdf", fig)
    return fig
end