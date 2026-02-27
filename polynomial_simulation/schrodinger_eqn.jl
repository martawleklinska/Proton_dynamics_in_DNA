using LinearAlgebra, SparseArrays
using Arpack
using CairoMakie
include("hermite.jl")

function model_at(x)
    a_can, b_can, c_can = 0.01548757014342916, 0.0544195348802887, 0.04817678375503837
    a_bar, b_bar, c_bar = -0.010377758242695038, 0.007613726939775605, 0.0310333936186076
    a_tau, b_tau, c_tau = 0.0125386460155263, -0.04565578211748337, 0.0619375921034291
    if x < -0.51
        return a_can .* x .^ 2 + b_can .* x + c_can
    elseif x < 1.21
        return a_bar .* x .^ 2 + b_bar .* x + c_bar
    else 
        return a_tau .* x .^ 2 + b_tau .* x + c_tau
    end
end

function model_gc(x)
    a_can, b_can, c_can = 0.006457467585167605, 0.03131130708309966, 0.03979932858576989
    a_bar, b_bar, c_bar = -0.006425438605566449, 0.0038910266591298437, 0.025223914926830745
    a_tau, b_tau, c_tau = 0.013406834311699567, -0.044575846347551726, 0.05473263795143025

    if x < -1.03
        return a_can .* x .^ 2 + b_can .* x + c_can
    elseif x < 1.15
        return a_bar .* x .^ 2 + b_bar .* x + c_bar
    else 
        return a_tau .* x .^ 2 + b_tau .* x + c_tau
    end
end
function get_potential_model(x; is_at=true)
    return is_at ? model_at.(x) : model_gc.(x)
end

""" 
TISE with FD returns: energies, wavefunctions, x
- takes:
    nstates: n/o of lowest eigenstates
    n: grid plot_at_instance
    xlims: (xmin, xmax)
- IMPORTANT
    This consideres the numerical solution of the whole system, so both wells 
    simultaneously. The analitical model however isn't at all solved nuemrically.
"""
function solve_schrodinger_sum_harmonic(nstates::Int64=10, n::Int64=1000, xlims::Tuple{Float64, Float64}=(-3.5, 3.0),
                            is_at::Bool=true)
    xmin, xmax = xlims
    x = range(xmin, stop = xmax, length = n)
    dx = step(x)
    nx = length(x)
    # potential matrix
    main = fill(1.0 / dx ^ 2 * 1/1836, nx)
    off = fill(-0.5 / dx ^ 2 * 1/1836, nx - 1)

    T = spdiagm(-1 => off, 0 => main, 1 => off)

    V_diag = get_potential_model(collect(x), is_at = is_at)
    V_diag .-= minimum(V_diag)  
    V = spdiagm(0 => V_diag)

    H = T + V

    nev = nstates
    vals, vecs = Arpack.eigs(H; nev=nev, which=:SR, maxiter=20000, tol=1e-22)

    idx = sortperm(real(vals))
    energies = real(vals[idx])
    wavefuncs = vecs[:, idx]
    
    # normalize
    for j in 1:nev
        ψ = wavefuncs[:, j]
        normfac = sqrt(sum(abs2, ψ) * dx)
        wavefuncs[:, j] .= ψ ./ normfac
        if imag(sum(conj.(wavefuncs[:, j]) .* wavefuncs[:, j])) ≈ 0
            wavefuncs[:, j] .= real(wavefuncs[:, j])
        end
    end

    return energies, wavefuncs, x
end

"""
Calculation of analytically known solutions of the translated single harmonic oscillator problem.
    Takes: 
    - x_range that is considered;
    - number of n states on the LHS well (canonical);
    - number of k states on the RHS well (tautomerical);
    - and the information whether its for A-T or G-C (bool).
    Returns:
    - left, right wavefunctions, left, right energies. 

"""
function get_wavefunctions_qho(x_range, n::Int, k::Int; is_at::Bool = true)
    a_can = is_at ? 0.01548757014342916 : 0.006457467585167605
    a_tau = is_at ? 0.0125386460155263 : 0.013406834311699567

    mass = 1836
    # 1/2 m ω^2=a-> mω=√2a/m * m=√2am
    dL = is_at ? -1.7737 : -2.42442
    dR = is_at ? 1.8963 : 1.787
    L = sqrt(a_can * mass)
    R = sqrt(a_tau * mass)
    
    left_wave_func = Matrix{Float64}(undef, length(x_range), n)
    right_wave_func = Matrix{Float64}(undef, length(x_range), k)
    
    for i in 1:n
        for (j, x) in enumerate(x_range)
            left_wave_func[j, i] = (2^(i-1) * factorial(i-1))^(-1/2) * (L^2/π)^(1/4) * 
                                   exp(-L^2 * (x - dL)^2/2) * eval_hermite(i-1, sqrt(L^2) * (x - dL))
        end
    end
    
    for i in 1:k
        for (j, x) in enumerate(x_range)
            right_wave_func[j, i] = (2^(i-1) * factorial(i-1))^(-1/2) * (R^2/π)^(1/4) * 
                                    exp(-R^2 * (x - dR)^2/2) * eval_hermite(i-1, sqrt(R^2) * (x - dR))
        end
    end
    
    energies_left = Vector{Float64}(undef, n)
    energies_right = Vector{Float64}(undef, k)
    ene_tau = is_at ? 0.021 : 0.0158
    
    for i in 1:n
        energies_left[i] = L/mass * ((i-1) + 0.5)
    end
    for i in 1:k 
        energies_right[i] = R/mass * ((i-1) + 0.5) + ene_tau
    end
    
    return left_wave_func, right_wave_func, energies_left, energies_right
end

"""
Simplified function from the above but used for changing R and L to determine the energies.
    Returns energies left and right for a given L, R param.
"""
function get_energy_levels(L::Float64, R::Float64, n_max::Int, k_max::Int, is_at::Bool = true)
    mass = 1836
    energies_left = Vector{Float64}(undef, n_max)
    energies_right = Vector{Float64}(undef, k_max)
    ene_tau = is_at ? 0.021 : 0.0158
    for i in 1:n_max
        energies_left[i] = L/mass * ((i-1) + 0.5)
    end
    for i in 1:k_max
        energies_right[i] = R/mass * ((i-1) + 0.5) + ene_tau
    end
    return energies_left, energies_right
end
##
function plot_harmonic_solutions(scale=0.0006, is_at::Bool = true)
    energies, wavefuncs, x = is_at ? solve_schrodinger_sum_harmonic(14, 1000,(-3.5, 3.0), true) : solve_schrodinger_sum_harmonic(15, 1000,(-4.3, 2.9), false)
    
    V = get_potential_model(x, is_at = is_at)

    fig = Figure(resolution = (800, 600))

    ax = Axis(fig[1, 1], xlabel = L"$x$ (a.u.)", ylabel = L"$E$ \text{ (a.u.)}",
        title = is_at ? L"\text{A-T: model harmoniczny}" : L"\text{G-C: model harmoniczny}",
        limits = is_at ? ((-3.3, 2.9), (-0.005, 0.045)) : ((-4.3, 2.5), (-0.002, 0.032)),
        ylabelsize = 35, xlabelsize = 35, titlesize = 35,
        xticklabelsize = 25, yticklabelsize = 25)
    ax_wf = Axis(fig[1, 1], ylabel = L"\psi(x) \text{ (arb. u.)}", ylabelsize=35,
        yticklabelsize = 25, yaxisposition = :right, limits = is_at ? ((-3.2, 3.0), (-1.2, 1.2)) : ((-4., 2.9), (-1.2, 1.2)))
    hidespines!(ax_wf)
    hidexdecorations!(ax_wf)
    cm = cgrad(:tab20c, 13)
    CairoMakie.lines!(ax, x, V, linewidth = 3.5, color = cm[1])
    
    cm = cgrad(:YlGnBu, 18)
    for (i, E) in enumerate(energies)
        ψ = wavefuncs[:, i]
        CairoMakie.lines!(ax, x, E .+ ψ .* scale, linewidth=2., color = cm[i+1], label="ψ[$(i-1)]")
        CairoMakie.lines!(ax, x, [E], linestyle=:dash, linewidth=1.5, label=false, color = cm[i+1], alpha = 0.5) 
    end
    
    # cm = cgrad(:PuRd, 5)
    # for (i, E) in enumerate(ene_right)
    #     ψ = right_wf[:, i]
    #     CairoMakie.lines!(ax, x, E .+ ψ .* scale, linewidth=2., color = cm[i+1], label="ψ[$(i-1)]")
    #     CairoMakie.lines!(ax, x, [E], linestyle=:dash, linewidth=1.5, label=false, color = cm[i+1], alpha = 0.5) 
    # end

    # display(fig)
    # filename = is_at ? "model_AT" : "model_GC"
    # save("graphics/model/$filename.pdf", fig)
end

function plot_solutions_with_density(scale=0.0007, is_at::Bool = true)
    x = is_at ? LinRange(-3.0, 3.0, 300) : LinRange(-4.0, 2.9, 300)
    left_wf, right_wf, ene_left, ene_right = is_at ? get_wavefunctions_qho(x, 12, 5) : get_wavefunctions_qho(x, 14, 4, is_at = false)
    V = get_potential_model(x, is_at = is_at)

    fig = Figure(resolution = (800, 600))

    ax = Axis(fig[1, 1], xlabel = L"$x$ (a.u.)", ylabel = L"\text{Energy (a.u.)}",
        title = is_at ? L"\text{A-T: harmonic model}" : L"\text{G-C: harmonic model}",
        limits = is_at ? ((-3.2, 3.0), (-0.005, 0.045)) : ((-4., 2.7), (-0.002, 0.032)),
        ylabelsize = 30, xlabelsize = 30, titlesize = 30,
        xticklabelsize = 20, yticklabelsize = 20)
    ax_wf = Axis(fig[1, 1], ylabel = L"|ψ(x)|^2 \text{ (arb. u.)}", ylabelsize=30,
        yticklabelsize = 20, yaxisposition = :right, limits = is_at ? ((-3.2, 3.0), (-1.2, 1.2)) : ((-4., 2.9), (-1.2, 1.2)))
    hidespines!(ax_wf)
    hidexdecorations!(ax_wf)
    cm = cgrad(:tab20c, 13)
    CairoMakie.lines!(ax, x, V, linewidth = 3.5, color = cm[1])
    
    cm = cgrad(:YlGnBu, 13)
    for (i, E) in enumerate(ene_left)
        ψ = left_wf[:, i]
        ρ = conj(ψ) .*ψ
        CairoMakie.lines!(ax, x, E .+ ρ .* scale, linewidth=2., color = cm[i+1], label="ψ[$(i-1)]")
        CairoMakie.lines!(ax, x, [E], linestyle=:dash, linewidth=1.5, label=false, color = cm[i+1], alpha = 0.5) 
    end
    
    cm = cgrad(:PuRd, 5)
    for (i, E) in enumerate(ene_right)
        ψ = right_wf[:, i]
        ρ = conj(ψ) .*ψ
        CairoMakie.lines!(ax, x, E .+ ρ .* scale, linewidth=2., color = cm[i+1], label="ψ[$(i-1)]")
        CairoMakie.lines!(ax, x, [E], linestyle=:dash, linewidth=1.5, label=false, color = cm[i+1], alpha = 0.5) 
    end

    display(fig)
    # filename = is_at ? "model_AT_den" : "model_GC_den"
    # save("graphics/model/$filename.pdf", fig)
end

# ene_at, wf_at, x_at = solve_schrodinger_sum_harmonic(13, 1000, (-3.5, 3.), true)
# plot_harmonic_solutions()

# g-C
# ene_gc, wf_gc, x_gc = solve_schrodinger_sum_harmonic(14, 1000, (-4., 2.9), false)
# plot_harmonic_solutions(0.0004, false)

##
# plot_harmonic_solutions()
# plot_solutions_with_density(0.0007, false)
