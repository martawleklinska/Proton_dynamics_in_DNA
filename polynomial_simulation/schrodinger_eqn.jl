using LinearAlgebra, SparseArrays
# using Arpack
using CairoMakie

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
function get_potential(x; is_at=true)
    return is_at ? model_at.(x) : model_gc.(x)
end

""" 
TISE with FD returns: energies, wavefunctions, x
- takes:
    nstates: n/o of lowest eigenstates
    n: grid plot_at_instance
    xlims: (xmin, xmax)

"""
function solve_schrodinger(nstates::Int64=10, n::Int64=1000, xlims::Tuple{Float64, Float64}=(-3.5, 3.0),
                            is_at::Bool=true)
    xmin, xmax = xlims
    x = range(xmin, stop = xmax, length = n)
    dx = step(x)
    nx = length(x)
    # potential matrix
    main = fill(1.0 / dx ^ 2 * 1/1836, nx)
    off = fill(-0.5 / dx ^ 2 * 1/1836, nx - 1)

    T = spdiagm(-1 => off, 0 => main, 1 => off)

    V_diag = get_potential(collect(x), is_at = is_at)
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


#

function plot_solutions(ene, wavefuncs, x; scale=0.01, is_at::Bool = true)
    V = get_potential(x, is_at = is_at)
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"$x$ (a.u.)", ylabel = L"\text{Energy (a.u.)}",
        title = is_at ? L"\text{A-T: harmonic model}" : L"\text{G-C: harmonic model}",
        limits = is_at ? ((-3.2, 3.0), (-0.005, 0.045)) : ((-4., 2.7), (-0.002, 0.035)),
        ylabelsize = 30, xlabelsize = 30, titlesize = 30,
        xticklabelsize = 20, yticklabelsize = 20)
    ax_wf = Axis(fig[1, 1], ylabel = L"\psi(x) \text{ (arb. u.)}", ylabelsize=30,
        yticklabelsize = 20, yaxisposition = :right, limits = is_at ? ((-3.2, 3.0), (-1.2, 1.2)) : ((-4., 2.9), (-1.2, 1.2)))
    hidespines!(ax_wf)
    hidexdecorations!(ax_wf)
    cm = cgrad(:tab20c, 13)
    CairoMakie.lines!(ax, x, V, linewidth = 2.5, color = cm[1])

    for (i, E) in enumerate(ene)
        ψ = wavefuncs[:, i]
        CairoMakie.lines!(ax, x, E .+ ψ .* scale, linewidth=2., color = cm[i+1], label="ψ[$(i-1)]")
        CairoMakie.lines!(ax, x, [E], linestyle=:dash, linewidth=1.5, label=false, color = cm[i+1]) 
    end

    # display(fig)
    filename = is_at ? "model_AT" : "model_GC"
    save("graphics/model/$filename.pdf", fig)
end

ene_at, wf_at, x_at = solve_schrodinger(13, 1000, (-3.5, 3.), true)
plot_solutions(ene_at, wf_at, x_at; scale=0.001, is_at = true)

## g-C
ene_gc, wf_gc, x_gc = solve_schrodinger(14, 1000, (-4., 2.9), false)
plot_solutions(ene_gc, wf_gc, x_gc; scale = 0.001, is_at = false)
