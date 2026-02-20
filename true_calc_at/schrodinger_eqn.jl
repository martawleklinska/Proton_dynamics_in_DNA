using LinearAlgebra, SparseArrays
using Arpack
using CairoMakie

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


function get_potential(x; is_at=true)
    return is_at ? fourth_order.(x) : double_morse.(x)
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


##

function plot_solutions(ene, wavefuncs, x; scale=0.01, is_at::Bool = true)
    V = get_potential(x, is_at = is_at)
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"$x$ (a.u.)", ylabel = L"$E$ \text{ (a.u.)}",
        title = is_at ? L"\text{A-T}" : L"\text{G-C}",
        limits = is_at ? ((-3.2, 3.0), (-0.005, 0.045)) : ((-4., 2.7), (-0.002, 0.035)),
        ylabelsize = 35, xlabelsize = 35, titlesize = 35,
        xticklabelsize = 25, yticklabelsize = 25)
    ax_wf = Axis(fig[1, 1], ylabel = L"\psi(x) \text{ (arb. u.)}", ylabelsize=35,
        yticklabelsize = 25, yaxisposition = :right, limits = is_at ? ((-3.2, 3.0), (-1.2, 1.2)) : ((-4., 2.9), (-1.2, 1.2)))
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
    filename = is_at ? "fourth_order_wave_funcs_AT" : "morse_wave_funcs_GC"
    save("graphics/true_sim/$filename.pdf", fig)
end

# ene_at, wf_at, x_at = solve_schrodinger(12, 1000, (-3.5, 3.), true)
# plot_solutions(ene_at, wf_at, x_at; scale=0.001, is_at = true)

# ## g-C
# ene_gc, wf_gc, x_gc = solve_schrodinger(12, 1000, (-4., 2.9), false)
# plot_solutions(ene_gc, wf_gc, x_gc; scale = 0.001, is_at = false)
