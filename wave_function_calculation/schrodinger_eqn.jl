using LinearAlgebra, SparseArrays
# using Arpack
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


function get_potential(x)
    return fourth_order.(x)
end

""" 
TISE with FD returns: energies, wavefunctions, x
- takes:
    nstates: n/o of lowest eigenstates
    n: grid plot_at_instance
    xlims: (xmin, xmax)

"""
function solve_schrodinger_fourth(nstates::Int64=10, n::Int64=1000, xlims::Tuple{Float64, Float64}=(-3.5, 3.0))
    xmin, xmax = xlims
    x = range(xmin, stop = xmax, length = n)
    dx = step(x)
    nx = length(x)
    # potential matrix
    main = fill(1.0 / dx ^ 2 * 1/1836, nx)
    off = fill(-0.5 / dx ^ 2 * 1/1836, nx - 1)

    T = spdiagm(-1 => off, 0 => main, 1 => off)

    V_diag = get_potential(collect(x))
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

ene, wf, x = solve_schrodinger_fourth(4, 1000, (-3.5, 3.0))
println("Lowest energies a.u.:")
println(ene)

##

function plot_solutions(ene, wavefuncs, x; scale=0.01)
    V = get_potential(x)
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"$x$ (a.u.)", ylabel = L"\text{Energy (a.u.)}",
        title = L"\text{A-T}",
        limits = ((-3.2, 3.0), (-0.005, 0.045)),
        ylabelsize = 30, xlabelsize = 30, titlesize = 30,
        xticklabelsize = 20, yticklabelsize = 20)
    ax_wf = Axis(fig[1, 1], ylabel = L"\psi(x) \text{ (arb. u.)}", ylabelsize=30,
        yticklabelsize = 20, yaxisposition = :right, limits = ((-3.2, 3.0), (-1.2, 1.2)))
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
    save("graphics/fourth_order_wave_funcs_AT.pdf", fig)
end

ene, wf, x = solve_schrodinger_fourth(12, 1000, (-3.5, 3.0))
plot_solutions(ene, wf, x; scale=0.001)
