include("schrodinger_eqn.jl")
using DelimitedFiles
using LinearAlgebra
using FFTW
using CairoMakie
"""
Computing the WDF out of wavefunctions
"""

function compute_discrete_wdf(ψ, x::AbstractVector{<:Real})
    N = length(x)
    @assert N == length(ψ) "psi must be equally length"
    dx = x[2] - x[1]
    @assert all(abs.(diff(x) .- dx) .< 1e-12) "x uniformed"

    pad = N 
    ψ_pad = vcat(zeros(ComplexF64, pad), ComplexF64.(ψ), zeros(ComplexF64, pad))
    offset = pad # we expand the wavefunction so that we can access x-y/2 arguments

    if N % 2 != 0
        println("Please use even N")
    end
    N2 = div(N, 2)
    n_range = -N2:(N2-1)
    M = length(n_range)

    ## wdf
    W = zeros(Float64, N, N)
    c = Array{ComplexF64}(undef, M)

    for j in 1:N
        b = offset + j 
        for (idx, n) in enumerate(n_range)
            left = ψ_pad[b-n]
            right = ψ_pad[b+n]
            c[idx] = conj(left) * right
        end
        ck = fft(c)
        constan = dx/π
        W[j, :] = constan .* real(fftshift(ck))
    end
    kvals = collect(-N2:(N2-1))
    p = (π /(N * dx)).* kvals
    return W, x, p
end

function fftshift(v::AbstractVector)
    N = length(v)
    p = fld(N, 2)
    return vcat(v[begin+p+1:end], v[begin:begin+p])
end


function get_wdf_map(idx_of_can_state::Int64 = 6, idx_of_tau_state::Int64 = 7, is_at::Bool = true)
    ψ0 = is_at ? wf_at[:, idx_of_can_state] : wf_gc[:, idx_of_can_state]
    ψ1 = is_at ? wf_at[:, idx_of_tau_state] : wf_gc[:, idx_of_tau_state]  
    x = is_at ? x_at : x_gc
    W, xout, p = compute_discrete_wdf(ψ0, x)
    W7, _, _ = compute_discrete_wdf(ψ1, x)
    W = W + W7
    fig = Figure(resolution=(700,500))
    title = is_at ? "wdf_at_godbeer" : "wdf_gc_slocombe" 
    ax = Axis(fig[1,1], xlabel=L"$x$ (a.u.)", ylabel=L"$p$ (a.u.)", title=is_at ? L"\text{WDF 5$\rightarrow$0 (A-T)}" : L"\text{WDF 5$\rightarrow$0 (G-C)}", 
        limits = is_at ? ((-3.2, 3.), (-12, 12)) : ((-4., 2.7), (-12, 12)), ylabelsize = 30, xlabelsize = 30, titlesize = 30,
            xticklabelsize = 20, yticklabelsize = 20)
    hm = CairoMakie.heatmap!(ax, xout, p, W, colormap = :seismic, colorrange = (-0.3, 0.3))
    CairoMakie.Colorbar(fig[1,2], hm, label = L"\varrho(x,\; p; \; t)", labelsize = 30, ticklabelsize = 20)
    save("graphics/$title.pdf", fig)
end

get_wdf_map(6, 7)
get_wdf_map(6, 7, false)