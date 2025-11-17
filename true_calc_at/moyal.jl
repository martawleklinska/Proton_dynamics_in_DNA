using FFTW, LinearAlgebra, Printf

include("wdf_definition.jl")

const ħ = 1.0
const im = complex(0,1)

using FFTW
using LinearAlgebra

Nx = 1000
dx = 0.05
x = (-Nx/2:Nx/2-1) .* dx

# momentum grid (FFT dual)
dp = π/(Nx*dx)
p = (-Nx/2:Nx/2-1) .* dp

ψ0 = wf_at[:, 6]          # for example, canonical state
W0, xgrid, pgrid = compute_discrete_wdf(ψ0, x)
W = copy(W0)                  # this will be evolved in time

# read a(t) series once (time-indexed series of coefficients)
using DelimitedFiles
const a_series = readdlm("data/L_series_at.txt") |> vec

function V_xt(x, t, a_val)
    # a_val is a scalar value of the time-dependent coefficient a(t)
    return 0.5 * a_val * x.^2
end

function get_a(idx::Int)
    # index-based accessor for the preloaded series
    idx = clamp(idx, 1, length(a_series))
    return a_series[idx]
end

# V = V_xt.(x, t, a)
# Vhat = fftshift(fft(ifftshift(V))) * dx   # Fourier convention for convolution

function stream_step!(W, p, dx, dt, m)
    # shift in x via FFT
    Wx = fft(W, 1)
    phase = exp.(-im .* (p .* dt ./ m) .* (2π .* (0:size(W,1)-1) / (size(W,1)*dx)))
    Wx .*= phase
    W .= real(ifft(Wx, 1))
end

function potential_step!(W, Vhat, dp, dt)
    # convolution in p via FFT in the second dimension
    Wp = fft(W, 2)
    Wp .*= exp.(-im .* Vhat .* dt)   # effective phase (simplified form)
    W .= real(ifft(Wp, 2))
end

m = 1836.0          # proton mass in a.u., change if you want KIE
dt = 0.0005
Nt = 1000

W = copy(W0)

for n in 1:Nt
    t = (n-1)*dt

    # --- build potential & FFT it ---
    a_val = get_a(n)
    V = V_xt.(x, t, a_val)
    Vhat = fftshift(fft(ifftshift(V))) * dx

    # --- Strang Splitting ---
    potential_step!(W, Vhat, dp, dt/2)
    stream_step!(W, p, dx, dt, m)
    potential_step!(W, Vhat, dp, dt/2)

    # save / visualize occasionally
    if n % 200 == 0
        println("t = $t")
    end
end

function run_wigner_sim(ψ0, x, p; dt=0.0005, Nt=5000, m=1836., a_series_in = nothing)
    W, _, _ = compute_discrete_wdf(ψ0, x)
    W = copy(W)

    dx = x[2] - x[1]
    dp = p[2] - p[1]

    # choose series: provided or fallback to global a_series
    a_s = isnothing(a_series_in) ? a_series : a_series_in

    for n in 1:Nt
        t = (n-1)*dt
        a_val = clamp(n, 1, length(a_s)) <= length(a_s) ? a_s[clamp(n,1,length(a_s))] : a_s[end]
        V = V_xt.(x, t, a_val)
        Vhat = fftshift(fft(ifftshift(V))) * dx

        potential_step!(W, Vhat, dp, dt/2)
        stream_step!(W, p, dx, dt, m)
        potential_step!(W, Vhat, dp, dt/2)
    end

    return W
end

ψ0 = wf_at[:, 6]
W_final = run_wigner_sim(ψ0, x, p; dt=0.0005, Nt=3000)
##
function get_wdf_ma(idx_of_can_state::Int64 = 6, idx_of_tau_state::Int64 = 7, is_at::Bool = true)
    x = is_at ? x_at : x_gc
    # W, xout, p = compute_discrete_wdf(ψ0, x)
    # W7, _, _ = compute_discrete_wdf(ψ1, x)
    W = W_final
    fig = Figure(resolution=(700,500))
    # title = is_at ? "wdf_at_godbeer" : "wdf_gc_slocombe" 
    ax = Axis(fig[1,1], xlabel=L"$x$ (a.u.)", ylabel=L"$p$ (a.u.)",# title=is_at ? L"\text{WDF 5$\rightarrow$0 (A-T)}" : L"\text{WDF 5$\rightarrow$0 (G-C)}", 
        limits = is_at ? ((-3.2, 3.), (-12, 12)) : ((-4., 2.7), (-12, 12)), ylabelsize = 30, xlabelsize = 30, titlesize = 30,
            xticklabelsize = 20, yticklabelsize = 20)
    hm = CairoMakie.heatmap!(ax, x, p, W, colormap = :seismic, colorrange = (-0.3, 0.3))
    CairoMakie.Colorbar(fig[1,2], hm, label = L"\varrho(x,\; p; \; t)", labelsize = 30, ticklabelsize = 20)
    # save("graphics/true_sim/$title.pdf", fig)
    return fig
end
get_wdf_ma()