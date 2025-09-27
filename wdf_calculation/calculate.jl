using DelimitedFiles

function eval_genlaguerre(n::Int, x::Float64)
    if n == 0
        return one(x)
    elseif n == 1
        return 1 - x
    else
        l0 = one(x)       
        l1 = 1 - x        
        for k in 1:(n-1)
            l2 = ((2k + 1 - x) * l1 - k * l0) / (k + 1)
            l0, l1 = l1, l2
        end
        return l1
    end
end


function WDF_func(n, x, p, R; d=0.0)
    x_shift = x .- d
    zeta = 2 * x_shift.^2 * R + 2 * p.^2 / R
    exponen = exp.(-R * x_shift.^2 .- (1 / R) * p.^2)
    Laguerre = eval_genlaguerre(n, zeta)
    return (-1)^n / π * exponen .* Laguerre
end

function compute_WDF_map(x_vals, p_vals, dL, dR, L, R; n=8, k=0)
    Z = zeros(length(x_vals), length(p_vals))

    for (i, x) in enumerate(x_vals)
        for (j, p) in enumerate(p_vals)
            Z[i, j] = WDF_func(n, x, p, L, d=dL) + WDF_func(k, x, p, R, d=dR)
        end
    end
    return Z
end

# getting WDF from the time-varying potential
dL, dR = -1.968, 1.968

function compute_WDF_for_time_varying_potential(x_vals::Vector{Float64}, p_vals::Vector{Float64}, 
    L_series::Vector{Float64}, R_series::Vector{Float64},
    n::Int64, k::Int64; title = "gc")

    nx, np, nframes = length(x_vals), length(p_vals), length(R_series)
    WDF_series = Array{Float64}(undef, nx, np, nframes)

    for (i, R) in enumerate(R_series)
        WDF_series[:, :, i] .= compute_WDF_map(x_vals, p_vals, dL, dR, L_series[i], R, n=n, k=k)
    end

    # writedlm("data/WDF_series_$title.txt", reshape(WDF_series, nx*np, nframes))
    return WDF_series
end

