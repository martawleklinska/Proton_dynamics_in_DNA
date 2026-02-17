using CairoMakie 
using LinearAlgebra 
using DelimitedFiles

include("../true_calc_at/schrodinger_eqn.jl")
include("schrodinger_eqn.jl")
kB = 3.167e-6  # Hartree/K

function get_quantum_model_Z(T::Float64; is_at::Bool = true)
    x = LinRange(-4.0, 3.0, 100)
    _, _, ene_left, ene_right = is_at ? get_wavefunctions_qho(x, 12, 5; is_at=true) : get_wavefunctions_qho(x, 14, 4; is_at=false)
    
    β = 1.0 / (kB * T)  
    Z = 0.0

    # partition func = \frac{kB T}{\hbar \omega0} exp(-β V_0)
    ωb = is_at ? sqrt(2*0.010377758242695038/1836) : sqrt(2*0.006425438605566449/1836)
    V0 = is_at ? 0.029636925 : 0.0246348381
    Z += 1/(β * ωb) * exp(-β * V0) 
    for ene in ene_left
        Z += exp(-β * ene)
    end
    
    for ene in ene_right
        Z += exp(-β * ene)
    end
    
    return Z
end

function get_quantum_real_Z(T::Float64; is_at::Bool = true)
    energies, _, _ = is_at ? solve_schrodinger(13, 1000, (-3.5, 3.0), true) : solve_schrodinger(12, 1000, (-4., 2.9), false)
    
    β = 1.0 / (kB * T)
    Z = 0.0
    
    for ene in energies
        Z += exp(-β * ene)
    end
    
    return Z
end

function get_temp_dependence_Z(;is_at::Bool = true)
    fig = Figure(size = (1000, 500))
    ax1 = Axis(fig[1,1], xlabel = L"$T$ [K]", ylabel = L"$Z$", title = L"\text{A-T}", 
    xlabelsize = 35, ylabelsize = 35, titlesize = 35, xticklabelsize = 30, yticklabelsize = 30)
    # ax2 = Axis(fig[1,2], xlabel = L"$T$ [K]", ylabel = L"$Z$", title = L"\text{Dokładne rozwiązanie}", 
    # xlabelsize = 25, ylabelsize = 25, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    
    T_range = LinRange(100, 400, 100)
    
    Z1 = [get_quantum_model_Z(T; is_at=is_at) for T in T_range]
    Z2 = [get_quantum_real_Z(T; is_at=is_at) for T in T_range]
    
    lines!(ax1, T_range, Z1, linewidth=3, label = L"\text{model harmoniczny}")
    lines!(ax1, T_range, Z2, linewidth=3, label = L"\text{dokładne rozwiązanie}")
    
    ax1.title = is_at ? L"\text{A-T}" : L"\text{G-C}"
    Legend(fig[2,1], ax1, framevisible = false, orientation = :horizontal, labelsize = 35)
    return fig
end

fig = get_temp_dependence_Z(is_at = false)
