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
    color1 = is_at ? :darkorange2 : :darkorchid4
    color2 = is_at ? :steelblue3 : :seagreen
    
    
    Z1 = [get_quantum_model_Z(T; is_at=is_at) for T in T_range]
    Z2 = [get_quantum_real_Z(T; is_at=is_at) for T in T_range]
    
    lines!(ax1, T_range, Z1, linewidth=5, label = L"\text{model harmoniczny}", color = color1)
    lines!(ax1, T_range, Z2, linewidth=5, label = L"\text{dokładne rozwiązanie}", color = color2)
    if is_at==false
        ax1.yticks = [0.0, 0.2, 0.4, 0.6]
        ax1.title = L"\text{G-C}"
    end

    Legend(fig[2,1], ax1, framevisible = false, orientation = :horizontal, labelsize = 35)
    save("graphics/part_func$is_at.pdf", fig)
    return fig
end

fig = get_temp_dependence_Z(is_at = false)

function get_USF(is_HO::Bool;is_at::Bool = true)
    T_range = LinRange(100, 400, 100)
    x = LinRange(-4.0, 3.0, 100)
    if is_HO
        _, _, ene_left, ene_right = is_at ? get_wavefunctions_qho(x, 12, 5; is_at=true) : get_wavefunctions_qho(x, 14, 4; is_at=false)
    else
        energies, _, _ = is_at ? solve_schrodinger(13, 1000, (-3.5, 3.0), true) : solve_schrodinger(12, 1000, (-4., 2.9), false)
    end
    # Calculate thermodynamic quantities for each temperature
    F = Float64[]  # Free energy
    U = Float64[]  # Internal energy  
    S = Float64[]  # Entropy
    
    for T in T_range
        β = 1.0 / (kB * T) 
        
        # Partition function
        Z = get_quantum_model_Z(T; is_at=is_at)
        
        # Free energy: F = -kB T ln(Z)
        free_energy = -(kB * T) * log(Z)
        push!(F, free_energy)
        
        # Internal energy: U = <E> = Σᵢ Eᵢ exp(-βEᵢ) / Z
        weighted_energy = 0.0
        if is_HO
            for ene in ene_left
                weighted_energy += ene * exp(-β * ene)
            end
            for ene in ene_right
                weighted_energy += ene * exp(-β * ene)
            end
        else
            for ene in energies
                weighted_energy += ene * exp(-β * ene)
            end  
        end
        internal_energy = weighted_energy / Z
        push!(U, internal_energy)
        
        entropy = kB * log(Z) + internal_energy / T
        push!(S, entropy)
    end
    
    return T_range, F, U, S
end

function get_other_ther_funcs(;is_at::Bool = true)
    fig = Figure(size = (1000, 800))
    ax1 = Axis(fig[1,1], 
    # xlabel = L"$T$ [K]",xlabelsize = 35, 
    ylabel = L"$Z$ [a.u.]", xticklabelsize = 0,
    # title = L"\text{Funkcja rozdziału}", titlesize = 35, 
    ylabelsize = 35, 
    yticklabelsize = 30)
    
    ax2 = Axis(fig[1,2], 
    # xlabel = L"$T$ [K]", xlabelsize = 35, 
    ylabel = L"$F$ [$10^{-3}$ a.u.]", xticklabelsize = 0,
    # title = L"\text{Energia swobodna}", titlesize = 35, 
    ylabelsize = 35, yticklabelsize = 30)

    ax3 = Axis(fig[2,1], xlabel = L"$T$ [K]", ylabel = L"$U$ [$10^{-3}$ a.u.]",
    # title = L"\text{Energia wewnętrzna}", titlesize = 35, 
    xlabelsize = 35, ylabelsize = 35, 
    xticklabelsize = 30, yticklabelsize = 30)

    ax4 = Axis(fig[2,2], xlabel = L"$T$ [K]", ylabel = L"$S$ [$10^{-6}$ a.u.]",
    # title = L"\text{Entropia}", titlesize = 35, 
    xlabelsize = 35, ylabelsize = 35, 
    xticklabelsize = 30, yticklabelsize = 30)

    T_range, F, U, S = get_USF(true; is_at=is_at)
    T_range2, F2, U2, S2 = get_USF(false; is_at=is_at)
    color1 = is_at ? :darkorange2 : :darkorchid4
    color2 = is_at ? :steelblue3 : :seagreen
    
    lines!(ax2, T_range, F*1e03, linewidth=10, color = color1)
    lines!(ax3, T_range, U*1e03, linewidth=5, color = color1)
    lines!(ax4, T_range, S*1e06, linewidth=5, color = color1)
    
    lines!(ax2, T_range2, F2*1e03, linewidth=5, color = color2)
    lines!(ax3, T_range2, U2*1e03, linewidth=5, color = color2)
    lines!(ax4, T_range2, S2*1e06, linewidth=5, color = color2)
    
    T_range = LinRange(100, 400, 100)
    
    Z1 = [get_quantum_model_Z(T; is_at=is_at) for T in T_range]
    Z2 = [get_quantum_real_Z(T; is_at=is_at) for T in T_range]
    
    lines!(ax1, T_range, Z1, linewidth=5, color = color1, label = L"\text{model harmoniczny}")
    lines!(ax1, T_range, Z2, linewidth=5, color = color2, label = L"\text{dokładne rozwiązanie}")
    
    # hidespines!(ax1)
    # hidexdecorations!(ax1)
    # hidespines!(ax2)
    # hidexdecorations!(ax2)
    if is_at==false
        ax1.yticks = [0.0, 0.2, 0.4, 0.6]
    end
    Legend(fig[3, 1:2], ax1, framevisible = false, orientation = :horizontal, labelsize = 35)
    save("graphics/term_funcs_at$is_at.pdf", fig)
    return fig
end
# get_other_ther_funcs(;is_at = false)