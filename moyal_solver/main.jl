using CairoMakie, DelimitedFiles, Printf
using StatsBase
using GLMakie

const dt = .1  # time step 

function create_wigner_animation(run_sim::Bool = true)
    output_paths = [
        # "build/output/",  
        # "masters/moyal_solver/build/output/",               
        # "moyal_solver/build/output/", 
        "moyal_solver/build/output/"             
    ]
    
    output_dir = nothing
    for path in output_paths
        if isdir(path)
            output_dir = path
            println("Found output directory at: $path")
            break
        end
    end
    
    wigner_files = filter(f -> startswith(f, "wigner_") && endswith(f, ".dat"), 
                         readdir(output_dir))
    sort!(wigner_files)
    
    
    first_file = joinpath(output_dir, wigner_files[1])
    data = readdlm(first_file)
    x_coords, p_coords = data[:, 1], data[:, 2]
    
    x_unique = sort(unique(x_coords))
    p_unique = sort(unique(p_coords))
    nx, np = length(x_unique), length(p_unique)
    
    w_min, w_max = Inf, -Inf
    sample_files = wigner_files[1:max(1, length(wigner_files)÷10):end]  
    
    for filename in sample_files
        data = readdlm(joinpath(output_dir, filename))
        w_vals = data[:, 3]
        w_min = min(w_min, minimum(w_vals))
        w_max = max(w_max, maximum(w_vals))
    end
    println("Wigner range: $(round(w_min, digits=4)) to $(round(w_max, digits=4))")
    
    animation_files = wigner_files[1:1:end]
    n_frames = length(animation_files)
        
    fig = Figure(size = (900, 700), fontsize = 16)
    
    time_obs = Observable("t = 0.0")
    wigner_obs = Observable(zeros(nx, np))
    
    if run_sim
        ax = Axis(fig[1, 1], 
                xlabel = L"\text{Położenie}\; x", 
                ylabel = L"\text{Pęd} \;p",
                title = time_obs,
                titlesize = 25,
                limits = ((-5., 5.), (-12., 12.)),
                xlabelsize = 25,
                ylabelsize = 25)
        
        hm = heatmap!(ax, x_unique, p_unique, wigner_obs,
                        colormap = :RdBu,
                        colorrange = (-0.2, 0.2))

        Colorbar(fig[1, 2], hm, label = L"\varrho(x,p; t)", labelsize = 25)
        
        gif_filename = "moyal_solver/graphics/GC/wdf_evolution.gif"
        record(fig, gif_filename, 1:n_frames; framerate = 8) do frame_idx
            filename = animation_files[frame_idx]
            
            try
                step_str = match(r"wigner_(\d+)\.dat", filename).captures[1]
                step = parse(Int, step_str)
                time_val = step * dt
                
                data = readdlm(joinpath(output_dir, filename))
                wigner_real = data[:, 3]
                
                W = reshape(wigner_real, np, nx)'
                W_vis = W
                
                time_obs[] = @sprintf("Funkcja Wignera (t = %.3f)", time_val)
                wigner_obs[] = W_vis
                
                if frame_idx % max(1, n_frames÷20) == 0
                    progress = round(100 * frame_idx / n_frames, digits=1)
                    println("Progress: $progress% (frame $frame_idx/$n_frames)")
                end
            catch e
                println("Warning: Error processing frame $frame_idx ($filename): $e")
            end
        end
    end
end

create_wigner_animation()
##
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
        return a_can .* x .^ 2 + b_can .* x + c_can .-0.0018
    elseif x < 1.15
        return a_bar .* x .^ 2 + b_bar .* x + c_bar .-0.0018
    else 
        return a_tau .* x .^ 2 + b_tau .* x + c_tau .-0.0018
    end
end
##
function plot_wigner_snapshots(; is_harmonic::Bool=false, is_gc::Bool=false, is_at::Bool=true)
    if is_at
        output_dir = "moyal_solver/build_at_short_sim/output/"
    elseif is_gc
        output_dir = "moyal_solver/build_gc_model_tautomeric/output/"
    end
    isdir(output_dir) || error("Output directory not found: $output_dir")

    wigner_files = filter(f -> startswith(f, "wigner_") && endswith(f, ".dat"), 
                         readdir(output_dir))
    sort!(wigner_files)

    # ======== wczytanie siatki z pierwszego pliku =================================
    data      = readdlm(joinpath(output_dir, wigner_files[1]))
    x_unique  = sort(unique(data[:, 1]))
    p_unique  = sort(unique(data[:, 2]))
    nx, np    = length(x_unique), length(p_unique)

    # ======== potencjał i poziomice Hamiltonianu ==============================
    if is_harmonic
        m      = 1.0
        Vx     = @. 0.5 * x_unique^2
        levels = range(0, 20, length=11)
        label  = "HO"
    elseif is_gc
        m  = 1836.0
        V1, V2      = 0.1617, 0.082
        a1_gc, a2_gc = 0.305, 0.755
        r1, r2      = -2.7, 2.1
        Vx     = @. V1*(exp(-2a1_gc*(x_unique-r1)) - 2exp(-a1_gc*(x_unique-r1))) +
                    V2*(exp(-2a2_gc*(r2-x_unique)) - 2exp(-a2_gc*(r2-x_unique))) + 0.166 + 0.00019
        levels = range(0, 0.1, length=15)
        label  = "GC"
        limits = ((-4., 2.7), (-11.5, 11.5))
        potential(x) = model_gc(x)
    elseif is_at
        m  = 1836.0
        levels = range(0, 0.1, length=15)
        label  = "AT"
        limits = ((-3.3, 2.85), (-12.5, 12.5))
        potential(x) = model_at(x)
    else
        error("Wybierz jeden potencjał: is_harmonic, is_gc lub is_at")
    end
    H = [(p^2)/(2m) + model_at(x) for p in p_unique, x in x_unique]

    # =========== zakres kolorów z próbki plików ==========================
    sample = wigner_files[1:max(1, length(wigner_files)÷10):end]
    w_min, w_max = mapreduce(
        f -> extrema(readdlm(joinpath(output_dir, f))[:, 3]),
        (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
        sample
    )
    println("Wigner range: $(round(w_min, digits=4)) to $(round(w_max, digits=4))")

    # =========== cztery pierwsze snapshoty ==========================
    wigner_files = filter(f -> startswith(f, "wigner_") && endswith(f, ".dat"), 
                         readdir(output_dir))
    sort!(wigner_files)
    # snap4        = wigner_files[1:40:min(1000, length(wigner_files))]
    snap4        = wigner_files[1:6:min(800, length(wigner_files))]

    # fig = Figure(size=(1000, 1200)) # for long sims 
    fig = Figure(size = (900, 800)) # short sims
    β = 0.01
    wigner_scale = ReversibleScale(
        w ->  asinh(w / β) / log(10),  
        w -> β * sinh(log(10) * w)     
    )

    for (idx, filename) in enumerate(snap4)
        row, col = divrem(idx - 1, 2) .+ 1          # (1,1) (1,2) (2,1) (2,2)
        # titles = [L"t = 0.0 \text{ a.u.}", L"t = 2\times 10^3 \text{ a.u.}", L"t = 4\times 10^3\text{ a.u.}", L"t = 6\times 10^3\text{ a.u.}", L"t = 8\times 10^3\text{ a.u.}", L"t = 10\times 10^3\text{ a.u.}"]
        titles = [L"t = 0.0 \text{ a.u.}", L"t = 200.0 \text{ a.u.}", L"t = 400.0\text{ a.u.}", L"t = 600.0\text{ a.u.}"]
        
        step     = parse(Int, match(r"wigner_(\d+)\.dat", filename).captures[1])
        t_val    = step * dt
        W        = reshape(readdlm(joinpath(output_dir, filename))[:, 3], np, nx)'

        ax = Axis(fig[row, col],
                  xlabel       = L"x \; (\text{a.u.})",
                  ylabel       = L"p \; (\text{a.u.})",
                  title        = titles[idx],
                  titlesize    = 30,
                  xlabelsize   = 30,
                  ylabelsize   = 30,
                  yticks = [-10., -5., 0., 5., 10.],
                  xticks = [-2., 0., 2.],
                  xticklabelsize = 25,
                  yticklabelsize = 25,
                  limits       = limits
                  )

        hm = heatmap!(ax, x_unique, p_unique, W;
        colormap   = :RdBu,
        colorscale = wigner_scale,
                      colorrange = (-0.32, 0.32))

        contour!(ax, x_unique, p_unique, H';
                 levels    = levels,      
                 linewidth = 1.2,
                 color     = :gray,
                 alpha     = 0.4)

        col == 2 && Colorbar(fig[row, col+1], hm;
                             label     = L"\varrho(x,p;\,t)\; (\text{a.u.})",
                             labelsize = 30, ticks = [-0.2, -0.1, 0.0, 0.1, 0.2], ticklabelsize = 25)
    end

    out_dir  = "moyal_solver/graphics/$label"
    mkpath(out_dir)
    # out_file = joinpath(out_dir, "wigner_snapshots_model.pdf")
    out_file = joinpath(out_dir, "wigner_snapshots_model_short_sim.pdf")
    save(out_file, fig; px_per_unit=2)
    println("Zapisano: $out_file")
end

plot_wigner_snapshots(; is_harmonic=false, is_gc=false, is_at=true)##
##
function create_nonclassicality_plot(;is_at::Bool = true, isnot_sole_nonclassicality::Bool = false)
    stats_file = is_at ? "moyal_solver/build_at_model/output/stats.dat" : "moyal_solver/build_gc_model/output/stats.dat"
    if !isfile(stats_file)
        println("Stats file not found, skipping nonclassicality plot")
        return nothing
    end
    nonclass_filename = is_at ? "moyal_solver/graphics/AT/nonclassicality.pdf" : "moyal_solver/graphics/GC/nonclassicality.pdf"
    
    data = readdlm(stats_file, skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]
    delta = size(data, 2) >= 8 ? data[:, 8] : zeros(length(t))  
    
    fig = Figure(size = (1000, 500))
    ## ====== axis delta ========
    ax = Axis(fig[1, 1],
              xlabel = L"t \;\; (10^3\text{ a.u.)} ",
              ylabel = L"\delta_{\mathrm{W}}(t)",
            #   title = L"\text{Ewolucja nieklasyczności stanu kwantowego}",
              titlesize = 30,
                      xlabelsize = 35,leftspinecolor = :purple, yticklabelcolor = :purple, ylabelcolor = :purple, ytickcolor = :purple,
                      ylabelsize = 35, xticklabelsize = 30, yticklabelsize = 30)
    if isnot_sole_nonclassicality
        ## ====== axis x ========
        ax1_color = :royalblue1
        ax1 = Axis(fig[1,1], xlabel = L"t\; \;(10^3\text{ a.u.})", ylabel = L"\langle x\rangle\; (\text{a.u.})", 
        xlabelsize = 35, ylabelsize = 35, xticklabelsize = 30, yticklabelsize = 30, yaxisposition = :right, 
        yticklabelcolor = ax1_color, ylabelcolor = ax1_color, ytickcolor = ax1_color, ylabelpadding  = 70)
        
        hidexdecorations!(ax1)

        ## ====== axis p ========
        ax2_color = :crimson
        ax2 = Axis(fig[1,1], ylabel = L"\langle p \rangle\; \text{(a.u.)}", ylabelsize = 35, titlesize = 30,
        yticklabelsize = 30, yaxisposition = :right, rightspinecolor = ax2_color, 
        yticklabelcolor = ax2_color, ylabelcolor = ax2_color, ytickcolor = ax2_color)
        hidexdecorations!(ax2)
        if is_at
            ylims!(ax1, -2.8, 0.0)
            ylims!(ax2, -7.5, 9.0) # for at
        else
            ylims!(ax1, -4.2, 0.0)
            ylims!(ax2, -9.5, 13.0) # for gc    
        end
            ## exp vals
        lines!(ax1, t/1e03, x, color = :royalblue1, linewidth = 4)
        lines!(ax2, t/1e03, p, color = :crimson, linewidth = 4)
        nonclass_filename = is_at ? "moyal_solver/graphics/AT/nonclassicality_exp.pdf" : "moyal_solver/graphics/GC/nonclassicality_exp.pdf"
    end
    ## ==== delta 
    lines!(ax, t/1e03, delta, linewidth = 6, color = :purple)
    hlines!(ax, [0], color = :black, linestyle = :dash, alpha = 0.5)
    
    # display(fig)
    
    nonclass_dir = dirname(nonclass_filename)
    if !isdir(nonclass_dir)
        mkpath(nonclass_dir)
        println("Created directory: $nonclass_dir")
    end
    save(nonclass_filename, fig)
    return fig
end
create_nonclassicality_plot(;is_at = true, isnot_sole_nonclassicality = true)
##
alpha = 1.963
a4 = 0.0207
a3 = -0.0053
a2 = -0.0414
a1 = 0.0158
a0 = 0.0312
  

function plot_Godbeer_AT_potential()
    x_unique = range(-3.25, 2.9, 300)
    p_unique = range(-15.5, 15.5, 200)
    # v_unique = x_unique./alpha
    
    t = 0.0  
    # Vx = @. a4 * v_unique^4 + a3 * v_unique^3 + a2 * v_unique^2 + a1 * v_unique + a0 
    m = 1836
    H = [(p^2)/(2m) + model_at(x) for p in p_unique, x in x_unique]
    
    Emin = 0
    Emax = 0 + 0.1  
    levels = range(Emin, Emax, length=15)
    
    fig = Figure(size=(1000, 400))
    
    ax1 = Axis(fig[1,1], 
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"U^{\mathrm{A-T}}(x) \; \text{(a.u.)}",
               title = L"\text{Przybliżenie harmoniczne A-T}",
               xlabelsize = 27, 
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20,
               limits = ((-3.25, 2.9), (-0.001, 0.035)))
    lines!(ax1, x_unique, model_at, linewidth=3, color=:blue)
    
    ax2 = Axis(fig[1,2], 
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"$p$ \text{ (a.u.)}",
               title = L"\text{Hamiltonian w przestrzeni fazowej}",
               xlabelsize = 27,limits = ((-2.95, 2.8), (-15.5, 15.5)),
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.5)
    
# bands symbolizing Ωs
    Ω_C = -3.3:0.1:-0.5
    Ω_B = -0.5:0.1:1.2
    Ω_T = 1.2:0.1:2.9
    lines!(ax1, Ω_C, [-0.0005], color = :lightskyblue, linewidth = 8.)
    lines!(ax1, Ω_T, [-0.0005], color = :springgreen4, linewidth = 8.)
    lines!(ax1, Ω_B, [-0.0005], color = :gray60, linewidth = 8.)
    
    text!(ax1, -1.9, 0.003; text =  L"\Omega_{\mathrm{C}}", fontsize = 25)
    text!(ax1, 0.2, 0.003; text =  L"\Omega_{\mathrm{B}}", fontsize = 25)
    text!(ax1, 1.85, 0.003; text =  L"\Omega_{\mathrm{T}}", fontsize = 25)
    
    display(fig)
    save("moyal_solver/graphics/hamiltonian_model_at.pdf", fig)
    return fig
end

plot_Godbeer_AT_potential()

##
function plot_Slocombe_GC_potential()
    V1 = 0.1617
    V2 = 0.082
    a1 = 0.305
    a2 = 0.755
    r1 = -2.7
    r2 = 2.1
    m = 1836
    x_unique = range(-4.4, 2.6, 300)
    p_unique = range(-15.5, 15.5, 200)
    
    t = 0.0  
    Vx = @. V1 * (exp(-2 * a1 * (x_unique - r1)) - 2 * exp(-a1 * (x_unique - r1))) + V2 * (exp(-2 * a2 * (r2 - x_unique)) - 2 * exp(-a2 * (r2 - x_unique))) + 0.166 + 0.00019
    H = [(p^2)/(2m) + model_gc(x) for p in p_unique, x in x_unique]
    
    Emin = 0
    Emax = 0 + 0.1  
    levels = range(Emin, Emax, length=15)
    
    fig = Figure(size=(1000, 400))
    
    ax1 = Axis(fig[1,1], 
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"U^{\mathrm{G-C}}(x) \; \text{(a.u.)}",
               title = L"\text{Przybliżenie harmoniczne G-C}",
               xlabelsize = 27,
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20,
               limits = ((-4.4, 2.6), (-0.001, 0.026)))
    lines!(ax1, x_unique, model_gc, linewidth=3, color=:blue)
    # text!(ax1, -3.7, 0.024; text = L"\text{(a)}", fontsize = 30)

    ax2 = Axis(fig[1,2], 
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"$p$ \text{ (a.u.)}",
               title = L"\text{Hamiltonian w przestrzeni fazowej}",
               xlabelsize = 27,limits = ((-3.8, 2.4), (-15.5, 15.5)),
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.5)
# bands symbolizing Ωs
    Ω_C = -4.4:0.1:-1.0
    Ω_B = -1.0:0.1:1.11
    Ω_T = 1.1:0.1:2.7
    lines!(ax1, Ω_C, [-0.0005], color = :orange, linewidth = 8.)
    lines!(ax1, Ω_T, [-0.0005], color = :plum2, linewidth = 8.)
    lines!(ax1, Ω_B, [-0.0005], color = :gray60, linewidth = 8.)
    
    text!(ax1, -2.7, 0.0015; text =  L"\Omega_{\mathrm{C}}", fontsize = 25)
    text!(ax1, -0.1, 0.0015; text =  L"\Omega_{\mathrm{B}}", fontsize = 25)
    text!(ax1, 1.7, 0.0015; text =  L"\Omega_{\mathrm{T}}", fontsize = 25)
    
    # display(fig)
    save("moyal_solver/graphics/hamiltonian_model_gc.pdf", fig)
    return fig
end
plot_Slocombe_GC_potential()
## exp values

function get_exp_vals()
    data = readdlm("moyal_solver/build/output/stats.dat", skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]
    
    fig = Figure(size = (1000, 500))
    ax1_color = :royalblue1
    ax = Axis(fig[1,1], xlabel = L"t\; \;(10^3\text{ a.u.})", ylabel = L"\langle x\rangle\; (\text{a.u.})", 
    xlabelsize = 40, ylabelsize = 40, xticklabelsize = 30, yticklabelsize = 30,
    leftspinecolor = ax1_color, yaxisposition = :left, 
    yticklabelcolor = ax1_color, ylabelcolor = ax1_color, ytickcolor = ax1_color)
    ax2_color = :crimson
    ax2 = Axis(fig[1,1], ylabel = L"\langle p \rangle\; \text{(a.u.)}", ylabelsize = 40, titlesize = 30,
    yticklabelsize = 30, yaxisposition = :right, rightspinecolor = ax2_color, 
    yticklabelcolor = ax2_color, ylabelcolor = ax2_color, ytickcolor = ax2_color)
    hidexdecorations!(ax2)
    # ylims!(ax2, 8, 10)
    
    lines!(ax, t/1e03, x, color = ax1_color, linewidth = 4)
    lines!(ax2, t/1e03, p, color = ax2_color, linewidth = 4)
    
    display(fig)
    # save("moyal_solver/graphics/GC/exp_vals.pdf", fig)
end

get_exp_vals()

## create trajectory of xp values
function get_traj_of_exp_vals(is_at::Bool = true)
    if is_at
        m = 1836
        x_unique = range(-3.8, 2.8, 300)
        p_unique = range(-12.5, 12.5, 200)
        
        t = 0.0  
        Emin = 0
        Emax = 0 + 0.1  
        levels = range(Emin, Emax, length=20)
        
    else
        x_unique = range(-4.15, 2.8, 300)
        p_unique = range(-12.5, 12.5, 200)
        
        t = 0.0  
        m = 1836
        
        Emin = 0
        Emax = 0 + 0.1  
        levels = range(Emin, Emax, length=20)
    end
    H = is_at ? [(p^2)/(2m) + model_at(x) for p in p_unique, x in x_unique] : [(p^2)/(2m) + model_gc(x) for p in p_unique, x in x_unique]
    data = is_at ? readdlm("moyal_solver/build_at_model/output/stats.dat", skipstart = 1) : readdlm("moyal_solver/build_gc_model/output/stats.dat", skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]
    
    contour_cmap = :viridis                  
    trajectory_color = (:midnightblue, 0.9)   
    start_point_color = :coral1              
    # -------------------------------

    fig = Figure(size=(800, 500))
    title = is_at ? L"\text{A-T}" : L"\text{G-C}"
    ax2 = Axis(fig[1,1], 
               xlabel = L"x \; \text{(a.u.)}", 
               ylabel = L"p \; (\text{a.u.})",
               title = title,
               xlabelsize = 30,
               ylabelsize = 30, titlesize = 30,
               xticklabelsize = 24, yticklabelsize = 24,
               limits = is_at ? ((-2.95, 2.8), (-12.5, 12.5)) :  ((-4.1, 2.4), (-12.5, 12.5))
               )
               
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.2, alpha = 0.4, colormap = contour_cmap)
    
    lines!(ax2, x, p, label = "trajektoria wartości oczekiwanych", color = trajectory_color, linewidth = 3.5)
    
    scatter!(ax2, [-1.1], [5.4], color=start_point_color, markersize=15, label="centrum początkowego gaussianu")
    
    Legend(fig[2,1], ax2, position=:lb, framevisible = false, labelsize = 24, orientation = :horizontal, nbanks = 2)
    
    display(fig)
    if is_at
        save("moyal_solver/graphics/AT/trajectory.pdf", fig)
    else
        save("moyal_solver/graphics/GC/trajectory.pdf", fig)
    end
    return fig
end
get_traj_of_exp_vals(true)
##
function get_traj_harmonic_oscillator()
    x_unique = range(-5.2, 5.2, 200)
    p_unique = range(-5.2, 5.2, 200)
    m = 1
    H = [(p.^2)/(2m) + 0.5 * x.^2 for p in p_unique, x in x_unique] 
    data = readdlm("moyal_solver/build2_HO/output/stats.dat", skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]
    Emin = 0
    Emax = 0 + 7.
    levels = range(Emin, Emax, length=15)
    contour_cmap = :viridis                  
    trajectory_color = (:midnightblue, 0.9)   
    start_point_color = :coral1  
    fig = Figure(size=(800, 500))
    ax2 = Axis(fig[1,1], 
               xlabel = L"x \; \text{(a.u.)}", 
               ylabel = L"p \; (\text{a.u.})",
            #    title = title,
               xlabelsize = 25,
               ylabelsize = 25, titlesize = 25,
               xticklabelsize = 20, yticklabelsize = 20,
               limits = ((-2.7, 2.7), (-2.7, 2.7))
               )
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.5, alpha = 0.5, colormap = contour_cmap)
    lines!(ax2, x, p, label = "trajekroria wartości oczekiwanych", color = trajectory_color, linewidth = 3.5)
    
    scatter!(ax2, [-1.0], [1.], color=start_point_color, markersize=15, label="centrum początkowego gaussianu")
    Legend(fig[2,1], ax2, position=:lb, framevisible = false, labelsize = 24, orientation = :horizontal, nbanks = 2)
    save("moyal_solver/graphics/HO_trajectory.pdf", fig)
    return fig
end
get_traj_harmonic_oscillator()
## 3 D WDF 

using GLMakie
using DelimitedFiles

function plot_wigner_3d(x_coords, p_coords, W;
                        scaling_factor=0.3,
                        colorrange=(-0.32, 0.32))
    
    α = scaling_factor
    W_scaled = sign.(W) .* asinh.(abs.(W) / α) * α

    GLMakie.activate!()
    
    β = 0.01
    wigner_scale = ReversibleScale(
        w ->  asinh(w / β) / log(10),   
        w -> β * sinh(log(10) * w)      
    )

    fig = Figure(size=(1200, 800))
    
    # Definiujemy limity jawnie, aby użyć ich do pozycjonowania ścianek
    x_lims = (-3.95, 2.45)
    p_lims = (-14.5, 14.5)
    z_lims = (-0.32, 0.32)

    # x_lims = (-2.95, 2.8)
    # p_lims = (-14.5, 14.5)
    # z_lims = (-0.32, 0.32)

    ax3d = Axis3(fig[1, 1],
                 xlabel=L"x \; (\text{a.u.})",
                 ylabel=L"p \; (\text{a.u.})",
                 zlabel=L"\varrho(x,p) \; (\text{a.u.})",
                 xticklabelsize = 25, yticklabelsize = 25, zticklabelsize = 25,
                 title=L"t = 9.25 \times 10^{4} \; \text{a.u.}",
                 titlesize=30,
                 zlabeloffset = 60,
                 xlabelsize=30,
                 ylabelsize=30,
                 zlabelsize=30,
                 limits = (x_lims, p_lims, z_lims),
                 )
    
    # Główny wykres powierzchniowy funkcji Wignera
    surf = surface!(ax3d, x_coords, p_coords, W_scaled,
                   colormap=:RdBu,
                   colorscale = wigner_scale,
                   colorrange=colorrange,
                   shading=true,
                   transparency=true,
                   alpha=0.8)
    
    # Mapowanie płaskie (podłoga)
    z_offset = z_lims[1]
    Z_flat = fill(z_offset, length(x_coords), length(p_coords))
    heatmap_surf = surface!(ax3d, x_coords, p_coords, Z_flat,
                           color=W_scaled,
                           colorscale = wigner_scale,
                           colormap=:RdBu,
                           colorrange=colorrange,
                           shading=false,
                           transparency=false)
    
    # --- OBLICZANIE ROZKŁADÓW BRZEGOWYCH (Całkowanie numeryczne) ---
    # KROKI siatki (zakładamy równomierną siatkę)
    dx = x_coords[2] - x_coords[1]
    dp = p_coords[2] - p_coords[1]
    
    # P(x) = \int W(x,p) dp (suma wzdłuż kolumn/wierszy pędu)
    # W ma wymiary (nx, np), więc sumujemy po drugim wymiarze (p)
    P_x = sum(W, dims=2)[:] .* dp
    
    # P(p) = \int W(x,p) dx
    # Sumujemy po pierwszym wymiarze (x)
    P_p = sum(W, dims=1)[:] .* dx

    max_z_display = 0.25
    P_x_plot = (P_x ./ maximum(P_x)) .* max_z_display .+ z_lims[1]
    P_p_plot = (P_p ./ maximum(P_p)) .* max_z_display .+ z_lims[1]

    # --- RYSOWANIE NA ŚCIANKACH ---
    # 1. Rozkład położenia P(x) na tylnej ściance (p = p_lims[2])
    p_wall = fill(p_lims[2], length(x_coords))
    lines!(ax3d, x_coords, p_wall, P_x_plot.+0.32, 
           color=:black, linewidth=4, label=L"P(x)")
    
    # 2. Rozkład pędu P(p) na lewej ściance (x = x_lims[1])
    x_wall = fill(x_lims[2], length(p_coords))
    lines!(ax3d, x_wall, p_coords, P_p_plot.+0.32, 
           color=:darkred, linewidth=4, label=L"P(p)")
    
    # ---
    Colorbar(fig[1, 2], surf;
             label=L"\varrho(x,p; t)",
             labelsize=30, ticklabelsize = 25)

    return fig, ax3d
end

function load_and_plot_wigner_3d(filename; kwargs...)
    data = readdlm(filename)
    x_coords = data[:, 1]
    p_coords = data[:, 2] 
    wigner_real = data[:, 3]
    
    x_unique = sort(unique(x_coords))
    p_unique = sort(unique(p_coords))
    nx, np = length(x_unique), length(p_unique)
    
    W = reshape(wigner_real, np, nx)'

    return plot_wigner_3d(x_unique, p_unique, W;
                          kwargs...)
end

fig, ax3d = load_and_plot_wigner_3d("moyal_solver/build_gc_model/output/wigner_00000096000.dat")
# fig, ax3d = load_and_plot_wigner_3d("moyal_solver/build_godbeer_higher_p0/output/wigner_00000096200.dat")
display(fig)
save("moyal_solver/graphics/GC/WDF_3d.png", fig)

## function that determines what percent of the WDF is in the canonical form / barrier / tautomerical form
function get_percent_of_WDF(; is_at::Bool = true, is_tau_init::Bool = false)
    if is_tau_init
        output_dir = is_at ?
            "moyal_solver/build_at_model_tautomeric/output/" :
            "moyal_solver/build_gc_model_tautomeric/output/"
    else
    output_dir = is_at ?
        "moyal_solver/build_at_model/output/" :
        "moyal_solver/build_gc_model/output/"
    end
    isdir(output_dir) || error("Output directory not found: $output_dir")

    wigner_files = filter(
        f -> startswith(f, "wigner_") && endswith(f, ".dat"),
        readdir(output_dir)
    )

    sort!(wigner_files)

    # take every 50th snapshot
    snap4 = wigner_files[1:1:end]

    # ============================================================
    # Read grid once
    # ============================================================

    data = readdlm(joinpath(output_dir, wigner_files[1]))

    x_unique = sort(unique(data[:, 1]))
    p_unique = sort(unique(data[:, 2]))

    nx = length(x_unique)
    np = length(p_unique)

    dx = x_unique[2] - x_unique[1]
    dp = p_unique[2] - p_unique[1]

    dA = dx * dp

    # ============================================================
    # Region masks
    # ============================================================

    can_mask = is_at ? x_unique .< -0.51 : x_unique .< -1.03
    bar_mask = is_at ? (-0.51 .<= x_unique) .& (x_unique .< 1.21) : (-1.03 .<= x_unique) .& (x_unique .< 1.15)
    tau_mask = is_at ? x_unique .>= 1.21 : x_unique .>= 1.15

    # ============================================================
    # Arrays
    # ============================================================

    # signed quasi-probability
    P_can = Float64[]
    P_bar = Float64[]
    P_tau = Float64[]

    # absolute localization
    L_can = Float64[]
    L_bar = Float64[]
    L_tau = Float64[]

    # negativity
    Neg_total = Float64[]

    # ============================================================
    # Main loop
    # ============================================================

    for filename in snap4

        raw = readdlm(joinpath(output_dir, filename))

        W = reshape(raw[:, 3], np, nx)'

        # --------------------------------------------------------
        # Signed populations
        # --------------------------------------------------------

        Pcan = sum(W[can_mask, :]) * dA
        Pbar = sum(W[bar_mask, :]) * dA
        Ptau = sum(W[tau_mask, :]) * dA

        push!(P_can, Pcan)
        push!(P_bar, Pbar)
        push!(P_tau, Ptau)

        # --------------------------------------------------------
        # Absolute localization
        # --------------------------------------------------------

        abs_total = sum(abs.(W)) * dA

        Lcan = sum(abs.(W[can_mask, :])) * dA / abs_total
        Lbar = sum(abs.(W[bar_mask, :])) * dA / abs_total
        Ltau = sum(abs.(W[tau_mask, :])) * dA / abs_total

        push!(L_can, Lcan)
        push!(L_bar, Lbar)
        push!(L_tau, Ltau)

        # --------------------------------------------------------
        # Negativity volume
        # --------------------------------------------------------

        neg = sum(abs.(W) .- W) * dA / 2

        push!(Neg_total, neg)
    end

    # ============================================================
    # Time axis
    # ============================================================

    t = collect(1:length(snap4))./20

    # ============================================================
    # Figure 1 — signed quasi-probabilities
    # ============================================================

    fig1 = Figure(size = (1000, 500))

    ax1 = Axis(
        fig1[1,1],
        xlabel = L"t \; (10^3\ \mathrm{a.u.})",
        ylabel = L"\mathcal{P}_{\Omega}(t)",
        xlabelsize = 32,
        ylabelsize = 32,
        xticklabelsize = 24,
        yticklabelsize = 24
    )

    lines!(ax1, t, P_can, label = L"\Omega_{\mathrm{C}}", linewidth = 5.0, color = is_at ? :lightskyblue : :orange)
    lines!(ax1, t, P_bar, label = L"\Omega_{\mathrm{B}}", linewidth = 5.0, color = :gray)
    lines!(ax1, t, P_tau, label = L"\Omega_{\mathrm{T}}", linewidth = 5.0, color = is_at ? :springgreen4 : :plum2)

    Legend(fig1[2,1], ax1, labelsize = 32, orientation = :horizontal, framevisible = false)

    # ============================================================
    # Figure 2 — localization measure
    # ============================================================

    fig2 = Figure(size = (1000, 500))

    ax2 = Axis(
        fig2[1,1],
        xlabel = L"t \; (10^3\ \mathrm{a.u.})",
        ylabel = L"\mathcal{P}_{\Omega}(t)",
        xlabelsize = 28,
        ylabelsize = 28,
        xticklabelsize = 22,
        yticklabelsize = 22
    )

    lines!(ax2, t, L_can, label = L"\Omega_{\mathrm{C}}", linewidth = 5.0, color = :limegreen)
    lines!(ax2, t, L_bar, label = L"\Omega_{\mathrm{B}}", linewidth = 5.0, color = :deepskyblue3)
    lines!(ax2, t, L_tau, label = L"\Omega_{\mathrm{T}}", linewidth = 5.0, color = :palevioletred3)

    Legend(fig2[2,1], ax2, labelsize = 32, orientation = :horizontal, framevisible = false)

    # ============================================================
    # Figure 3 — negativity
    # ============================================================

    fig3 = Figure(size = (1000, 500))

    ax3 = Axis(
        fig3[1,1],
        xlabel = L"t \; (10^3\ \mathrm{a.u.})",
        ylabel = L"\mathcal N",
        xlabelsize = 28,
        ylabelsize = 28,
        xticklabelsize = 22,
        yticklabelsize = 22
    )

    lines!(ax3, t, Neg_total, label = "nieklasyczność")

    Legend(fig3[2,1], ax3, labelsize = 28, orientation = :horizontal)

    if is_at
        if is_tau_init
            save("moyal_solver/graphics/AT/fraction_wdf_not_normalized_tau.pdf", fig1)
            save("moyal_solver/graphics/AT/fraction_wdf_tau.pdf", fig2)
            save("moyal_solver/graphics/AT/negativity_wdf_tau.pdf", fig3)
        else
            save("moyal_solver/graphics/AT/fraction_wdf_not_normalized.pdf", fig1)
            save("moyal_solver/graphics/AT/fraction_wdf.pdf", fig2)
            save("moyal_solver/graphics/AT/negativity_wdf.pdf", fig3)    
        end
    else
        if is_tau_init
            save("moyal_solver/graphics/GC/fraction_wdf_not_normalized_tau.pdf", fig1)
            save("moyal_solver/graphics/GC/fraction_wdf_tau.pdf", fig2)
            save("moyal_solver/graphics/GC/negativity_wdf_tau.pdf", fig3)
        else
            save("moyal_solver/graphics/GC/fraction_wdf_not_normalized.pdf", fig1)
            save("moyal_solver/graphics/GC/fraction_wdf.pdf", fig2)
            save("moyal_solver/graphics/GC/negativity_wdf.pdf", fig3)    
        end
    end
    return (
        P_can, P_bar, P_tau,
        L_can, L_bar, L_tau,
        Neg_total
    )
end
get_percent_of_WDF(is_at = false, is_tau_init = true)

## FFT delta
using FFTW
function find_spectral_peaks(freqs, power;
        n_peaks        = 3,
        min_freq       = nothing,
        min_prominence = 0.05)

    # Próg dolny: jeśli nie podano, pomijamy DC (f = 0)
    f_min = isnothing(min_freq) ? freqs[2] : min_freq
    mask  = freqs .>= f_min

    f_masked = freqs[mask]
    p_masked = power[mask]

    threshold = min_prominence * maximum(p_masked)

    # Wykrywanie lokalnych maksimów (sąsiad z obu stron mniejszy)
    is_peak = falses(length(p_masked))
    for i in 2:(length(p_masked) - 1)
        if p_masked[i] > p_masked[i-1] &&
           p_masked[i] > p_masked[i+1] &&
           p_masked[i] >= threshold
            is_peak[i] = true
        end
    end

    peak_freqs  = f_masked[is_peak]
    peak_powers = p_masked[is_peak]

    # Sortuj malejąco po amplitudzie, zwróć n_peaks pierwszych
    order       = sortperm(peak_powers; rev = true)
    n           = min(n_peaks, length(order))

    return peak_freqs[order[1:n]], peak_powers[order[1:n]]
end

function get_fft_delta(; is_at::Bool = true, isnot_sole_nonclassicality::Bool = false)
    stats_file = is_at ? "moyal_solver/build_at_model/output/stats.dat" :
                         "moyal_solver/build_gc_model/output/stats.dat"
    if !isfile(stats_file)
        println("Stats file not found, skipping nonclassicality plot")
        return nothing
    end

    nonclass_filename = is_at ? "moyal_solver/graphics/AT/nonclassicality_fft.pdf" :
                                "moyal_solver/graphics/GC/nonclassicality_fftw.pdf"

    data = readdlm(stats_file, skipstart = 1)
    t     = data[40:end, 2]
    x     = data[:, 3]
    p     = data[:, 4]
    delta = size(data, 2) >= 8 ? data[40:end, 8] : zeros(length(t))

    # ======== FFT =====================================

    N  = length(t)
    dt = (t[end] - t[1]) / (N - 1)        # zakładamy równomierne próbkowanie

    # Transformata i odpowiadające jej częstotliwości
    delta_fft  = fft(delta)                # zespolone widmo
    freqs      = fftfreq(N, 1 / dt)       # wektor częstotliwości [1/jednostka czasu]

    # Jednostronne widmo amplitudowe (f >= 0)
    half      = freqs .>= 0         
    freqs_pos = freqs[half]
    power     = abs.(delta_fft[half])
    power[2:end] .*= 2               # kompensacja jednostronności

    _, peak_idx = findmax(power[2:end])
    peak_freq   = freqs_pos[peak_idx + 1]
    peak_freqs, peak_powers = find_spectral_peaks(freqs_pos, power;
        n_peaks        = 3,
        min_prominence = 0.05,
    )

    fig = Figure(size = (900, 650))

    ax1 = Axis(fig[1, 1];
        xlabel = L"t \; \text{(a.u.)}",
        ylabel = L"\delta(t)",
        title  = is_at ? L"\text{Parametr nieklasyczności model AT}" : "Nieklasyczność – model GC",
        titlesize  = 20,
        xlabelsize = 20,
        ylabelsize = 20,
    )
    lines!(ax1, t, delta; color = :steelblue, linewidth = 1.5)

    ax2 = Axis(fig[2, 1];
        xlabel = L"\omega \; \text{(a.u.)}",
        ylabel = L"|\mathcal{F}_{t\to\omega}(\delta)|",
        title  = L"\text{Widmo amplitudowe} \; \delta",
        titlesize  = 20,
        xlabelsize = 20,
        ylabelsize = 20,
    )
    xlims!(ax2, 0.0, 0.01)
    # ylims!(ax2, 0.0, 30)
    lines!(ax2, freqs_pos, power; color = :crimson, linewidth = 1.5)

    vlines!(ax2, [peak_freq]; color = :orange, linestyle = :dash, linewidth = 1.5,
            label = "f_peak = $(round(peak_freq, digits=4))")
    colors = [:orange, :green, :purple]
    for (i, (f, a)) in enumerate(zip(peak_freqs, peak_powers))
        vlines!(ax2, [f];
            color     = colors[i],
            linestyle = :dash,
            linewidth = 1.5,
            label     = "f$(i) = $(round(f, digits=4))",
        )
        # Etykieta nad pikiem
        text!(ax2, f, a;
            text   = " f$(i)",
            color  = colors[i],
            fontsize = 11,
        )
    end
    axislegend(ax2; position = :rt, labelsize = 11)

    ax2.yscale = log10

    save(nonclass_filename, fig)
    # display(fig)
    println("Zapisano: $nonclass_filename")
    return fig
end
get_fft_delta()