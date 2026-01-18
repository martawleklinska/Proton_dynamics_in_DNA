using CairoMakie, DelimitedFiles, Printf
using StatsBase

const dt = 1.  # time step 

function create_wigner_animation()
    output_paths = [
        "build/output/",  
        "masters/moyal_solver/build/output/",               
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
    
    if output_dir === nothing
        error("Output directory not found! Run the C++ program first.")
    end

    wigner_files = filter(f -> startswith(f, "wigner_") && endswith(f, ".dat"), 
                         readdir(output_dir))
    sort!(wigner_files)
    
    if isempty(wigner_files)
        error("No Wigner function files found in $output_dir")
    end
    
    println("Found $(length(wigner_files)) Wigner function files")
    
    first_file = joinpath(output_dir, wigner_files[1])
    data = readdlm(first_file)
    x_coords, p_coords = data[:, 1], data[:, 2]
    
    x_unique = sort(unique(x_coords))
    p_unique = sort(unique(p_coords))
    nx, np = length(x_unique), length(p_unique)
    
    println("=== Grid size: $nx × $np")
    println("=== x range: $(round(minimum(x_unique), digits=2)) to $(round(maximum(x_unique), digits=2))")
    println("=== p range: $(round(minimum(p_unique), digits=2)) to $(round(maximum(p_unique), digits=2))")
    
    println("Calculating Wigner function range...")
    w_min, w_max = Inf, -Inf
    sample_files = wigner_files[1:max(1, length(wigner_files)÷10):end]  
    
    for filename in sample_files
        data = readdlm(joinpath(output_dir, filename))
        w_vals = data[:, 3]
        w_min = min(w_min, minimum(w_vals))
        w_max = max(w_max, maximum(w_vals))
    end
    println("Wigner range: $(round(w_min, digits=4)) to $(round(w_max, digits=4))")
    
    animation_files = wigner_files[1:4:end]
    n_frames = length(animation_files)
        
    fig = Figure(size = (900, 700), fontsize = 16)
    
    time_obs = Observable("t = 0.0")
    wigner_obs = Observable(zeros(nx, np))
    
    ax = Axis(fig[1, 1], 
              xlabel = L"\text{Położenie}\; x", 
              ylabel = L"\text{Pęd} \;p",
              title = time_obs,
              titlesize = 25,
              limits = ((-10., 10.), (-20., 20.)),
              xlabelsize = 25,
              ylabelsize = 25)
    
    hm = heatmap!(ax, x_unique, p_unique, wigner_obs,
                      colormap = :RdBu,
                      colorrange = (w_min, w_max))

    try
        Colorbar(fig[1, 2], hm, label = L"\varrho(x,p,t)", labelsize = 25)
    catch e
        println("Warning: Could not create colorbar: $e")
    end
    
    gif_filename = "moyal_solver/graphics/wigner_evolution_GC.gif"
    
    gif_dir = dirname(gif_filename)
    if !isdir(gif_dir)
        mkpath(gif_dir)
        println("Created directory: $gif_dir")
    end
    
    record(fig, gif_filename, 1:n_frames; framerate = 8) do frame_idx
        filename = animation_files[frame_idx]
        
        try
            step_str = match(r"wigner_(\d+)\.dat", filename).captures[1]
            step = parse(Int, step_str)
            time_val = step * dt
            
            data = readdlm(joinpath(output_dir, filename))
            wigner_real = data[:, 3]
            
            W = reshape(wigner_real, np, nx)'
            W_vis = sign.(W) .* abs.(W).^(1/3)  
            
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
    
    snapshot_files = wigner_files[1:max(1, length(wigner_files)÷10):end]  
    
    for (i, filename) in enumerate(snapshot_files)
        step_str = match(r"wigner_(\d+)\.dat", filename).captures[1]
        step = parse(Int, step_str)
        time_val = step * dt
        
        data = readdlm(joinpath(output_dir, filename))
        wigner_real = data[:, 3]
        
        W = reshape(wigner_real, np, nx)'
        W_vis = sign.(W) .* abs.(W).^(1/3)
        
        fig_snap = Figure(size = (1000, 600))
        ax_snap = Axis(fig_snap[1, 1],
                      xlabel = L"\text{Położenie } x",
                      ylabel = L"\text{Pęd } p", 
                      title = @sprintf("Funkcja Wignera (t = %.3f)", time_val),
                      titlesize = 35,
                      limits = ((-10., 10.), (-20., 20.)),
                      xlabelsize = 35,
                      ylabelsize = 35)
        
        hm_snap = try
            heatmap!(ax_snap, x_unique, p_unique, W_vis,
                    colormap = :RdBu,
                    colorrange = (w_min, w_max))
        catch e
            println("Warning: Using fallback colormap for snapshot $i: $e")
            heatmap!(ax_snap, x_unique, p_unique, W_vis,
                    colormap = :RuBu,
                    colorrange = (w_min, w_max))
        end
        
        try
            Colorbar(fig_snap[1, 2], hm_snap, label = L"\varrho(x,p,t)", labelsize = 35)
        catch e
            println("Warning: Could not create colorbar for snapshot $i: $e")
        end
        
        png_filename = @sprintf("moyal_solver/graphics/GC/wigner_snapshot_t%.3f.png", time_val)
        

        png_dir = dirname(png_filename)
        if !isdir(png_dir)
            mkpath(png_dir)
            println("Created directory: $png_dir")
        end
        
        save(png_filename, fig_snap, px_per_unit = 2)  

    end
    
    println("\nVisualization complete!")
    return fig, gif_filename
end

create_wigner_animation()
##

function create_nonclassicality_plot()
    stats_file = "masters/moyal_solver/build/output/stats.dat"
    if !isfile(stats_file)
        println("Stats file not found, skipping nonclassicality plot")
        return nothing
    end
    
    data = readdlm(stats_file, skipstart = 1)
    t = data[:, 2]
    delta = size(data, 2) >= 8 ? data[:, 8] : zeros(length(t))  
    
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1],
              xlabel = L"\text{Czas } t",
              ylabel = L"\text{Parametr nieklasyczności } \delta(t)",
              title = L"\text{Ewolucja nieklasyczności stanu kwantowego}",
              titlesize = 25,
                      xlabelsize = 25,
                      ylabelsize = 25, xticklabelsize = 20, yticklabelsize = 20)
    
    lines!(ax, t, delta, linewidth = 3, color = :purple)
    hlines!(ax, [0], color = :black, linestyle = :dash, alpha = 0.5)
    
    # Create directory if it doesn't exist
    nonclass_filename = "moyal_solver/graphics/nonclassicality.png"
    nonclass_dir = dirname(nonclass_filename)
    if !isdir(nonclass_dir)
        mkpath(nonclass_dir)
        println("Created directory: $nonclass_dir")
    end
    
    save(nonclass_filename, fig)
    return fig
end
create_nonclassicality_plot()
##
alpha = 1.963
a4 = 0.0207
a3 = -0.0053
a2 = -0.0414
a1 = 0.0158
a0 = 0.0312
        # double alpha = 1.963;
        # double a4 = 0.0207;
        # double a3 = -0.0053;
        # double a2 = -0.0414;
        # double a1 = 0.0158;
        # double a0 = 0.0312;

        # double operator()(double v, double t = 0.0) const override {
        #     double x = v/alpha;
        #     return a4*x * x * x * x + a3*x * x * x + a2*x * x + a1*x + a0;
        # }

function plot_duffing_potential()
    x_unique = range(-3, 3, 300)
    p_unique = range(-15.5, 15.5, 200)
    v_unique = x_unique./alpha
    
    t = 0.0  
    Vx = @. a4 * v_unique^4 + a3 * v_unique^3 + a2 * v_unique^2 + a1 * v_unique + a0 
    m = 1836
    H = [(p^2)/(2m) + V for p in p_unique, V in Vx]
    
    Emin = 0
    Emax = 0 + 0.1  
    levels = range(Emin, Emax, length=25)
    
    fig = Figure(size=(1000, 400))
    
    ax1 = Axis(fig[1,1], 
               xlabel = L"x", 
               ylabel = L"V(x)",
               title = L"\text{Potencjał 4 stopnia}",
               xlabelsize = 20,
               ylabelsize = 20, titlesize = 18)
    lines!(ax1, x_unique, Vx, linewidth=3, color=:blue)

    ax2 = Axis(fig[1,2], 
               xlabel = L"x", 
               ylabel = L"p",
               title = L"\text{Hamiltonian w przestrzeni fazowej}",
               xlabelsize = 20,
               ylabelsize = 20, titlesize = 18)
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.5)
    
    # scatter!(ax2, [x_min, -x_min], [0, 0], color=:red, markersize=12, label="studnie potencjału")
    
    # scatter!(ax2, [-4.0], [1.15], color=:green, markersize=15, label="warunek początkowy")
    # axislegend(ax2, position=:lb, framevisible = false)
    

    # display(fig)
    save("moyal_solver/graphics/hamiltonian_godbeer.pdf", fig)
    return fig
end

plot_duffing_potential()

## exp values

function get_exp_vals()
    data = readdlm("moyal_solver/build/output/stats.dat", skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]
    
    fig = Figure(size = (1000, 500))
    ax1_color = :royalblue1
    ax = Axis(fig[1,1], xlabel = L"t\; (\text{a.u.})", ylabel = L"\langle x\rangle\; (\text{a.u.})", 
    xlabelsize = 40, ylabelsize = 40, xticklabelsize = 30, yticklabelsize = 30,
    leftspinecolor = ax1_color, yaxisposition = :left, 
    yticklabelcolor = ax1_color, ylabelcolor = ax1_color, ytickcolor = ax1_color)
    ax2_color = :crimson
    ax2 = Axis(fig[1,1], ylabel = L"\langle p \rangle\; \text{(a.u.)}", ylabelsize = 40, titlesize = 30,
    yticklabelsize = 30, yaxisposition = :right, rightspinecolor = ax2_color, 
    yticklabelcolor = ax2_color, ylabelcolor = ax2_color, ytickcolor = ax2_color)
    hidexdecorations!(ax2)
    # ylims!(ax2, 8, 10)
    
    lines!(ax, t, x, color = ax1_color, linewidth = 4)
    lines!(ax2, t, p, color = ax2_color, linewidth = 4)
    
    # display(fig)
    save("moyal_solver/graphics/xp_exp_val.pdf", fig)
end

get_exp_vals()
## create trajectory of xp values
function get_traj_of_exp_vals()
    alpha = -5e-1
    beta  =  1e-2
    gamma =  0.5
    omega =  1.0e-1
    m = 1
    x_unique = range(-10, 12, 300)
    p_unique = range(-5.5, 5.5, 200)
    
    t = 0.0  
    Vx = @. 0.5*alpha * x_unique^2 + 
            0.25*beta * x_unique^4 - 
            gamma * x_unique * cos(omega * t)
    
    H = [(p^2)/(2m) + V for p in p_unique, V in Vx]

    data = readdlm("moyal_solver/build/output/stats.dat", skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]

    x_min = sqrt(-alpha/beta)  
    V_min = 0.5*alpha*x_min^2 + 0.25*beta*x_min^4
    
    Emin = minimum(H)
    Emax = V_min + 10.0  
    levels = range(Emin, Emax, length=25)
    
    fig = Figure(size=(800, 500))
    
    ax2 = Axis(fig[1,1], 
               xlabel = L"x", 
               ylabel = L"p",
               title = L"\text{Hamiltonian w przestrzeni fazowej}",
               xlabelsize = 25,
               ylabelsize = 25, titlesize = 25,
               xticklabelsize = 20, yticklabelsize = 20)
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.5)
    lines!(ax2, x, p, label = "trajekroria wartości oczekiwanych")
    
    scatter!(ax2, [-4.0], [1.15], color=:green, markersize=15, label="warunek początkowy")
    axislegend(ax2, position=:lb, framevisible = false, labelsize = 20)
    

    # display(fig)
    save("moyal_solver/graphics/trajectory.pdf", fig)
    return fig
end
get_traj_of_exp_vals()
