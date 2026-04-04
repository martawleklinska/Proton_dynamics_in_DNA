using CairoMakie, DelimitedFiles, Printf
using StatsBase

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
function plot_wigner_snapshots(; is_harmonic::Bool=false, is_gc::Bool=false, is_at::Bool=true)

    output_dir = "moyal_solver/build/output/"
    isdir(output_dir) || error("Output directory not found: $output_dir")

    wigner_files = filter(f -> startswith(f, "wigner_") && endswith(f, ".dat"), 
                         readdir(output_dir))
    sort!(wigner_files)

    # ── wczytanie siatki z pierwszego pliku ──────────────────────────────────
    data      = readdlm(joinpath(output_dir, wigner_files[1]))
    x_unique  = sort(unique(data[:, 1]))
    p_unique  = sort(unique(data[:, 2]))
    nx, np    = length(x_unique), length(p_unique)

    # ── potencjał i poziomice Hamiltonianu ───────────────────────────────────
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
        levels = range(0, 0.1, length=25)
        label  = "GC"
    elseif is_at
        m  = 1836.0
        alpha = 1.963
        a4 = 0.0207
        a3 = -0.0053
        a2 = -0.0414
        a1 = 0.0158
        a0 = 0.0312
        v  = x_unique ./ alpha
        Vx = @. a4*v^4 + a3*v^3 + a2*v^2 + a1*v + a0
        levels = range(0, 0.1, length=21)
        label  = "AT"
    else
        error("Wybierz jeden potencjał: is_harmonic, is_gc lub is_at")
    end

    H = [(p^2)/(2m) + V for p in p_unique, V in Vx]

    # ── zakres kolorów z próbki plików ───────────────────────────────────────
    sample = wigner_files[1:max(1, length(wigner_files)÷10):end]
    w_min, w_max = mapreduce(
        f -> extrema(readdlm(joinpath(output_dir, f))[:, 3]),
        (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
        sample
    )
    println("Wigner range: $(round(w_min, digits=4)) to $(round(w_max, digits=4))")

    # ── cztery pierwsze snapshoty ─────────────────────────────────────────────
    wigner_files = filter(f -> startswith(f, "wigner_") && endswith(f, ".dat"), 
                         readdir(output_dir))
    sort!(wigner_files)
    snap4        = wigner_files[1:25:min(600, length(wigner_files))]
    # snap4        = wigner_files[1:200:min(800, length(wigner_files))]

    fig = Figure(size=(1000, 1200))
    β = 0.01
    wigner_scale = ReversibleScale(
        w ->  asinh(w / β) / log(10),   # forward:  W  → kolor
        w -> β * sinh(log(10) * w)       # inverse:  kolor → W (dla colorbar)
    )

    for (idx, filename) in enumerate(snap4)
        row, col = divrem(idx - 1, 2) .+ 1          # (1,1) (1,2) (2,1) (2,2)
        titles = [L"t = 0.0 \text{ a.u.}", L"t = 2\times 10^3 \text{ a.u.}", L"t = 4\times 10^3\text{ a.u.}", L"t = 6\times 10^3\text{ a.u.}", L"t = 8\times 10^3\text{ a.u.}", L"t = 10\times 10^3\text{ a.u.}"]
        # titles = [L"t = 0.0 \text{ a.u.}", L"t = 200.0 \text{ a.u.}", L"t = 400.0\text{ a.u.}", L"t = 600.0\text{ a.u.}"]
        
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
                  limits       = ((-3., 3.), (-12., 12.)))

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
    out_file = joinpath(out_dir, "wigner_snapshots_2x2_sqrt.pdf")
    # out_file = joinpath(out_dir, "wigner_snapshots_2x2_long_sim.pdf")
    save(out_file, fig; px_per_unit=2)
    println("Zapisano: $out_file")
end

plot_wigner_snapshots(; is_harmonic=false, is_gc=true, is_at=false)##
##
function create_nonclassicality_plot()
    stats_file = "moyal_solver/build/output/stats.dat"
    if !isfile(stats_file)
        println("Stats file not found, skipping nonclassicality plot")
        return nothing
    end
    
    data = readdlm(stats_file, skipstart = 1)
    t = data[:, 2]
    delta = size(data, 2) >= 8 ? data[:, 8] : zeros(length(t))  
    
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1],
              xlabel = L"t \;\; (10^3\text{ a.u.)} ",
              ylabel = L"\delta_{\mathrm{W}}(t)",
            #   title = L"\text{Ewolucja nieklasyczności stanu kwantowego}",
              titlesize = 30,
                      xlabelsize = 30,
                      ylabelsize = 30, xticklabelsize = 25, yticklabelsize = 25)
    
    lines!(ax, t/1e03, delta, linewidth = 6, color = :purple)
    hlines!(ax, [0], color = :black, linestyle = :dash, alpha = 0.5)
    
    # Create directory if it doesn't exist
    nonclass_filename = "moyal_solver/graphics/GC/nonclassicality.pdf"
    nonclass_dir = dirname(nonclass_filename)
    if !isdir(nonclass_dir)
        mkpath(nonclass_dir)
        println("Created directory: $nonclass_dir")
    end
    
    # save(nonclass_filename, fig)
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
  

function plot_Godbeer_AT_potential()
    x_unique = range(-2.95, 2.8, 300)
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
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"U^{\mathrm{A-T}}(x) \; \text{(a.u.)}",
               title = L"\text{Potencjał 4 stopnia}",
               xlabelsize = 27, 
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    lines!(ax1, x_unique, Vx, linewidth=3, color=:blue)
    text!(ax1, -2.8, 0.034; text = L"\text{(a)}", fontsize = 30)

    ax2 = Axis(fig[1,2], 
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"$p$ \text{ (a.u.)}",
               title = L"\text{Hamiltonian w przestrzeni fazowej}",
               xlabelsize = 27,limits = ((-2.95, 2.8), (-15.5, 15.5)),
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.5)
    text!(ax2, -2.8, 10.2; text = L"\text{(b)}", fontsize = 30)
    
    display(fig)
    save("moyal_solver/graphics/hamiltonian_godbeer.pdf", fig)
    return fig
end

# plot_Godbeer_AT_potential()
##

function plot_Slocombe_GC_potential()
    V1 = 0.1617
    V2 = 0.082
    a1 = 0.305
    a2 = 0.755
    r1 = -2.7
    r2 = 2.1
    m = 1836
    x_unique = range(-3.8, 2.4, 300)
    p_unique = range(-15.5, 15.5, 200)
    
    t = 0.0  
    Vx = @. V1 * (exp(-2 * a1 * (x_unique - r1)) - 2 * exp(-a1 * (x_unique - r1))) + V2 * (exp(-2 * a2 * (r2 - x_unique)) - 2 * exp(-a2 * (r2 - x_unique))) + 0.166 + 0.00019
    H = [(p^2)/(2m) + V for p in p_unique, V in Vx]
    
    Emin = 0
    Emax = 0 + 0.1  
    levels = range(Emin, Emax, length=25)
    
    fig = Figure(size=(1000, 400))
    
    ax1 = Axis(fig[1,1], 
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"U^{\mathrm{G-C}}(x) \; \text{(a.u.)}",
               title = L"\text{Podwójny potencjał Morse'a}",
               xlabelsize = 27,
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    lines!(ax1, x_unique, Vx, linewidth=3, color=:blue)
    text!(ax1, -3.7, 0.024; text = L"\text{(a)}", fontsize = 30)

    ax2 = Axis(fig[1,2], 
               xlabel = L"$x$ \text{ (a.u.)}", 
               ylabel = L"$p$ \text{ (a.u.)}",
               title = L"\text{Hamiltonian w przestrzeni fazowej}",
               xlabelsize = 27,limits = ((-3.8, 2.4), (-15.5, 15.5)),
               ylabelsize = 27, titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    contour!(ax2, x_unique, p_unique, H', levels=levels, linewidth=1.5)
 
    text!(ax2, -3.7, 10.2; text = L"\text{(b)}", fontsize = 30)
   
    # display(fig)
    save("moyal_solver/graphics/hamiltonian_slocombe.pdf", fig)
    return fig
end
# plot_Slocombe_GC_potential()
## exp values

function get_exp_vals()
    data = readdlm("moyal_solver/build/output/stats.dat", skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]
    
    fig = Figure(size = (1000, 500))
    ax1_color = :royalblue1
    ax = Axis(fig[1,1], xlabel = L"t\; \;(10^3\text{a.u.})", ylabel = L"\langle x\rangle\; (\text{a.u.})", 
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
    
    # display(fig)
    save("moyal_solver/graphics/GC/exp_vals.pdf", fig)
end

get_exp_vals()
