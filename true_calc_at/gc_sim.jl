using DelimitedFiles
include("utils.jl")

using CairoMakie

function double_morse(q; V1=0.1617, V2=0.082, a1=0.305, a2=0.755, r1=-2.7, r2=2.1)
    const_term = 0.166 + 0.00019
    exp1 = exp.(-2 .* a1 .* (q .- r1))
    exp2 = exp.(-a1 .* (q .- r1))
    exp3 = exp.(-2 .* a2 .* (r2 .- q))
    exp4 = exp.(-a2 .* (r2 .- q))
    return V1 .* (exp1 .- 2 .* exp2) .+ V2 .* (exp3 .- 2 .* exp4) .+ const_term
end

function simulate_double_morse(; T=100000, dt=1000)
    times = 0:dt:T
    r2_series = [2.1 + 0.3*sin(0.0001*t) for t in times] 
    return times, r2_series
end

function make_gif_double_morse(; filename="/home/marta/Documents/studia/masters/graphics/true_sim/double_morse.gif")
    times, r2_series = simulate_double_morse()
    xs = LinRange(-4, 3, 400)
    fig = Figure(resolution=(800,600))
    ax = Axis(fig[1,1], xlabel=L"$x$", ylabel=L"$U(x)$", title=L"\text{Double Morse G-C}", 
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    lines!(ax, xs, double_morse(xs; r2=2.1), color=:blue)
    CairoMakie.ylims!(ax, -0.002, 0.033)
    ys = Observable(double_morse(xs; r2=2.1))
    lineplot = lines!(ax, xs, ys, color=:red, label="ewolucja w czasie")
    axislegend(ax, position = :rb)
    record(fig, filename, 1:length(r2_series); framerate=20) do i
        ys[] = double_morse(xs; r2=r2_series[i])
    end
end

make_gif_double_morse()
