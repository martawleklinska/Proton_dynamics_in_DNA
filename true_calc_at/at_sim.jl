include("utils.jl")

function evolve_at(is_at::Bool = true)
    idx = is_at ? 4 : 10
    a_coefs = hermite_coefficients(idx; is_at)
    T, dt = 10000, 100.
    times = 0:dt:T

    a_series = Float64[]
    b_series = Float64[]
    for t in times
        a_coefs[5] = 0.0207 + 0.005 * sin(0.00095 * t)
        a_coefs[3] = -0.0414 + 0.005 * sin(0.00095 * t + 500)

        push!(a_series, a_coefs[5])
        push!(b_series, a_coefs[3])
    end
    return a_series, b_series
end

function make_gif_at(a_series, b_series; filename="/home/marta/Documents/studia/masters/graphics/true_sim/fourth_at.gif", is_at::Bool=true)
    nframes = length(a_series)
    title = Observable("AT ")

    fig = Figure(size=(600,400))
    ax = Axis(fig[1,1], xlabel=L"$x$", ylabel=L"$U(x)$", title=title, 
    xlabelsize = 25, ylabelsize = 25, titlesize = 20, xticklabelsize = 17, yticklabelsize = 17)
    xs = Observable(LinRange(-3.7, 4.0, 200))
    ys = Observable(zeros(length(xs[])))

    lines!(ax, xs, ys, color=:blue)
    CairoMakie.ylims!(ax, -0.005, 0.045)

    record(fig, filename, 1:nframes; framerate=20) do i
        a4 = a_series[i]
        a2 = b_series[i]
        a3, a1, a0 = -0.0053, 0.0158, 0.0312
        alpha = 1.963
        ys[] = [a4*(x/alpha)^4 + a3*(x/alpha)^3 + a2*(x/alpha)^2 + a1*(x/alpha) + a0 for x in xs[]]
        title[] = "A-T t = $i"
    end
end

a_s, b_s = evolve_at()
make_gif_at(a_s, b_s)
