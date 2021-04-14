

module Fertility

export U, ipopt_solve, theory_solve, solve_s, plot_solutions

using JuMP
using Ipopt
using Plots
using Roots

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)
@variable(model, 0 <= N₁)
@variable(model, 0 <= N₂)
@variable(model, 0 <= s <= 1)
@NLparameter(model, h == 1)
@NLparameter(model, a == 0.1)
@NLparameter(model, b == 0.1)
@NLparameter(model, σ == 0.8)

@NLobjective(model, Max, 
    (1 - s - b * N₁)^(1 - σ)/(1 - σ) +
    ((s * h) * (1 - b * N₂))^(1 - σ)/(1 - σ) + 
    a * (N₁ + N₂)
)


function U(N₁::Real, N₂::Real, s::Real, h::Real, a::Real, b::Real, σ::Real)
    (1 - s - b * N₁)^(1 - σ)/(1 - σ) +
    ((s * h) * (1 - b * N₂))^(1 - σ)/(1 - σ) + 
    a * (N₁ + N₂)
end


function ipopt_solve(H::Real, A::Real, B::Real, S::Real)

    set_value(h, H)
    set_value(a, A)
    set_value(b, B)
    set_value(σ, S)

    optimize!(model::JuMP.Model)

    [value(N₁), value(N₂), value(s)]
end

"solve for s when N₁ == 0"
function solve_s(h, a, b, σ)
    
    rhs = (a / (b*h)) ^ ((1 - σ)/σ^2)
    f(e) = (1 - e) * e ^ ((1 - 2σ)/σ^2) - rhs

    ss = find_zeros(f, 0.001, 1)
    if length(ss) != 1
        @warn "Found $(length(ss)) values for e"
        @warn "h $h; a $a; b $b; σ $σ"
    end

    ss[1]
end

function theory_solve(h::Real, a::Real, b::Real, σ::Real)

    N₁ = 1/b * (1 - (b/a) ^ (1/(2σ - 1)) * h ^ ((1 - σ)/(2σ - 1))) - 
         (1/b) * (b/a) ^ (1/σ)

    N₁ = max(0, N₁)

    s = if N₁ == 0
        solve_s(h, a, b, σ)
    else
        (b/a) ^ (1/(2σ-1)) * h ^ ((1 - σ)/(2σ - 1))
    end

    N₂ = 1/b * (1 - (b/a) ^ (1/σ) * (s * h) ^ ((1 - σ)/σ))
    N₂ = max(0, N₂)

    if N₂ == 0
        s = 1/(1 + h ^ ((σ-1)/σ))
    end

    [N₁, N₂, s]
end

function plot_solutions(hmin, hmax, a, b, σ; theory = true, ipopt = true)
    len = 100
    hs = range(hmin, hmax, length = len)

    ipopt_result = Array{Float64}(undef, 3, len)
    theory_result = copy(ipopt_result)

    for (ix, h) = zip(1:len, hs) 
        ipopt && (ipopt_result[:, ix]  = ipopt_solve(h, a, b, σ))
        theory && (theory_result[:, ix] = theory_solve(h, a, b, σ))
    end

    if ipopt
        p = plot(hs, ipopt_result[1, :], label = "N1 ipopt")
        plot!(p, hs, ipopt_result[2, :], label = "N2 ipopt")
        plot!(p, hs, ipopt_result[3, :], label = "s ipopt")
    end
    if theory
        if ! ipopt
            p = plot(hs, theory_result[1, :], label = "N₁", 
                    linestyle = :dash, lw = 2, xlabel = "h", lc = :black)
        else
            plot!(p, hs, theory_result[1, :], label = "N₁", 
                linestyle = :dash, lw = 2, lc = :black)
        end
        plot!(p, hs, theory_result[2, :], label = "N₂", 
            linestyle = :dashdot, lw = 2, lc = :black)
        plot!(p, hs, theory_result[1, :] + theory_result[2, :] .+ 0.1, label = "N", 
            linestyle = :solid, lw = 2, lc = :black)
        # plot!(p, hs, theory_result[3, :], label = "s",  
        #     linestyle = :dash, lw = 2)
    end

    display(p)
end

function dN2dh(h, a, b ,σ)

    s = solve_s(h, a, b, σ)

    N1, N2, _ = theory_solve(h, a, b, σ) 
    N1 == 0 || (@warn "N1 not 0")
    N2 > 0  || (@warn "N2 not positive")

    W = (b/a) ^ ((σ-1)/(σ^2)) * (σ-1)/(σ^2) 
    X = (1-2σ)/σ^2 * s ^ ((1-2σ)/σ^2 - 1) - (1 - σ)^2/σ^2 * s ^ ((1-2σ)/σ^2)
    W /= X

    dsdh = W * h ^ ((σ-1)/(σ^2) - 1)

    dN2dh = - 1/b * (b/a)^(1/σ) * (1 - σ)/σ *  (s * h)^((1 - 2σ)/σ) * (s + h * dsdh)

    dN2dh
end

end