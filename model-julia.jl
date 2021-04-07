
# let's learn julia

# ==== using Ipopt and JuMP ====

using JuMP
using Ipopt
using Plots
using Roots

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)
@variable(model, 0 <= N₁)
@variable(model, 0 <= N₂)
@variable(model, 0 <= e)
@NLparameter(model, h == 1)
@NLparameter(model, a == 0.1)
@NLparameter(model, b == 0.1)
@NLparameter(model, σ == 0.8)

@NLobjective(model, Max, 
    (1 - e - b * N₁)^(1 - σ)/(1 - σ) +
    ((e * h) * (1 - b * N₂))^(1 - σ)/(1 - σ) + 
    a * (N₁ + N₂)
)

function ipopt_solve(H::Float64, A::Float64, B::Float64, S::Float64)

    set_value(h, H)
    set_value(a, A)
    set_value(b, B)
    set_value(σ, S)

    optimize!(model)

    [value(N₁), value(N₂), value(e)]
end

"solve for e when N₁ == 0"
function solve_e(h, a, b, σ)
    
    rhs = (a / (b*h)) ^ ((1 - σ)/σ^2)
    f(e) = (1 - e) * e ^ ((1 - 2σ)/σ^2) - rhs

    es = find_zeros(f, 0.001, 1)
    if size(es) != (1,)
        @warn "Found more or less than 1 value for e"
        @warn "h $h; a $a; b $b; σ $σ"
    end

    es[1]
end

function theory_solve(h::Float64, a::Float64, b::Float64, σ::Float64)

    N₁ = 1/b * (1 - (b/a) ^ (1/(2σ - 1)) * h ^ ((1 - σ)/(2σ - 1))) - 
         (1/b) * (b/a) ^ (1/σ)

    N₁ = max(0, N₁)

    e = if N₁ == 0
        solve_e(h, a, b, σ)
    else
        (b/a) ^ (1/(2σ-1)) * h ^ ((1 - σ)/(2σ - 1))
    end

    N₂ = 1/b * (1 - (b/a) ^ (1/σ) * (e * h) ^ ((1 - σ)/σ))
    N₂ = max(0, N₂)

    if N₂ == 0
        e = 1/(1 + h ^ ((σ-1)/σ))
    end

    [N₁, N₂, e]
end

function plot_solutions(hmin, hmax, a, b, σ; theory = true)
    len = 100
    hs = range(hmin, hmax, length = len)

    ipopt_result = Array{Float64}(undef, 3, len)
    theory_result = copy(ipopt_result)

    for (ix, h) = zip(1:len, hs) 
        ipopt_result[:, ix]  = ipopt_solve(h, a, b, σ)
        theory_result[:, ix] = theory_solve(h, a, b, σ)
    end

    p = plot(hs, ipopt_result[1, :], label = "N1")
    plot!(p, hs, ipopt_result[2, :], label = "N2")
    plot!(p, hs, ipopt_result[3, :], label = "e")
    if theory
        plot!(p, hs, theory_result[1, :], label = "N1 theory", 
            linestyle = :dash, lw = 2)
        plot!(p, hs, theory_result[2, :], label = "N2 theory", 
            linestyle = :dash, lw = 2)
        plot!(p, hs, theory_result[3, :], label = "e theory",  
            linestyle = :dash, lw = 2)
    end

    display(p)
end