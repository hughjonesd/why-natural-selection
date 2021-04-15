

using JuMP, Ipopt, Plots

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)

@variable(model, 0 <= N₁)
@variable(model, 0 <= N₂)
@variable(model, 0 <= s <= 1)

@NLparameter(model, h == 1)
@NLparameter(model, a == 0.1)
@NLparameter(model, b == 0.1)
@NLparameter(model, γ == 0.8)

@NLconstraint(model, 1 - b * N₂ >= 0)
@NLconstraint(model, 1 - s - b * N₁ >= 0)

@NLobjective(model, Max,
          - exp(-γ * (1 - s - b * N₁)) 
          - exp(-γ * (s * h) * (1 - b * N₂)) 
          + a * (N₁ + N₂)
       )


hhs = 0.1:0.02:6
ipopt_result = Array{Float64}(undef, 3, length(hhs))


set_value(a, 0.2)
set_value(b, 0.8)
set_value(γ, 1.5)

for (ix, hh) in enumerate(hhs) 
    set_value(h, hh)
    optimize!(model)
    ipopt_result[:, ix] = [value(N₁), value(N₂), value(s)]
end


plot(hhs, ipopt_result[1,:] + ipopt_result[2,:], label = "N");
plot!(hhs, ipopt_result[1,:], label = "N1", linestyle = :dash);
plot!(hhs, ipopt_result[2,:], label = "N2", linestyle = :dash);
plot!(hhs, ipopt_result[3,:], label = "s", color = "black", 
    linestyle = :dot)


# not CARA but with a money cost

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)
@variable(model, 0 <= N₁)
@variable(model, 0 <= N₂)
@variable(model, 0 <= s <= 1)
@NLparameter(model, h == 1)
@NLparameter(model, a == 0.1)
@NLparameter(model, b == 0.1)
@NLparameter(model, m == 0.1)
@NLparameter(model, σ == 0.8)

@NLconstraint(model, 1 - s - b * N₁ - m * N₁ >= 0)
@NLconstraint(model, (s * h) * (1 - b * N₂) - m * N₂ >= 0)

@NLobjective(model, Max, 
    (1 - s - b * N₁ - m * N₁)^(1 - σ)/(1 - σ) +
    ((s * h) * (1 - b * N₂) - m * N₂)^(1 - σ)/(1 - σ) + 
    a * (N₁ + N₂)
)

hhs = 0.1:0.01:0.8
ipopt_result = Array{Float64}(undef, 3, length(hhs))


set_value(a, 0.4)
set_value(b, 0.175)
set_value(σ, 0.7)
set_value(m, 0.075)

for (ix, hh) in enumerate(hhs) 
    set_value(h, hh)
    optimize!(model)
    ipopt_result[:, ix] = [value(N₁), value(N₂), value(s)]
end


plot(hhs, ipopt_result[1,:] + ipopt_result[2,:] .+ 0.03, label = "N", ylims = (0,2), color = "black", xlabel = "h");
plot!(hhs, ipopt_result[1,:], label = "N1", linestyle = :dash,
    color = "black");
plot!(hhs, ipopt_result[2,:], label = "N2", linestyle = :dashdot,
    color = "black")
 plot!(hhs, ipopt_result[3,:], label = "s", color = "black", 
     linestyle = :dot)
