

using JuMP, Ipopt, Plots

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)

@variable(model, 0 <= N₁)
@variable(model, 0 <= N₂)
@variable(model, 0 <= s <= 1)

@NLparameter(model, h == 1)
@NLparameter(model, a == 0.1)
@NLparameter(model, b == 0.1)
@NLparameter(model, σ == 0.8)

@NLconstraint(model, 1 - b * N₂ >= 0)
@NLconstraint(model, 1 - s - b * N₁ >= 0)

@NLobjective(model, Max,
          - exp(-σ*(1-s-b*N₁)) 
          - exp(-σ * (s*h)*(1-b*N₂)) 
          + a*(N₁ + N₂)
       )

hhs = 0.1:0.1:2
ipopt_result = Array{Float64}(undef, 3, length(hhs))

for (ix, hh) in enumerate(hhs) 
    set_value(h, hh)
    optimize!(model)
    ipopt_result[:, ix] = [value(N₁), value(N₂), value(s)]
end


plot(hhs, ipopt_result[1,:] + ipopt_result[2,:], label = "N")
plot!(hhs, ipopt_result[1,:], label = "N1", linestyle = :dash)
plot!(hhs, ipopt_result[2,:], label = "N2", linestyle = :dash)
plot!(hhs, ipopt_result[3,:], label = "s", color = "black", 
    linestyle = :dot)