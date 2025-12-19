"Box model with variable number of modes, fixed thresholds, & hydrodynamic kernel"

using OrdinaryDiffEq

include("../utils/box_model_helpers.jl")
include("../utils/netcdf_helpers.jl")

FT = Float64

"""
Initializes a set of gamma distributions.
    - `init_moments` - initial moments for the first mode
    - `num_modes` - how many gamma modes to use 
"""
function init_conditions_gamma(init_moments, num_modes)
    if num_modes > 1
        thresholds = 10 .^ range(-9, -3; length=num_modes)
        thresholds[end] = Inf
        thresholds = Tuple(thresholds)
    else
        thresholds = (Inf,)
    end
    
    mom_init = zeros(FT, num_modes*3)
    mom_init[1:3] .= init_moments[1:3]

    # initialization of dists is arbitrary, as it gets updated in the time stepping
    dist_init = Vector{GammaPrimitiveParticleDistribution}(undef, num_modes)
    for i in 1:num_modes
        dist_init[i] = GammaPrimitiveParticleDistribution(FT(0), FT(1), FT(1))
        dist_init[i] = update_dist_from_moments(dist_init[i], Tuple(mom_init[3*(i-1)+1:3*i]))
    end
    dist_init = Tuple(dist_init)

    NProgMoms = map(dist_init) do dist
        nparams(dist)
    end
    return (dist_init, NProgMoms, mom_init, thresholds)
end


# Common parameters
moments_init = [1e8, 1e-2, 2e-12] # single gamma initial condition
kernel_func = HydrodynamicKernelFunction(1e2 * π) # 1e2 π m^3/kg^(4/3)/s
kernel = CoalescenceTensor(kernel_func, 4, FT(1e-6))
tspan = (FT(0), FT(240))
norms = (1e6, 1e-9) # 1e6/m^3; 1e-9 kg

n_dist_list = (1, 2, 4)
for nd in n_dist_list
    print("computing for $(nd) distributions \n")
    (dist_init, NProgMoms, mom_init, thresholds) = init_conditions_gamma(moments_init, nd)
    coal_data = CoalescenceData(kernel, NProgMoms, thresholds, norms)
    rhs = make_box_model_rhs(AnalyticalCoalStyle())
    ODE_parameters = (; pdists = dist_init, coal_data = coal_data, NProgMoms = NProgMoms, norms = norms, dt = FT(10))

    prob = ODEProblem(rhs, mom_init, tspan, ODE_parameters)
    sol = solve(prob, SSPRK33(), dt = ODE_parameters.dt)

    box_output(sol, ODE_parameters, "results/fig5_box_$(nd)_gamma_hydro.nc", FT)
end