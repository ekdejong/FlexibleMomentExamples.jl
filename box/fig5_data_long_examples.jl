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
function init_conditions_gamma(init_moments, num_modes; moving=false)
    if num_modes > 1
        if moving
            thresholds = 0.99 * ones(num_modes)
            thresholds[end] = 1.0
        else
            thresholds = 5 * 10 .^ range(-10, -6; length=num_modes)
            thresholds[end] = Inf
        end
        thresholds = Tuple(thresholds)
    else
        if moving
            thresholds = (1.0,)
        else
            thresholds = (Inf,)
        end
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

"""
Initializes a set of exponential distributions.
    - `init_moments` - initial moments for the first mode
    - `num_modes` - how many exp modes to use 
"""
function init_conditions_exp(init_moments, num_modes)
    if num_modes > 1
        thresholds = 10 .^ range(-10, -5; length=num_modes)
        thresholds[end] = Inf
        thresholds = Tuple(thresholds)
    else
        thresholds = (Inf,)
    end
    
    mom_init = zeros(FT, num_modes*2)
    mom_init[1:2] .= init_moments[1:2]

    # initialization of dists is arbitrary, as it gets updated in the time stepping
    dist_init = Vector{ExponentialPrimitiveParticleDistribution}(undef, num_modes)
    for i in 1:num_modes
        dist_init[i] = ExponentialPrimitiveParticleDistribution(FT(0), FT(1))
        dist_init[i] = update_dist_from_moments(dist_init[i], Tuple(mom_init[2*(i-1)+1:2*i]))
    end
    dist_init = Tuple(dist_init)

    NProgMoms = map(dist_init) do dist
        nparams(dist)
    end
    return (dist_init, NProgMoms, mom_init, thresholds)
end


# Common parameters
moments_init = [1e8, 1e-2, 2e-12] # single gamma initial condition
kernel_func = LongKernelFunction(5e-10, 9.44e9, 5.78) # 5.236e-10 kg; 9.44e9 m^3/kg^2/s; 5.78 m^3/kg/s
tspan = (FT(0), FT(120))
norms = (1e6, 1e-9) # 1e6/m^3; 1e-9 kg

n_dist_list = (2, 4, 8)
for nd in n_dist_list
    print("computing for $(nd) distributions \n")
    (dist_init, NProgMoms, mom_init, thresholds) = init_conditions_gamma(moments_init, nd)
    matrix_of_kernels = ntuple(nd) do i
        ntuple(nd) do j
            if thresholds[i] < 5e-10 && thresholds[j] < 5e-10
                CoalescenceTensor(kernel_func, 2, FT(5e-10))
            else
                CoalescenceTensor(kernel_func, 2, FT(1e-6), FT(5e-10))
            end
        end
    end
    coal_data = CoalescenceData(matrix_of_kernels, NProgMoms, thresholds, norms)
    rhs = make_box_model_rhs(AnalyticalCoalStyle())
    ODE_parameters = (; pdists = dist_init, coal_data = coal_data, NProgMoms = NProgMoms, norms = norms, dt = FT(1))

    prob = ODEProblem(rhs, mom_init, tspan, ODE_parameters)
    sol = solve(prob, SSPRK33(), dt = ODE_parameters.dt)

    # plot_params!(sol, ODE_parameters; file_name = "figures/fig5/box_$(nd)_gamma_long_params.pdf", yscale=:identity)
    # plot_moments!(sol, ODE_parameters; file_name = "figures/fig5/box_$(nd)_gamma_long_moments.pdf")
    # plot_spectra!(sol, ODE_parameters; file_name = "figures/fig5/box_$(nd)_gamma_long_spectra.pdf", logxrange = (-12, -3))
    box_output(sol, ODE_parameters, "results/fig5_box_$(nd)_gamma_long.nc", FT)
end