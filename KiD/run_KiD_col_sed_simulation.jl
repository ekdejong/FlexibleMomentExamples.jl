import OrdinaryDiffEq as ODE
import ClimaCore as CC
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import KinematicDriver
import KinematicDriver.Common as CO
import KinematicDriver.K1DModel as K1D

include("create_parameters.jl")
include("plotting_utils.jl")

function run_KiD_col_sed_simulation(::Type{FT}, opts) where {FT}

    # Equations to solve for precipitation variables
    precipitation_choice = opts["precipitation_choice"]
    rain_formation_choice = opts["rain_formation_choice"]
    sedimentation_choice = opts["sedimentation_choice"]
    @info precipitation_choice, rain_formation_choice, sedimentation_choice

    # Decide the output flder name based on options
    output_folder = opts["output_folder"]
    path = joinpath(opts["root_path"], output_folder)
    mkpath(path)

    # Overwrite the defaults parameters based on options
    default_toml_dict = CP.create_toml_dict(FT)
    toml_dict = override_toml_dict(
        default_toml_dict;
        precip_sources = 1,
        precip_sinks = 0,
        prescribed_Nd = FT(opts["prescribed_Nd"]),
    )
    toml_dict.data["SB2006_cloud_gamma_distribution_coeff_nu"]["value"] = opts["k"]
    # Create Thermodynamics.jl and KinematicDriver model parameters
    # (some of the CloudMicrophysics.jl parameters structs are created later based on model choices)
    common_params = create_common_parameters(toml_dict)
    kid_params = create_kid_parameters(toml_dict)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    air_params = CMP.AirProperties(toml_dict)
    activation_params = CMP.AerosolActivationParameters(toml_dict)

    moisture_choice =
        precipitation_choice == "CloudyPrecip" ? "CloudyMoisture" : "NonEquilibriumMoisture"
    moisture = CO.get_moisture_type(moisture_choice, toml_dict)
    precip = CO.get_precipitation_type(
        precipitation_choice,
        toml_dict,
        rain_formation_choice = rain_formation_choice,
        sedimentation_choice = sedimentation_choice,
    )

    @info "Initialising"
    # Initialize the timestepping struct
    TS = CO.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

    # Create the coordinates
    space, face_space =
        K1D.make_function_space(FT, FT(opts["z_min"]), FT(opts["z_max"]), opts["n_elem"])
    coord = CC.Fields.coordinate_field(space)
    face_coord = CC.Fields.coordinate_field(face_space)

    # Initialize the netcdf output Stats struct
    fname = joinpath(path, opts["output_nc_file"])
    Stats = CO.NetCDFIO_Stats(
        fname,
        1.0,
        parent(face_coord),
        parent(coord),
        output_profiles = Dict(
            :ρ => "density",
            :q_tot => "q_tot",
            :q_liq => "q_liq",
            :q_rai => "q_rai",
            :N_liq => "N_liq",
            :N_rai => "N_rai",
        ),
    )

    # Create the initial condition profiles
    ic_0d = initial_condition_0d(
        FT,
        thermo_params,
        opts["qt"],
        opts["prescribed_Nd"],
        opts["k"],
        opts["rhod"],
    )
    ic_zero = initial_condition_0d(
        FT,
        thermo_params,
        FT(0),
        eps(FT),
        opts["k"],
        opts["rhod"],
    )
    z_bot = FT(opts["z_bot"])
    z_top = FT(opts["z_top"])
    if precipitation_choice == "CloudyPrecip"
        cloudy_disttypes = determine_cloudy_disttypes(opts["num_moments"], false)
        cloudy_params, cloudy_pdists =
            create_cloudy_parameters(FT, cloudy_disttypes, opts["kernel"])
        ic_cloudy = CO.cloudy_initial_condition(cloudy_pdists, ic_0d, opts["k"])
        ic_cloudy_zero = CO.cloudy_initial_condition(cloudy_pdists, ic_zero, opts["k"])
        init = map(coord) do c
            c.z >= z_bot && c.z <= z_top ? ic_cloudy : ic_cloudy_zero
        end
    else
        cloudy_params = nothing
        init = map(coord) do c
            c.z >= z_bot && c.z <= z_top ? ic_0d : ic_zero
        end
    end

    # Create aux vector and apply initial condition
    aux = K1D.initialise_aux(
        FT,
        init,
        common_params,
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        TS,
        Stats,
        face_space,
        moisture,
        precip,
#        "ShipwayHill2012",
        cloudy_params,
    )

    # Create state vector and apply initial condition
    Y = CO.initialise_state(moisture, precip, init)

    # Output the initial condition
    CO.simulation_output(aux, 0.0)

    # Define callbacks for output
    callback_io = ODE.DiscreteCallback(
        CO.condition_io,
        CO.affect_io!;
        save_positions = (false, false),
    )
    callbacks = ODE.CallbackSet(callback_io)

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = K1D.make_rhs_function_col_sed(moisture, precip)

    # Solve the ODE operator
    problem = ODE.ODEProblem(ode_rhs!, Y, (FT(opts["t_ini"]), FT(opts["t_end"])), aux)
    @info "Solving"
    solver = ODE.solve(
        problem,
        ODE.SSPRK33(),
        dt = TS.dt,
        saveat = TS.dt_io,
        callback = callbacks,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )

    # Some basic plots
    @info "Plotting"
    plot_timeheight(fname, output = path)
end

opts_common = Dict(
    "qt" => 1e-3,
    "prescribed_Nd" => 1e7,
    "k" => 2.0,
    "rhod" => 1.0,
    "precipitation_choice" => nothing,
    "rain_formation_choice" => nothing,
    "sedimentation_choice" => nothing,
    "num_moments" => 6,
    "z_min" => 0.0,
    "z_max" => 3000.0,
    "z_bot" => 1500.0,
    "z_top" => 2250.0,
    "n_elem" => 60,
    "dt" => 1.0,
    "dt_output" => 10.0,
    "t_ini" => 0.0,
    "t_end" => 1500.0,
    "root_path" => joinpath(@__DIR__, "Output_KiD_col_sed"),
    "output_folder" => "output",
    "output_nc_file" => "Output.nc",
)
opts_cloudy4 = merge(
    opts_common,
    Dict(
        "precipitation_choice" => "CloudyPrecip",
        "num_moments" => 4,
        "kernel" => "Long",
        "output_folder" => "CloudyPrecip_4",
        "output_nc_file" => "Flexible_4M,_Long.nc",
    ),
)
opts_cloudy6 = merge(
    opts_common,
    Dict(
        "precipitation_choice" => "CloudyPrecip",
        "num_moments" => 6,
        "kernel" => "Long",
        "output_folder" => "CloudyPrecip_6",
        "output_nc_file" => "Flexible_6M,_Long.nc",
    ),
)
opts_cloudy7 = merge(
    opts_common,
    Dict(
        "precipitation_choice" => "CloudyPrecip",
        "num_moments" => 7,
        "kernel" => "Long",
        "output_folder" => "CloudyPrecip_7",
        "output_nc_file" => "Flexible_7M,_Long.nc",
    ),
)
opts_cloudy6_gol = merge(
    opts_common,
    Dict(
        "precipitation_choice" => "CloudyPrecip",
        "num_moments" => 6,
        "kernel" => "Golovin",
        "output_folder" => "CloudyPrecip_6_gol",
        "output_nc_file" => "Flexible_6M,_Golovin.nc",
    ),
)
opts_1m = merge(
    opts_common,
    Dict(
        "precipitation_choice" => "Precipitation1M",
        "rain_formation_choice" => "CliMA_1M",
        "sedimentation_choice" => "CliMA_1M",
        "output_folder" => "Precipitation1M",
        "output_nc_file" => "1M,_Kessler.nc",
    ),
)
opts_2m = merge(
    opts_common,
    Dict(
        "precipitation_choice" => "Precipitation2M",
        "rain_formation_choice" => "SB2006",
        "sedimentation_choice" => "SB2006",
        "output_folder" => "Precipitation2M",
        "output_nc_file" => "2M,_SB2006.nc",
    ),
)

multiple_models_opts = (opts_cloudy6,) #(opts_cloudy4, opts_cloudy6, opts_cloudy7, opts_1m, opts_2m)
for opts in multiple_models_opts
    run_KiD_col_sed_simulation(Float64, opts)
end

# plot cwp, rwp and rr 
output_nc_path(_opts) =
    joinpath(_opts["root_path"], _opts["output_folder"], _opts["output_nc_file"])
data_files = [
    output_nc_path(_opts) for
    _opts in multiple_models_opts
]
pysdm_file = joinpath(
    opts_common["root_path"],
    "../../results/pysdm/pysdm_colSed_partDomain_N0=10.nc",
)
plot_timeheight(pysdm_file, output = opts_common["root_path"], pysdm = true)
plot_cwp_rwp_rr(
    [data_files..., pysdm_file],
    output = opts_common["root_path"],
    is_the_last_pysdm = true,
)
