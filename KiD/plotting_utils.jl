"""
Plotting utilities
"""

import NCDatasets as NC
import CloudMicrophysics.PrecipitationSusceptibility as CMPS

ENV["GKSwstype"] = "nul"
import ClimaCorePlots, Plots
Plots.GRBackend()
using CairoMakie
using Dates

function plot_profiles_in_time(nc_data_file; output = "output", n = 10)

    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    profs = ds.group["profiles"]

    z_plt = Array(profs["zc"])
    t_plt = Array(profs["t"])
    nt = length(t_plt)

    # Select n evenly distributed time indices
    Δi = round(Int, (nt - 2) / n)
    time_indices = [1; Δi:Δi:(nt-Δi); nt]  # explicitly include first and last indices
    ni = length(time_indices)
    selected_times = t_plt[time_indices]

    kg_to_g = 1e3
    m⁻³_to_cm⁻³ = 1e-6
    q_tot = profs["q_tot"][:, time_indices] .* kg_to_g
    q_liq = profs["q_liq"][:, time_indices] .* kg_to_g
    q_ice = profs["q_ice"][:, time_indices] .* kg_to_g
    q_rai = profs["q_rai"][:, time_indices] .* kg_to_g
    q_sno = profs["q_sno"][:, time_indices] .* kg_to_g
    N_liq = profs["N_liq"][:, time_indices] .* m⁻³_to_cm⁻³
    N_rai = profs["N_rai"][:, time_indices] .* m⁻³_to_cm⁻³
    SN_liq_act = profs["SN_liq_act"][:, time_indices] .* m⁻³_to_cm⁻³

    close(ds)

    # Create colormap for time progression
    # colormap = :viridis
    colormap = cgrad(:darkrainbow, (0:ni) / ni, categorical = true)
    # colors = get(colorschemes[colormap], range(0, 1, length = n))
    ymax = z_plt[end] + (z_plt[end] - z_plt[end-1])  # extend y-axis a bit beyond the last point
    lims(x) = iszero(x) ? (nothing, (0, ymax)) : ((0, 1.1 * maximum(x)), (0, ymax))

    fig = Figure(size = (2000, 1500))

    # Create axes
    ax_q_tot =
        Axis(fig[1, 1]; xlabel = "q_tot [g/kg]", limits = lims(q_tot), ylabel = "z [m]")
    ax_q_liq = Axis(fig[1, 2]; xlabel = "q_liq [g/kg]", limits = lims(q_liq))
    ax_N_liq = Axis(fig[1, 3]; xlabel = "N_liq [1/cm³]", limits = lims(N_liq))
    ax_q_rai =
        Axis(fig[2, 1]; xlabel = "q_rai [g/kg]", limits = lims(q_rai), ylabel = "z [m]")
    ax_N_rai = Axis(fig[2, 2]; xlabel = "N_rai [1/cm³]", limits = lims(N_rai))
    ax_q_ice = Axis(fig[2, 3]; xlabel = "q_ice [g/kg]", limits = lims(q_ice))
    ax_q_sno =
        Axis(fig[3, 1]; xlabel = "q_sno [g/kg]", limits = lims(q_sno), ylabel = "z [m]")
    ax_SN_liq_act =
        Axis(fig[3, 2]; xlabel = "SN_liq_act [1/cm³]", limits = lims(SN_liq_act))

    # Link y-axes
    axs = [ax_q_tot, ax_q_liq, ax_N_liq, ax_q_rai, ax_N_rai, ax_q_ice, ax_q_sno]
    linkyaxes!(axs...)

    # Plot lines for each selected time
    args = (; colormap, colorrange = (0, 1))
    for (i, t) in enumerate(selected_times)
        color = (i - 0.5) / ni  # Normalize color for colormap
        lines!(ax_q_tot, q_tot[:, i], z_plt; args..., color)
        lines!(ax_q_liq, q_liq[:, i], z_plt; args..., color)
        lines!(ax_N_liq, N_liq[:, i], z_plt; args..., color)
        lines!(ax_q_rai, q_rai[:, i], z_plt; args..., color)
        lines!(ax_N_rai, N_rai[:, i], z_plt; args..., color)
        lines!(ax_q_ice, q_ice[:, i], z_plt; args..., color)
        lines!(ax_q_sno, q_sno[:, i], z_plt; args..., color)
        lines!(ax_SN_liq_act, SN_liq_act[:, i], z_plt; args..., color)
    end

    # Add colorbar
    tloc = ((1:ni) .- 0.5) / ni  # locations offset by 0.5 to center label on the color
    tlab = begin
        # Canonicalize time labels
        can_times = @. canonicalize(Second(selected_times))
        # shorten time labels
        map(can_times) do t
            s = string(t)
            s == "empty period" && return "0 s"
            s = replace(
                s,
                "seconds" => "s",
                "minutes" => "m",
                "hours" => "h",
                "days" => "d",
            )
            s = replace(s, "second" => "s", "minute" => "m", "hour" => "h", "day" => "d")
            s = replace(s, "milli" => "m", "micro" => "μ", "nano" => "n")
        end
    end
    Colorbar(fig[:, 4]; colormap, colorrange = (0, 1), width = 20, ticks = (tloc, tlab))

    # Add title
    Label(fig[0, :], "Vertical Profiles at Multiple Time Steps", fontsize = 24)

    # Save figure
    save(joinpath(path, "profiles_multitime.png"), fig)

    return fig
end

function plot_timeheight(nc_data_file; output = "output", pysdm = false)
    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    if pysdm
        t_plt = Array(ds["time"])
        z_plt = Array(ds["height"])
        q_liq_plt = transpose(Array(ds["cloud water mixing ratio"])) / 1e3
        q_rai_plt = transpose(Array(ds["rain water mixing ratio"])) / 1e3
        q_vap = transpose(Array(ds["water_vapour_mixing_ratio"]))
        q_tot_plt = q_vap + q_liq_plt + q_rai_plt
    else
        t_plt = Array(ds.group["profiles"]["t"])
        z_plt = Array(ds.group["profiles"]["zc"])
        q_tot_plt = Array(ds.group["profiles"]["q_tot"])
        q_liq_plt = Array(ds.group["profiles"]["q_liq"])
        q_rai_plt = Array(ds.group["profiles"]["q_rai"])
    end

    p1 = Plots.heatmap(
        t_plt,
        z_plt,
        q_tot_plt .* 1e3,
        title = "q_tot [g/kg]",
        xlabel = "time [s]",
        ylabel = "z [m]",
        color = :BuPu,
        clims = (0, 1),
    )
    p2 = Plots.heatmap(
        t_plt,
        z_plt,
        q_liq_plt .* 1e3,
        title = "q_liq [g/kg]",
        xlabel = "time [s]",
        ylabel = "z [m]",
        color = :BuPu,
        clims = (0, 1),
    )
    p3 = Plots.heatmap(
        t_plt,
        z_plt,
        q_rai_plt .* 1e3,
        title = "q_rai [g/kg]",
        xlabel = "time [s]",
        ylabel = "z [m]",
        color = :BuPu,
        clims = (0, 0.25),
    )
    p = Plots.plot(
        p1,
        p2,
        p3,
        size = (1200.0, 300.0),
        bottom_margin = 30.0 * Plots.PlotMeasures.px,
        left_margin = 30.0 * Plots.PlotMeasures.px,
        layout = (1, 3),
    )
    Plots.png(p, joinpath(path, "timeheight.png"))
end

function plot_cwp_rwp_rr(nc_data_files; output = "output", is_the_last_pysdm = false)
    path = joinpath(@__DIR__, output)
    mkpath(path)

    nc_data_files isa AbstractVector ||
        error("nc_data_files must be an array of NetCDF file names.")

    # Return profile arrays in shape (nz, nt)
    function orient_profile(A, nz, nt)
        A = Array(A)
        if size(A) == (nz, nt)
            return A
        elseif size(A) == (nt, nz)
            return permutedims(A)
        else
            error(
                "Unexpected profile shape $(size(A)); expected ($(nz), $(nt)) or ($(nt), $(nz)).",
            )
        end
    end

    # Reconstruct layer thicknesses from cell centers
    # (used for PySDM if interface heights are not available)
    function thickness_from_centers(zc)
        nz = length(zc)
        Δz = similar(zc)
        if nz == 1
            Δz[1] = one(eltype(zc))
            return Δz
        end
        Δz[1] = zc[2] - zc[1]
        for k = 2:(nz-1)
            Δz[k] = (zc[k+1] - zc[k-1]) / 2
        end
        Δz[nz] = zc[nz] - zc[nz-1]
        return Δz
    end

    function read_one_file(nc_data_file; pysdm = false)
        ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))

        if pysdm
            t = Array(ds["time"])
            zc = Array(ds["height"])

            nz = length(zc)
            nt = length(t)

            q_liq = orient_profile(ds["cloud water mixing ratio"], nz, nt) / 1e3
            q_rai = orient_profile(ds["rain water mixing ratio"], nz, nt) / 1e3
            q_vap = orient_profile(ds["water_vapour_mixing_ratio"], nz, nt)
            ρ = orient_profile(ds["rhod"], nz, nt)
            ρ .*= (1 .+ q_liq .+ q_rai .+ q_vap)

            q_tot = q_vap + q_liq + q_rai

            # PySDM file does not provide zf here, so approximate Δz from centers
            Δz = thickness_from_centers(zc)
        else
            prof = ds.group["profiles"]

            t = Array(prof["t"])
            zc = Array(prof["zc"])
            zf = Array(prof["zf"])

            nz = length(zc)
            nt = length(t)

            q_tot = orient_profile(prof["q_tot"], nz, nt)
            q_liq = orient_profile(prof["q_liq"], nz, nt)
            q_rai = orient_profile(prof["q_rai"], nz, nt)
            ρ = orient_profile(prof["density"], nz, nt)

            length(zf) == nz + 1 || error(
                "Expected length(zf) = nz + 1, but got length(zf) = $(length(zf)) and nz = $(nz).",
            )

            Δz = zf[2:end] .- zf[1:end-1]
        end

        NC.close(ds)

        cwp = vec(sum(ρ .* q_liq .* Δz, dims = 1))
        rwp = vec(sum(ρ .* q_rai .* Δz, dims = 1))
        twp = vec(sum(ρ .* q_tot .* Δz, dims = 1))
        rr = max.(twp[1] .- twp, zero(eltype(twp))) # accumulated precipitation
        t_min = t ./ 60

        return (; t_min, cwp, rwp, rr)
    end

    results = Vector{NamedTuple}(undef, length(nc_data_files))
    labels = [replace(splitext(basename(f))[1], "_" => " ") for f in nc_data_files]

    for i in eachindex(nc_data_files)
        pysdm = is_the_last_pysdm && (i == lastindex(nc_data_files))
        results[i] = read_one_file(nc_data_files[i]; pysdm = pysdm)
        if pysdm
            labels[i] = "PySDM"
        end
    end

    fig = CairoMakie.Figure(size = (1200, 350))

    ax1 = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "time [min]",
        ylabel = "CWP [kg m⁻²]",
        title = "cloud water path",
    )

    ax2 = CairoMakie.Axis(
        fig[1, 2],
        xlabel = "time [min]",
        ylabel = "RWP [kg m⁻²]",
        title = "rain water path",
    )

    ax3 = CairoMakie.Axis(
        fig[1, 3],
        xlabel = "time [min]",
        ylabel = "surface precipitation [kg m⁻²]",
        title = "accumulated surface precipitation",
    )

    for i in eachindex(results)
        r = results[i]
        CairoMakie.lines!(ax1, r.t_min, r.cwp, label = labels[i])
        CairoMakie.lines!(ax2, r.t_min, r.rwp, label = labels[i])
        CairoMakie.lines!(ax3, r.t_min, r.rr, label = labels[i])
    end

    CairoMakie.axislegend(ax1, position = :rc, backgroundcolor = (:white, 0.5))

    CairoMakie.save(joinpath(path, "cwp_rwp_rr.png"), fig)
end
