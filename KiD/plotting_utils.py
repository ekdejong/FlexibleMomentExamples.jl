# Improved plotting code for the notebook

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# ============================================================================
# 1. IMPROVED VERTICAL PROFILE PLOT
# ============================================================================

def plot_vertical_profiles(datasets, t_plt=[0, 300, 900], figsize=(10, 8), lw=None, linestyles=None, ylim=3000):
    """
    Plot vertical profiles of qc, qr, Nc, Nr at multiple times.

    Parameters:
    -----------
    datasets : dict
        Dictionary of xarray datasets with keys as dataset names
    t_plt : list
        List of times (in seconds) to plot
    """
    # Naming conventions
    pv = {
        "qc": "cloud water mixing ratio",
        "qr": "rain water mixing ratio",
        "Nc": "nc",
        "Nr": "nr",
    }
    kv = {
        "qc": "q_liq",
        "qr": "q_rai",
        "Nc": "N_liq",
        "Nr": "N_rai",
    }
    factors = {
        "qc": 1000,  # kg/kg -> g/kg
        "qr": 1000,  # kg/kg -> g/kg
        "Nc": 1e-6,  # m^-3 -> cm^-3
        "Nr": 1e-6,  # m^-3 -> cm^-3
    }
    labels = {
        "qc": "Cloud water mixing ratio (g/kg)",
        "qr": "Rain water mixing ratio (g/kg)",
        "Nc": "Cloud droplet number (cm$^{-3}$)",
        "Nr": "Rain drop number (cm$^{-3}$)",
    }

    # Color and line style configuration
    colors = ['k', "tab:blue", "tab:green", "tab:orange", "tab:red", "tab:purple"]
    if linestyles is None:
        linestyles = ['-.', '-', '--']
    if lw is None:
        lw = [2, 3, 3]

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=figsize, sharey=True)
    axes = axes.flatten()

    # Plot each variable
    for i, (key, value) in enumerate(pv.items()):
        ax = axes[i]
        ax.set_xlabel(labels[key], fontsize=11)
        ax.grid(True, alpha=0.3, linestyle='--')

        for j, t in enumerate(t_plt):
            for k, (dataset_name, dataset) in enumerate(datasets.items()):
                if dataset_name == "PySDM":
                    data = dataset[pv[key]].sel(time=t, method='nearest') * factors[key]
                    z = dataset['height']
                    label = f't={t}s' if (i == 0 and k == 0) else (dataset_name if (j == 0) else None)
                else:
                    data = dataset[kv[key]].sel(t=t, method='nearest') * factors[key]
                    z = dataset['zc']
                    label = dataset_name if (i == 0 and j == 0 and k > 0) else None

                ax.plot(data, z, linestyles[j], color=colors[k],
                       label=label, linewidth=lw[j])
        ax.set_xlim(left=0)
        ax.set_ylim(0, ylim)

    # Set y-label only for left column
    axes[0].set_ylabel('Altitude (m)', fontsize=11)
    axes[2].set_ylabel('Altitude (m)', fontsize=11)

    # Create legends
    # Top-left: time labels
    time_lines = [plt.Line2D([0], [0], color='k', linestyle=linestyles[j],
                             label=f't={t}s') for j, t in enumerate(t_plt)]
    axes[0].legend(handles=time_lines, loc='best', framealpha=0.9, fontsize=9)

    # Top-right: dataset labels
    dataset_lines = [plt.Line2D([0], [0], color=colors[k], linestyle='-',
                                label=name)
                     for k, name in enumerate(datasets.keys())]
    axes[1].legend(handles=dataset_lines, loc='best', framealpha=0.9, fontsize=9)

    plt.tight_layout()
    return fig, axes


# ============================================================================
# 2. TIME SERIES PLOT (CWP/RWP, RAIN RATE, EFFECTIVE RADIUS)
# ============================================================================

def compute_path_integrated_quantities(dataset, is_pysdm=False, use_flux=False):
    """
    Compute column-integrated quantities from a dataset.

    Parameters:
    -----------
    dataset : xarray.Dataset
        Dataset containing the simulation data
    is_pysdm : bool
        Whether this is PySDM data (different coordinate convention)
    use_flux : bool
        If True, compute rain rate as downward flux through bottom of domain.
        If False, compute rain rate as -d(TWP)/dt (only valid when TWP is constant).
        Falls back to -d(TWP)/dt if flux data is not available.

    Returns:
    --------
    dict with keys: t, cwp, rwp, rr (rain rate), Nc_col, Nr_col (column number)
    """
    if is_pysdm:
        t = dataset['time'].values
        z = dataset['height'].values

        # Get mixing ratios (assuming already in kg/kg from notebook preprocessing)
        q_liq = dataset['cloud water mixing ratio'].values  # kg/kg
        q_rai = dataset['rain water mixing ratio'].values  # kg/kg
        q_vap = dataset['water_vapour_mixing_ratio'].values

        # Get number concentrations
        N_liq = dataset['nc'].values  # m^-3
        N_rai = dataset['nr'].values  # m^-3

        # Get density
        rhod = dataset['rhod'].values
        rho = rhod * (1 + q_liq + q_rai + q_vap)

        # Approximate layer thickness from cell centers
        dz = np.diff(z)
        dz = np.concatenate([[dz[0]], (dz[:-1] + dz[1:]) / 2, [dz[-1]]])

    else:
        prof = dataset
        t = prof['t'].values
        zc = prof['zc'].values
        zf = prof['zf'].values

        # Get mixing ratios
        q_liq = prof['q_liq'].values
        q_rai = prof['q_rai'].values
        q_tot = prof['q_tot'].values

        # Get number concentrations
        N_liq = prof['N_liq'].values  # m^-3
        N_rai = prof['N_rai'].values  # m^-3

        # Get density
        rho = prof['density'].values

        # Layer thickness
        dz = np.diff(zf)

    # Compute column water paths (kg/m^2) and domain-averaged number concentrations (cm^-3)
    # Note: data shape is (time, z) for KiD, (z, time) for PySDM
    if is_pysdm:
        # PySDM arrays are (z, time), so we sum along axis=0 (z-axis)
        # dz needs to be broadcast as (z, 1)
        dz_broadcast = dz[:, np.newaxis]
        cwp = np.sum(rho * q_liq * dz_broadcast, axis=0)
        rwp = np.sum(rho * q_rai * dz_broadcast, axis=0)

        # Domain-averaged number concentrations (weighted by dz)
        total_height = np.sum(dz)
        Nc_avg = np.sum(N_liq * dz_broadcast, axis=0) / total_height * 1e-6  # m^-3 -> cm^-3
        Nr_avg = np.sum(N_rai * dz_broadcast, axis=0) / total_height * 1e-6  # m^-3 -> cm^-3

        if use_flux and 'surface precipitation' in dataset:
            # Use PySDM's surface precipitation (units: m/s)
            # Convert to kg/m^2/s by multiplying by water density
            rho_water = 1000.0  # kg/m^3
            rr = dataset['surface precipitation'].values * rho_water  # kg m^-2 s^-1

            # Apply moving average to smooth noise
            window_size = 10
            rr_smoothed = np.convolve(rr, np.ones(window_size)/window_size, mode='same')
            rr_smoothed[:window_size//2] = rr[:window_size//2]
            rr_smoothed[-window_size//2:] = rr[-window_size//2:]
            rr = rr_smoothed
        else:
            # Compute rain rate as -d(TWP)/dt
            q_tot = q_vap + q_liq + q_rai
            twp = np.sum(rho * q_tot * dz_broadcast, axis=0)
            dt = np.diff(t)
            dtwp = -np.diff(twp)
            rr = dtwp / dt
            rr = np.concatenate([[0], rr])

            # Apply moving average to smooth noise
            window_size = 10
            rr_smoothed = np.convolve(rr, np.ones(window_size)/window_size, mode='same')
            rr_smoothed[:window_size//2] = rr[:window_size//2]
            rr_smoothed[-window_size//2:] = rr[-window_size//2:]
            rr = rr_smoothed
    else:
        # KiD arrays are (time, z), so we sum along axis=1 (z-axis)
        # dz needs to be broadcast as (1, z)
        cwp = np.sum(rho * q_liq * dz[np.newaxis, :], axis=1)
        rwp = np.sum(rho * q_rai * dz[np.newaxis, :], axis=1)

        # Domain-averaged number concentrations (weighted by dz)
        total_height = np.sum(dz)
        Nc_avg = np.sum(N_liq * dz[np.newaxis, :], axis=1) / total_height * 1e-6  # m^-3 -> cm^-3
        Nr_avg = np.sum(N_rai * dz[np.newaxis, :], axis=1) / total_height * 1e-6  # m^-3 -> cm^-3

        if use_flux:
            # Try to compute rain rate from sedimentation source terms
            # Rain rate at surface = -integral(rho * (Sq_liq_prc + Sq_rai_prc) * dz)
            has_flux = False

            if 'Sq_liq_prc' in prof and 'Sq_rai_prc' in prof:
                Sq_liq = prof['Sq_liq_prc'].values  # (time, z), kg/kg/s
                Sq_rai = prof['Sq_rai_prc'].values  # (time, z), kg/kg/s

                # Check if source terms are valid (not all fill values)
                if np.max(np.abs(Sq_liq)) < 1e20 and np.max(np.abs(Sq_rai)) < 1e20:
                    # Sum of source terms
                    Sq_tot = Sq_liq + Sq_rai
                    # Rain rate = -integral(rho * Sq_tot * dz) over the column
                    # Negative because Sq_tot < 0 for loss due to sedimentation
                    rr = -np.sum(rho * Sq_tot * dz[np.newaxis, :], axis=1)  # kg m^-2 s^-1
                    has_flux = True

            # Fall back to -d(TWP)/dt if flux data not available
            if not has_flux:
                twp = np.sum(rho * q_tot * dz[np.newaxis, :], axis=1)
                dt = np.diff(t)
                dtwp = -np.diff(twp)
                rr = dtwp / dt
                rr = np.concatenate([[0], rr])
        else:
            # Compute rain rate as -d(TWP)/dt
            twp = np.sum(rho * q_tot * dz[np.newaxis, :], axis=1)
            dt = np.diff(t)
            dtwp = -np.diff(twp)
            rr = dtwp / dt
            rr = np.concatenate([[0], rr])

    return {
        't': t,  # Time in seconds
        'cwp': cwp,
        'rwp': rwp,
        'rr': rr * 3600,  # Rain rate in kg m^-2 hr^-1 (equivalent to mm/hr)
        'Nc_avg': Nc_avg,  # Domain-averaged cloud number concentration (cm^-3)
        'Nr_avg': Nr_avg,  # Domain-averaged rain number concentration (cm^-3)
    }


def compute_mean_effective_radius(dataset, is_pysdm=False):
    """
    Compute domain-averaged effective radius for cloud droplets.

    r_eff = 3 * q / (4 * pi * rho_water * N)

    Returns None if N_liq data is invalid (e.g., all NaN or zero)
    """
    rho_water = 1000.0  # kg/m^3

    if is_pysdm:
        t = dataset['time'].values
        q_liq = dataset['cloud water mixing ratio'].values  # kg/kg (already converted in notebook)
        N_liq = dataset['nc'].values  # m^-3
        rhod = dataset['rhod'].values
        q_vap = dataset['water_vapour_mixing_ratio'].values
        q_rai = dataset['rain water mixing ratio'].values  # kg/kg (already converted in notebook)
        N_rai = dataset['nr'].values

        # Total density
        rho = rhod * (1 + q_liq + q_rai + q_vap)

        # PySDM arrays are (z, time), transpose to (time, z) for consistency
        q_liq = q_liq.T + q_rai.T  # Total liquid mixing ratio
        N_liq = N_liq.T + N_rai.T  # Total liquid number concentration
        rho = rho.T
        sum_axis = 1  # Sum over z-axis after transpose

    else:
        t = dataset['t'].values
        q_liq = dataset['q_liq'].values + dataset['q_rai'].values  # (time, z)
        N_liq = dataset['N_liq'].values + dataset['N_rai'].values  # (time, z)
        rho = dataset['density'].values  # (time, z)
        sum_axis = 1  # Sum over z-axis

    # Check if N_liq data is valid
    if np.all(np.isnan(N_liq)) or np.all(N_liq == 0) or np.max(N_liq) < 1e3:
        return None  # Invalid N_liq data, skip this dataset

    # Compute effective radius where cloud exists
    # r_eff = [(3 * rho * q_liq) / (4 * pi * rho_water * N_liq)]^(1/3)
    mask = (q_liq > 1e-8) & (N_liq > 1e3)  # Only where cloud exists

    r_eff = np.zeros_like(q_liq)
    # Calculate volume per droplet, then take cube root to get radius
    r_eff[mask] = np.cbrt((3 * rho[mask] * q_liq[mask]) / (4 * np.pi * rho_water * N_liq[mask]))

    # Convert to micrometers
    r_eff *= 1e6

    # Compute mass-weighted mean effective radius
    mass = rho * q_liq
    mass[~mask] = 0
    r_eff[~mask] = 0

    mean_r_eff = np.sum(mass * r_eff, axis=sum_axis) / np.maximum(np.sum(mass, axis=sum_axis), 1e-20)

    return (t, mean_r_eff)  # time in seconds, radius in micrometers


def plot_time_series(datasets, figsize=(12, 9), tmax=1000):
    """
    Plot time series of CWP/RWP, domain-averaged number concentrations, and rain rate.

    Parameters:
    -----------
    datasets : dict
        Dictionary of xarray datasets with keys as dataset names
    """
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True, layout='constrained')

    colors = ['k', "tab:blue", "tab:green", "tab:orange", "tab:red", "tab:purple"]

    # Track handles and labels for custom legend
    cwp_handles = []
    rwp_handles = []
    Nc_handles = []
    Nr_handles = []
    dataset_labels = []

    # Create twin axis for rain number concentration
    ax1_twin = axes[1].twinx()

    # Plot for each dataset
    for k, (dataset_name, dataset) in enumerate(datasets.items()):
        is_pysdm = (dataset_name == "PySDM")

        # Compute quantities (use flux method for rain rate)
        quantities = compute_path_integrated_quantities(dataset, is_pysdm=is_pysdm, use_flux=False)
        t_min = quantities['t']

        # Panel 1: CWP and RWP
        ax = axes[0]
        line1, = ax.plot(t_min, quantities['cwp'], '-', color=colors[k],
                         linewidth=2, alpha=0.7)
        line2, = ax.plot(t_min, quantities['rwp'], '--', color=colors[k],
                         linewidth=2, alpha=0.7)

        # Store handles for legend
        if is_pysdm:
            # For PySDM, show both CWP and RWP
            cwp_handles.append(line1)
            rwp_handles.append(line2)
            dataset_labels.append(dataset_name)
        else:
            # For others, only show dataset name once
            cwp_handles.append(line1)
            dataset_labels.append(dataset_name)

        # Panel 2: Domain-averaged number concentrations (cloud on left, rain on right)
        ax = axes[1]
        line3, = ax.plot(t_min, quantities['Nc_avg'], '-', color=colors[k],
                         linewidth=2, alpha=0.7)
        line4, = ax1_twin.plot(t_min, quantities['Nr_avg'], '--', color=colors[k],
                              linewidth=2, alpha=0.7)

        Nc_handles.append(line3)
        Nr_handles.append(line4)

        # Panel 3: Rain rate
        ax = axes[2]
        ax.plot(t_min, quantities['rr'], '-', color=colors[k],
                linewidth=2, label=dataset_name)

    # Configure axes
    axes[0].set_ylabel('Water Path (kg m$^{-2}$)', fontsize=11)
    axes[0].set_ylim(bottom=0)  # Set ymin to 0
    axes[0].set_xlim(0, tmax)

    # Custom legend for panel 1 - dataset labels
    legend_handles = []
    legend_labels = []
    for i, (handle, label) in enumerate(zip(cwp_handles, dataset_labels)):
            # Create a compound handle showing both solid and dashed
            legend_handles.append((handle, plt.Line2D([0], [0], color=handle.get_color(),
                                                       linestyle='--', alpha=0.7)))
            legend_labels.append(label)

    # Add second legend for cloud/rain line styles
    cloud_line = plt.Line2D([0], [0], color='k', linestyle='-', linewidth=2)
    rain_line = plt.Line2D([0], [0], color='k', linestyle='--', linewidth=2)
    leg2 = axes[0].legend([cloud_line, rain_line], ['Cloud', 'Rain'],
                          loc='center left', bbox_to_anchor=(1, 0.2),
                          framealpha=0.9, fontsize=9, )#title='Category')

    leg1 = axes[0].legend(legend_handles, legend_labels,
                          loc='center left', bbox_to_anchor=(1, 0.7),
                          framealpha=0.9, fontsize=9, )#title='Dataset')
    
    axes[0].add_artist(leg2)  # Add the first legend back
    axes[0].grid(True, alpha=0.3)


    # Configure panel 2
    axes[1].set_ylabel('Total $N_c$ (cm$^{-3}$)', fontsize=11)
    axes[1].set_ylim(bottom=0)
    axes[1].set_xlim(0, tmax)
    axes[1].grid(True, alpha=0.3)

    # Configure right y-axis for rain
    ax1_twin.set_ylabel('Total $N_r$ (cm$^{-3}$)', fontsize=11)
    ax1_twin.set_ylim(bottom=0)

    # Configure panel 3
    axes[2].set_ylabel('Rain Rate (mm hr$^{-1}$)', fontsize=11)
    axes[2].set_xlabel('Time (s)', fontsize=11)
    axes[2].set_ylim(bottom=0)
    axes[2].set_xlim(0, tmax)
    axes[2].grid(True, alpha=0.3)

    return fig, axes


# ============================================================================
# 3. TIMEHEIGHT PLOTS FOR ALL DATASETS (APPENDIX)
# ============================================================================

def plot_timeheight_comparison(datasets, figsize=(15, 10), vmax_q=1.0, vmax_N=100):
    """
    Create timeheight plots for total liquid concentration and total number
    concentration for all datasets (for appendix).

    Parameters:
    -----------
    datasets : dict
        Dictionary of xarray datasets with keys as dataset names
    figsize : tuple
        Figure size (width, height)
    vmax_q : float
        Maximum value for mixing ratio colorbar (g/kg)
    vmax_N : float
        Maximum value for number concentration colorbar (cm^-3)

    Returns:
    --------
    fig, axes : matplotlib figure and axes
    """
    n_datasets = len(datasets)
    fig, axes = plt.subplots(n_datasets, 2, figsize=figsize, sharex=True)

    # Handle case of single dataset
    if n_datasets == 1:
        axes = axes.reshape(1, -1)

    for i, (dataset_name, dataset) in enumerate(datasets.items()):
        is_pysdm = (dataset_name == "PySDM")

        if is_pysdm:
            # PySDM data
            t = dataset['time'].values
            z = dataset['height'].values

            # Get mixing ratios (convert from )
            q_liq = dataset['cloud water mixing ratio'].values * 1000  # kg/kg -> g/kg
            q_rai = dataset['rain water mixing ratio'].values * 1000  # kg/kg -> g/kg
            q_tot = (q_liq + q_rai)  # Total liquid in g/kg

            # Get number concentrations (convert from m^-3 to cm^-3)
            N_liq = dataset['nc'].values * 1e-6  # m^-3 -> cm^-3
            N_rai = dataset['nr'].values * 1e-6  # m^-3 -> cm^-3
            N_tot = N_liq + N_rai  # Total number in cm^-3

            # Transpose for plotting (PySDM is (z, time), need (time, z))
            q_tot = q_tot.T
            N_tot = N_tot.T

        else:
            # KiD data
            t = dataset['t'].values
            z = dataset['zc'].values

            # Get mixing ratios (already in kg/kg, convert to g/kg)
            q_liq = dataset['q_liq'].values * 1000  # kg/kg -> g/kg
            q_rai = dataset['q_rai'].values * 1000  # kg/kg -> g/kg
            q_tot = q_liq + q_rai  # Total liquid in g/kg

            # Get number concentrations (convert from m^-3 to cm^-3)
            N_liq = dataset['N_liq'].values * 1e-6  # m^-3 -> cm^-3
            N_rai = dataset['N_rai'].values * 1e-6  # m^-3 -> cm^-3
            N_tot = N_liq + N_rai  # Total number in cm^-3

        # Create meshgrid for pcolormesh
        T, Z = np.meshgrid(t, z)

        # Plot total mixing ratio
        ax_q = axes[i, 0]
        im_q = ax_q.pcolormesh(T, Z, q_tot.T, shading='auto',
                                cmap='BuPu', vmin=0, vmax=vmax_q)
        ax_q.set_ylabel(f'{dataset_name}\nAltitude (m)', fontsize=10)
        if i == n_datasets - 1:
            ax_q.set_xlabel('Time (min)', fontsize=10)

        # Add colorbar for mixing ratio
        cbar_q = plt.colorbar(im_q, ax=ax_q, pad=0.02)
        cbar_q.set_label('$q_c + q_r$ (g kg$^{-1}$)', fontsize=9)

        # Plot total number concentration
        ax_N = axes[i, 1]
        im_N = ax_N.pcolormesh(T, Z, N_tot.T, shading='auto',
                                cmap='BuPu', vmin=0, vmax=vmax_N)
        if i == n_datasets - 1:
            ax_N.set_xlabel('Time (min)', fontsize=10)

        # Add colorbar for number concentration
        cbar_N = plt.colorbar(im_N, ax=ax_N, pad=0.02)
        cbar_N.set_label('$N_c + N_r$ (cm$^{-3}$)', fontsize=9)

    # Add column titles
    axes[0, 0].set_title('Total Liquid Mixing Ratio', fontsize=12, fontweight='bold', pad=10)
    axes[0, 1].set_title('Total Liquid Number Concentration', fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()
    return fig, axes


def plot_timeheight_individual(dataset, dataset_name, is_pysdm=False,
                                figsize=(12, 4), vmax_q=1.0, vmax_N=100):
    """
    Create a single timeheight plot for one dataset showing total liquid
    concentration and total number concentration.

    Parameters:
    -----------
    dataset : xarray.Dataset
        Dataset to plot
    dataset_name : str
        Name of the dataset (for title)
    is_pysdm : bool
        Whether this is PySDM data (different coordinate convention)
    figsize : tuple
        Figure size (width, height)
    vmax_q : float
        Maximum value for mixing ratio colorbar (g/kg)
    vmax_N : float
        Maximum value for number concentration colorbar (cm^-3)

    Returns:
    --------
    fig, axes : matplotlib figure and axes
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    if is_pysdm:
        # PySDM data
        t = dataset['time'].values / 60  # Convert to minutes
        z = dataset['height'].values

        # Get mixing ratios
        q_liq = dataset['cloud water mixing ratio'].values / 1000  # g/kg -> kg/kg
        q_rai = dataset['rain water mixing ratio'].values / 1000  # g/kg -> kg/kg
        q_tot = (q_liq + q_rai) * 1000  # Total liquid in g/kg

        # Get number concentrations
        N_liq = dataset['nc'].values * 1e-6  # m^-3 -> cm^-3
        N_rai = dataset['nr'].values * 1e-6  # m^-3 -> cm^-3
        N_tot = N_liq + N_rai

        # Transpose for plotting
        q_tot = q_tot.T
        N_tot = N_tot.T

    else:
        # KiD data
        t = dataset['t'].values / 60  # Convert to minutes
        z = dataset['zc'].values

        # Get mixing ratios
        q_liq = dataset['q_liq'].values * 1000  # kg/kg -> g/kg
        q_rai = dataset['q_rai'].values * 1000  # kg/kg -> g/kg
        q_tot = q_liq + q_rai

        # Get number concentrations
        N_liq = dataset['N_liq'].values * 1e-6  # m^-3 -> cm^-3
        N_rai = dataset['N_rai'].values * 1e-6  # m^-3 -> cm^-3
        N_tot = N_liq + N_rai

    # Create meshgrid
    T, Z = np.meshgrid(t, z)

    # Plot total mixing ratio
    ax = axes[0]
    im_q = ax.pcolormesh(T, Z, q_tot.T, shading='auto',
                          cmap='BuPu', vmin=0, vmax=vmax_q)
    ax.set_xlabel('Time (min)', fontsize=11)
    ax.set_ylabel('Altitude (m)', fontsize=11)
    ax.set_title('Total Liquid Mixing Ratio', fontsize=12, fontweight='bold')
    cbar_q = plt.colorbar(im_q, ax=ax)
    cbar_q.set_label('$q_c + q_r$ (g kg$^{-1}$)', fontsize=10)

    # Plot total number concentration
    ax = axes[1]
    im_N = ax.pcolormesh(T, Z, N_tot.T, shading='auto',
                          cmap='BuPu', vmin=0, vmax=vmax_N)
    ax.set_xlabel('Time (min)', fontsize=11)
    ax.set_title('Total Liquid Number Concentration', fontsize=12, fontweight='bold')
    cbar_N = plt.colorbar(im_N, ax=ax)
    cbar_N.set_label('$N_c + N_r$ (cm$^{-3}$)', fontsize=10)

    fig.suptitle(dataset_name, fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    return fig, axes
