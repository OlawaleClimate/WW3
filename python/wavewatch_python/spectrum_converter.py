"""
Spectrum conversion utilities for WaveWatch III.

Converts between different spectral representations:
- Action density to energy density (frequency-direction)
- Action density to 1D frequency spectrum
- Action density to 1D directional spectrum

This implementation mirrors the Fortran w3iogomd.F90 and w3gridmd.F90 algorithms,
properly handling the coordinate transformation from action density A(k,θ) in
wavenumber-direction space to energy density E(f,θ) in frequency-direction space.
"""

import numpy as np
from .constants import GRAV, TPI, SMALL


def build_ww3_grid(fr1, xfr, nk, nth):
    """
    Build WW3 spectral grid arrays (mirrors w3gridmd.F90).

    Constructs the frequency and directional grids with proper frequency
    bandwidth (DSII) computation including special cases for first and last bins.

    Arguments:
        fr1 (float): First frequency [Hz]
        xfr (float): Frequency increment factor (geometric ratio)
        nk (int): Number of frequency bins
        nth (int): Number of directional bins

    Returns:
        dict with keys:
            'freq' : Frequency array [Hz]
            'sigma' : Angular frequency array SIG [rad/s]
            'theta' : Direction array [radians]
            'dth' : Directional bin width [radians]
            'dsii' : Frequency bandwidth array [rad/s]
            'dden' : DDEN = DTH × DSII × SIG (conversion factor)
            'fte' : Tail energy factor
    """
    dth = TPI / nth

    # Directions: evenly spaced, centred
    theta = np.array([(ith - 0.5) * dth for ith in range(1, nth + 1)])

    # Frequencies: geometric series (mimics Fortran w3gridmd.F90)
    sigma_start = fr1 * TPI / (xfr**2)
    sigma = np.zeros(nk)
    for ik in range(nk):
        sigma_start *= xfr
        sigma[ik] = sigma_start

    freq = sigma / TPI

    # Frequency bandwidths DSII
    sxfr = 0.5 * (xfr - 1.0 / xfr)
    dsii = sigma * sxfr

    # Special cases for first and last bins (w3gridmd.F90 lines 1305-1306)
    dsii[0] = 0.5 * sigma[0] * (xfr - 1.0)
    dsii[-1] = 0.5 * sigma[-1] * (xfr - 1.0) / xfr

    # DDEN: combined factor for discrete integration
    dden = dth * dsii * sigma

    # Tail energy factor
    fte = 0.25 * sigma[-1] * dth * sigma[-1]

    return {
        'freq': freq,
        'sigma': sigma,
        'theta': theta,
        'dth': dth,
        'dsii': dsii,
        'dden': dden,
        'fte': fte,
    }


def action_to_energy_2d(action, omega, group_velocity, dintegral=None,
                        directions=None, ddir=None, dsii=None):
    """
    Convert 2D action density spectrum to 2D energy density spectrum.

    Converts from action density A(k,θ) to energy density E(f,θ) using:
      E(f,θ) = A(k,θ) × σ × (∂k/∂f)
              = A(k,θ) × SIG × (2π / CG)

    Where:
      - A(k,θ) = F(k,θ) / σ (action = energy / intrinsic frequency)
      - σ = SIG = intrinsic (angular) frequency [rad/s]
      - ∂k/∂f = 2π / CG = Jacobian of (k,θ) → (f,θ) transformation
      - SIG is the angular frequency array from WW3 w3gridmd.F90

    This is a pure spectral conversion with NO bin-width factors included.
    Bin widths (DTH, DSII) are applied only when integrating to 1D spectra.

    Arguments:
        action (ndarray): Action density spectrum (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s], length nfreq
        group_velocity (ndarray): Group velocity [m/s], length nfreq
        dintegral (ndarray, optional): DDEN factors (not used in core conversion)
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]
        dsii (ndarray, optional): Frequency bandwidths [rad/s]

    Returns:
        tuple:
            - energy_2d (ndarray): Energy density (ndir, nfreq) [m²/Hz/rad]
            - frequencies (ndarray): Frequency array [Hz]
            - directions (ndarray): Direction array [radians]
            - grid_params (dict): Grid parameters including dth, dsii, dden
    """

    # Validate and reshape action
    action = np.atleast_2d(np.asarray(action, dtype=float))
    if action.shape[0] < action.shape[1]:
        action = action.T

    ndir, nfreq = action.shape
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    group_velocity = np.atleast_1d(np.asarray(group_velocity, dtype=float))

    # Validate dimensions
    if len(omega) != nfreq:
        raise ValueError(f"omega length ({len(omega)}) != nfreq ({nfreq})")
    if len(group_velocity) != nfreq:
        raise ValueError(f"group_velocity length ({len(group_velocity)}) != nfreq ({nfreq})")

    # Generate or use provided direction grid
    if directions is None:
        directions = np.linspace(0, 2*np.pi, ndir, endpoint=False)
    else:
        directions = np.atleast_1d(np.asarray(directions, dtype=float))

    # Compute or use provided directional spacing
    if ddir is None:
        ddir = 2.0 * np.pi / ndir

    # Compute or use provided DSII (frequency bandwidths)
    if dsii is None:
        sxfr = 0.5 * (1.1 - 1.0/1.1)  # Default XFR = 1.1
        dsii_array = omega * sxfr
    else:
        dsii_array = np.atleast_1d(np.asarray(dsii, dtype=float))

    # DDEN factor for later use
    dden = ddir * dsii_array * omega

    # CORRECT CONVERSION FORMULA (no bin widths):
    # E(f,θ) = A(k,θ) × σ × (2π / CG)
    # Where: σ = SIG = intrinsic (angular) frequency = omega
    # This combines the action-to-energy relationship (×σ) with the Jacobian (×2π/CG)
    conversion_factor = omega * (2.0 * np.pi) / (group_velocity + SMALL)

    # Apply conversion factor to all directions at each frequency
    energy_2d = action * conversion_factor[np.newaxis, :]

    # Convert omega to frequency
    frequencies = omega / (2.0 * np.pi)

    # Return grid parameters for use in 1D integration
    grid_params = {
        'dth': ddir,
        'dsii': dsii_array,
        'dden': dden,
    }

    return energy_2d, frequencies, directions, grid_params


def action_to_frequency_spectrum_1d(action, omega, group_velocity, dintegral=None,
                                     directions=None, ddir=None, dsii=None):
    """
    Convert 2D action density to 1D frequency spectrum.

    Integrates action density over all directions to produce a frequency spectrum E(f).
    Follows w3iogomd.F90 frequency spectrum computation.

    The conversion is:
      E(f) = ∑_θ E(f,θ) × Δθ = ∑_θ [A(θ,f) × SIG × (2π/CG)] × DTH

    Arguments:
        action (ndarray): Action density (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        dintegral (ndarray, optional): DDEN factors (w3gridmd DDEN array)
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]
        dsii (ndarray, optional): Frequency bandwidths [rad/s]

    Returns:
        tuple:
            - freq (ndarray): Frequency array [Hz]
            - e_freq (ndarray): 1D frequency spectrum [m²/Hz]
    """

    # Get 2D energy spectrum and grid parameters
    energy_2d, freq, dirs, grid_params = action_to_energy_2d(
        action, omega, group_velocity, dintegral, directions, ddir, dsii
    )

    dth = grid_params['dth']

    # Integrate over directions with directional bin width
    # e_freq = ∑_θ E(f,θ) × DTH
    e_freq = np.sum(energy_2d, axis=0) * dth

    return freq, e_freq


def action_to_directional_spectrum_1d(action, omega, group_velocity, dintegral=None,
                                       directions=None, ddir=None, dsii=None):
    """
    Convert 2D action density to 1D directional spectrum.

    Integrates action density over all frequencies to produce a directional spectrum E(θ).
    Follows w3iogomd.F90 directional spectrum computation.

    The conversion is:
      E(θ) = ∑_f [A(θ,f) × SIG × (2π/CG) × DSII] = ∑_f [A(θ,f) × DSII × SIG / CG]

    Arguments:
        action (ndarray): Action density (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        dintegral (ndarray, optional): DDEN factors (w3gridmd DDEN array)
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]
        dsii (ndarray, optional): Frequency bandwidths [rad/s]

    Returns:
        tuple:
            - dirs (ndarray): Direction array [radians]
            - e_dir (ndarray): 1D directional spectrum [m²/rad]
    """

    # Validate and reshape action
    action = np.atleast_2d(np.asarray(action, dtype=float))
    if action.shape[0] < action.shape[1]:
        action = action.T

    ndir, nfreq = action.shape
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    group_velocity = np.atleast_1d(np.asarray(group_velocity, dtype=float))

    # Validate dimensions
    if len(omega) != nfreq:
        raise ValueError(f"omega length ({len(omega)}) != nfreq ({nfreq})")
    if len(group_velocity) != nfreq:
        raise ValueError(f"group_velocity length ({len(group_velocity)}) != nfreq ({nfreq})")

    # Generate or use provided direction grid
    if directions is None:
        directions = np.linspace(0, 2*np.pi, ndir, endpoint=False)
    else:
        directions = np.atleast_1d(np.asarray(directions, dtype=float))

    # Compute or use provided DSII (frequency bandwidths)
    if dsii is None:
        sxfr = 0.5 * (1.1 - 1.0/1.1)  # Default XFR = 1.1
        dsii_array = omega * sxfr
    else:
        dsii_array = np.atleast_1d(np.asarray(dsii, dtype=float))

    # CORRECT CONVERSION FORMULA with frequency bin width:
    # E(θ) = ∑_f [A(θ,f) × SIG × (2π/CG) × DSII]
    conversion_factor = omega * (2.0 * np.pi) / (group_velocity + SMALL)

    # Integrate over frequencies with frequency bin width
    # e_dir = ∑_f [action * conversion_factor * DSII]
    e_dir = np.sum(action * conversion_factor[np.newaxis, :] * dsii_array[np.newaxis, :], axis=1)

    return directions, e_dir


def normalize_spectrum(energy_2d, target_energy=None):
    """
    Normalize 2D energy spectrum to specific total energy.

    Arguments:
        energy_2d (ndarray): 2D energy spectrum (ndir, nfreq)
        target_energy (float, optional): Target total energy [m²]
                                        If None, returns normalized spectrum

    Returns:
        ndarray: Normalized spectrum
    """

    total_energy = np.sum(energy_2d)

    if total_energy < SMALL:
        return energy_2d

    normalized = energy_2d / total_energy

    if target_energy is not None:
        normalized = normalized * target_energy

    return normalized


def get_peak_frequency(action_2d, omega, group_velocity, directions=None, ddir=None):
    """
    Get peak frequency from action density spectrum.

    Finds the frequency with maximum energy in the frequency spectrum.

    Arguments:
        action_2d (ndarray): 2D action spectrum (ndir, nfreq)
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]

    Returns:
        float: Peak frequency [Hz]
    """

    # Get 1D frequency spectrum
    freq, e_freq = action_to_frequency_spectrum_1d(
        action_2d, omega, group_velocity, directions=directions, ddir=ddir
    )

    # Find peak
    peak_idx = np.argmax(e_freq)

    return freq[peak_idx]


def get_peak_direction(action_2d, omega, group_velocity, directions=None, ddir=None):
    """
    Get peak direction from action density spectrum.

    Finds the direction with maximum energy in the directional spectrum.

    Arguments:
        action_2d (ndarray): 2D action spectrum (ndir, nfreq)
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]

    Returns:
        float: Peak direction [radians]
    """

    # Get 1D directional spectrum
    dirs, e_dir = action_to_directional_spectrum_1d(
        action_2d, omega, group_velocity, directions=directions, ddir=ddir
    )

    # Find peak
    peak_idx = np.argmax(e_dir)

    return dirs[peak_idx]
