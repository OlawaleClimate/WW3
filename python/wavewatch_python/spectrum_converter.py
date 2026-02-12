"""
Spectrum conversion utilities for WaveWatch III.

Converts between different spectral representations:
- Action density to energy density (frequency-direction)
- Action density to 1D frequency spectrum
- Action density to 1D directional spectrum
"""

import numpy as np
from .constants import GRAV, TPI, SMALL


def action_to_energy_2d(action, omega, group_velocity, dintegral=None,
                        directions=None, ddir=None):
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

    Reference: Wave action balance with coordinate transformation

    Arguments:
        action (ndarray): Action density spectrum (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s], length nfreq
        group_velocity (ndarray): Group velocity [m/s], length nfreq
        dintegral (ndarray, optional): Integration factors DDEN = DTH × DSII × SIG
                                       If not provided, computed from omega and ddir
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]
                               If not provided, computed from ndir

    Returns:
        tuple:
            - energy_2d (ndarray): Energy density (ndir, nfreq) [m²/Hz/rad]
            - frequencies (ndarray): Frequency array [Hz]
            - directions (ndarray): Direction array [radians]
            - dintegral_used (ndarray): Integration factors used
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

    # CORRECT CONVERSION FORMULA:
    # E(f,θ) = A(k,θ) × σ × (2π / CG)
    # Where: σ = SIG = intrinsic (angular) frequency = omega
    # This combines the action-to-energy relationship (×σ) with the Jacobian (×2π/CG)
    conversion_factor = omega * (2.0 * np.pi) / (group_velocity + SMALL)

    # Apply conversion factor to all directions at each frequency
    energy_2d = action * conversion_factor[np.newaxis, :]

    # Convert omega to frequency
    frequencies = omega / (2.0 * np.pi)

    return energy_2d, frequencies, directions, group_velocity


def action_to_frequency_spectrum_1d(action, omega, group_velocity, dintegral=None,
                                     directions=None, ddir=None):
    """
    Convert 2D action density to 1D frequency spectrum.

    Integrates action density over all directions to produce a frequency spectrum.

    Note: Since DDEN already includes DTH (directional binning), we only sum
    over directions without additional multiplication by ddir.

    Arguments:
        action (ndarray): Action density (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        dintegral (ndarray, optional): Integration factors
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]

    Returns:
        tuple:
            - freq (ndarray): Frequency array [Hz]
            - e_freq (ndarray): 1D energy spectrum [m²/Hz]
    """

    # Get 2D energy spectrum first
    energy_2d, freq, _, _ = action_to_energy_2d(
        action, omega, group_velocity, dintegral, directions, ddir
    )

    # Integrate over directions (sum all directions)
    # DTH is already included in DDEN, so we just sum without extra multiplication
    e_freq = np.sum(energy_2d, axis=0)

    return freq, e_freq


def action_to_directional_spectrum_1d(action, omega, group_velocity, dintegral=None,
                                       directions=None, ddir=None):
    """
    Convert 2D action density to 1D directional spectrum.

    Integrates action density over all frequencies to produce a directional spectrum.

    Uses the same conversion formula as 2D energy:
      E(θ) = ∑_f [A(θ,f) × SIG × (2π / CG)]

    Arguments:
        action (ndarray): Action density (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        dintegral (ndarray, optional): Deprecated - not used
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]

    Returns:
        tuple:
            - dirs (ndarray): Direction array [radians]
            - e_dir (ndarray): 1D energy spectrum [m²/rad]
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

    # CORRECT CONVERSION FORMULA:
    # E(θ) = ∑_f [A(θ,f) × σ × (2π / CG)]
    # Where: σ = SIG = intrinsic (angular) frequency = omega
    conversion_factor = omega * (2.0 * np.pi) / (group_velocity + SMALL)

    # Integrate over frequencies (sum all frequencies with conversion applied)
    e_dir = np.sum(action * conversion_factor[np.newaxis, :], axis=1)

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


def get_peak_frequency(action_2d, frequencies):
    """
    Get peak frequency from action density spectrum.

    Arguments:
        action_2d (ndarray): 2D action spectrum (ndir, nfreq)
        frequencies (ndarray): Frequency array [Hz]

    Returns:
        float: Peak frequency [Hz]
    """

    # Integrate over directions to get 1D frequency spectrum
    freq_spectrum = np.sum(action_2d, axis=0)

    # Find peak
    peak_idx = np.argmax(freq_spectrum)

    return frequencies[peak_idx]


def get_peak_direction(action_2d, directions):
    """
    Get peak direction from action density spectrum.

    Arguments:
        action_2d (ndarray): 2D action spectrum (ndir, nfreq)
        directions (ndarray): Direction array [radians]

    Returns:
        float: Peak direction [radians]
    """

    # Integrate over frequencies to get 1D directional spectrum
    dir_spectrum = np.sum(action_2d, axis=1)

    # Find peak
    peak_idx = np.argmax(dir_spectrum)

    return directions[peak_idx]
