"""
Dispersion relation solver for converting frequency to wavenumber.

Based on WW3's WAVNU1 subroutine.
Reference: w3dispmd.F90
"""

import numpy as np
from .constants import GRAV, SMALL

def solve_dispersion(sigma, depth, max_iterations=100, tolerance=1e-10):
    """
    Solve dispersion relation for wavenumber and group velocity.

    Uses Newton-Raphson iteration to solve:
        σ² = g·k·tanh(k·h)

    Arguments:
        sigma (float or array): Angular frequency [rad/s]
        depth (float or array): Water depth [m]
        max_iterations (int): Maximum Newton-Raphson iterations
        tolerance (float): Convergence tolerance

    Returns:
        wn (float or array): Wavenumber [1/m]
        cg (float or array): Group velocity [m/s]

    References:
        - WW3 w3dispmd.F90 (WAVNU1 subroutine)
        - Eckart, C. (1953). The laboratory generation and propagation of ocean waves
    """

    # Handle scalar and array inputs
    sigma = np.atleast_1d(np.asarray(sigma, dtype=float))
    depth = np.atleast_1d(np.asarray(depth, dtype=float))

    # Broadcast to common shape if needed
    if sigma.shape != depth.shape:
        sigma, depth = np.broadcast_arrays(sigma, depth)

    # Initialize wavenumber with shallow water approximation
    # k ≈ σ/√(g·h)
    wn = sigma / np.sqrt(GRAV * depth)

    # Newton-Raphson iteration
    for iteration in range(max_iterations):
        # Compute tanh(k·h)
        kh = wn * depth

        # Limit kh to avoid numerical issues
        kh_limited = np.minimum(kh, 100.0)

        # Compute tanh and its derivative
        sinh_kh = np.sinh(kh_limited)
        cosh_kh = np.cosh(kh_limited)
        tanh_kh = sinh_kh / cosh_kh

        # Derivative of tanh: 1 - tanh²
        dtanh_dkh = 1.0 - tanh_kh**2

        # Dispersion relation residual: σ² - g·k·tanh(k·h)
        sigma_squared = sigma**2
        f = sigma_squared - GRAV * wn * tanh_kh

        # Jacobian: df/dk = -g·(tanh + k·h·(1-tanh²))
        df_dwn = -GRAV * (tanh_kh + wn * depth * dtanh_dkh)

        # Avoid division by zero
        df_dwn = np.where(np.abs(df_dwn) < SMALL, SMALL, df_dwn)

        # Newton-Raphson update
        dwn = f / df_dwn
        wn_new = wn - dwn

        # Check convergence
        max_change = np.max(np.abs(dwn / (wn + SMALL)))
        if max_change < tolerance:
            break

        wn = wn_new

    # Compute group velocity
    # In general: cg = ∂σ/∂k = (σ/2k) × (1 + 2kh/sinh(2kh))
    # Simplified for deep/shallow: cg = σ/(2k) for deep, cg = √(gh) for shallow

    kh = wn * depth
    kh_limited = np.minimum(kh, 100.0)

    # Deep water limit (kh > π): cg = σ/(2k)
    # Shallow water limit (kh < π/20): cg = √(gh)

    sinh_2kh = np.sinh(2 * kh_limited)

    # Full formula for group velocity
    cg = (sigma / (2 * wn)) * (1.0 + (2.0 * kh_limited / sinh_2kh))

    # Return scalars if inputs were scalar
    if wn.size == 1:
        wn = float(wn[0])
        cg = float(cg[0])

    return wn, cg


def compute_wavelength(wn):
    """
    Compute wavelength from wavenumber.

    λ = 2π/k

    Arguments:
        wn (float or array): Wavenumber [1/m]

    Returns:
        wavelength (float or array): Wavelength [m]
    """
    from .constants import TPI
    return TPI / wn


def compute_group_velocity_deep_water(sigma):
    """
    Compute group velocity in deep water (simplified).

    cg = σ/(2k) = √(g/k) = √(gσ/2π)

    Arguments:
        sigma (float or array): Angular frequency [rad/s]

    Returns:
        cg (float or array): Group velocity [m/s]
    """
    from .constants import TPI
    return np.sqrt(GRAV * sigma / TPI) / 2.0
