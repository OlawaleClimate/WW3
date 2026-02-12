"""
WaveWatch Python Package

A Python implementation of wave parameter computation from action density spectrum,
based on WW3 (WaveWatch III) model algorithms.

This package provides tools to:
  - Compute spectral moments from action density spectrum
  - Solve dispersion relations for wavenumber and group velocity
  - Calculate wave parameters: HS, T01, T02, T0M1, mean direction, etc.
  - Handle Pierson-Moskowitz tail frequency extension

Basic Usage:
    >>> import numpy as np
    >>> from wavewatch import WaveSpectrum
    >>>
    >>> # Create action density spectrum (ndir, nfreq)
    >>> action = np.random.rand(36, 30) * 0.01
    >>> depth = 100.0  # Water depth [m]
    >>>
    >>> # Initialize spectrum
    >>> spectrum = WaveSpectrum(action, depth=depth)
    >>>
    >>> # Compute wave parameters
    >>> params = spectrum.compute_parameters()
    >>>
    >>> print(f"Significant Wave Height: {params['hs']:.2f} m")
    >>> print(f"Mean Period: {params['t01']:.2f} s")
    >>> print(f"Mean Direction: {np.degrees(params['thm']):.1f}°")

References:
    - WaveWatch III (WW3) User Manual
    - WW3 source code: w3iogomd.F90, w3dispmd.F90
    - Tolman, H. L., et al. (2002). Development and implementation of windwave models...

Author: Claude (based on WW3 algorithms)
License: MIT
"""

__version__ = "0.1.0"
__author__ = "Claude"

from .spectrum import WaveSpectrum
from .dispersion import solve_dispersion, compute_wavelength, compute_group_velocity_deep_water
from .spectrum_converter import (
    build_ww3_grid,
    action_to_energy_2d,
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d,
    get_peak_frequency,
    get_peak_direction,
    normalize_spectrum
)
from .unified_converter import SpectrumConverter
from .constants import (
    GRAV, TPI, TPIINV, RADE, DEGRAD,
    DEFAULT_NK, DEFAULT_NTH, DEFAULT_FR1, DEFAULT_XFR
)

__all__ = [
    'SpectrumConverter',
    'WaveSpectrum',
    'solve_dispersion',
    'compute_wavelength',
    'compute_group_velocity_deep_water',
    'build_ww3_grid',
    'action_to_energy_2d',
    'action_to_frequency_spectrum_1d',
    'action_to_directional_spectrum_1d',
    'get_peak_frequency',
    'get_peak_direction',
    'normalize_spectrum',
    'GRAV',
    'TPI',
    'TPIINV',
    'RADE',
    'DEGRAD',
    'DEFAULT_NK',
    'DEFAULT_NTH',
    'DEFAULT_FR1',
    'DEFAULT_XFR',
]
