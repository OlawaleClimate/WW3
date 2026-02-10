"""
Physical and mathematical constants for wave spectrum calculations.

Based on WW3 (WaveWatch III) model constants.
"""

import numpy as np

# Physical Constants
GRAV = 9.81                    # Gravitational acceleration (m/s²)
DWAT = 1025.0                  # Water density (kg/m³)
DAIR = 1.225                   # Air density (kg/m³)

# Mathematical Constants
PI = np.pi
TPI = 2.0 * PI                 # 2π
TPIINV = 1.0 / TPI             # 1/(2π)
RADE = 180.0 / PI              # Radians to degrees conversion
DEGRAD = PI / 180.0            # Degrees to radians conversion

# Small numbers for numerical stability
SMALL = 1.0e-15
UNDEF = -999.0                 # Undefined value marker

# Default model parameters
DEFAULT_NK = 30                # Number of frequency bins
DEFAULT_NTH = 36               # Number of directional bins
DEFAULT_FR1 = 0.04             # First frequency (Hz)
DEFAULT_XFR = 1.1              # Frequency ratio
DEFAULT_DMIN = 0.5             # Minimum water depth (m)
