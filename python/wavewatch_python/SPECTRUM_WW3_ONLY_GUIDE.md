# WaveWatch Python: WW3-Only Streamlined Spectrum Module

## Overview

`spectrum_ww3_only.py` is a clean, production-ready implementation of the `WaveSpectrum` class that **only accepts pre-computed WW3 grid parameters**. This version removes all auto-generation logic and simplifies the API for users working directly with WW3 simulations.

## When to Use

### Use `spectrum_ww3_only.py` if:
- ✅ You're extracting action density from actual WW3 simulations
- ✅ You have access to WW3 grid parameters (SIG, DDEN, WN, CG, tail factors)
- ✅ You want a clean, minimal dependency implementation
- ✅ You value code clarity and maintainability over flexibility
- ✅ You're building production systems that integrate with WW3

### Use the standard `spectrum.py` if:
- ✅ You want auto-generation of grid parameters
- ✅ You're working with synthetic spectra (testing, research)
- ✅ You prefer a more flexible, backwards-compatible API
- ✅ You need both WW3 and non-WW3 modes in one module

## API Reference

### Constructor

```python
WaveSpectrum(action, depth, omega, dintegral, fte, fttr, ftwl,
             wavenumber, group_velocity, directions=None)
```

#### Required Arguments

| Parameter | Type | Source | Description |
|-----------|------|--------|-------------|
| `action` | ndarray | WW3 output | Action density (ndir, nfreq) [m²·s·rad⁻¹] |
| `depth` | float | Simulation | Water depth [m] |
| `omega` | ndarray | w3gridmd.F90 | Angular frequencies (SIG array) [rad/s] |
| `dintegral` | ndarray | w3gridmd.F90 | Integration factors DDEN = DTH × DSII × SIG |
| `fte` | float | w3gridmd.F90 | Energy tail factor (Pierson-Moskowitz) |
| `fttr` | float | w3gridmd.F90 | Period tail factor |
| `ftwl` | float | w3gridmd.F90 | Wavelength tail factor |
| `wavenumber` | ndarray | WAVNU1 subroutine | Wavenumber WN [1/m] |
| `group_velocity` | ndarray | WAVNU1 subroutine | Group velocity CG [m/s] |

#### Optional Arguments

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `directions` | ndarray | Uniform | Direction grid [radians] |

#### Raises

- `ValueError`: If array dimensions don't match
- `TypeError`: If required parameters are missing

### Methods

#### `integrate_moments()`

Computes spectral moments from action density.

```python
moments = spectrum.integrate_moments()
```

**Returns:** Dictionary with keys:
- `m0`: Total energy [m²]
- `m1`: First moment (frequency-weighted) [m²·s]
- `m2`: Second moment (frequency²-weighted) [m²·s²]
- `m_1`: Inverse moment [m²·s]
- `mom_x`: E-W directional moment (cosine)
- `mom_y`: N-S directional moment (sine)
- `mom_wn`: Wavelength moment

#### `compute_parameters()`

Derives wave parameters from spectral moments.

```python
params = spectrum.compute_parameters()
```

**Returns:** Dictionary with keys:
- `hs`: Significant wave height [m]
- `t01`: Mean period [s]
- `t02`: Zero-crossing period [s]
- `t0m1`: Energy period [s]
- `tp`: Peak period [s]
- `thm`: Mean direction [radians]
- `ths`: Directional spread [radians]
- `wlm`: Mean wavelength [m]
- `width`: Spectral width (dimensionless)
- `moments`: Spectral moments (dictionary)
- `depth`: Water depth [m]

#### `get_spectrum_1d()`

Extracts 1D frequency spectrum by integrating over directions.

```python
freq, e_freq = spectrum.get_spectrum_1d()
```

**Returns:**
- `freq`: Frequency array [Hz]
- `e_freq`: 1D energy spectrum [m²/Hz]

#### `get_spectrum_directional_1d()`

Extracts 1D directional spectrum by integrating over frequencies.

```python
dirs, e_dir = spectrum.get_spectrum_directional_1d()
```

**Returns:**
- `dirs`: Direction array [radians]
- `e_dir`: 1D energy spectrum [m²/rad]

## Usage Examples

### Basic Usage

```python
import numpy as np
from wavewatch_python.spectrum_ww3_only import WaveSpectrum
from wavewatch_python.dispersion import solve_dispersion

# Extract action density from WW3 simulation
action = np.array([...])  # shape (36, 30) from WW3 output
depth = 100.0

# Get pre-computed WW3 parameters
omega_ww3 = np.array([...])  # from w3gridmd.F90 SIG
dden_ww3 = np.array([...])   # from w3gridmd.F90 DDEN
fte_ww3 = 0.25 * omega_ww3[-1]**2 * (2*np.pi/36)
fttr_ww3 = 0.20 * (2*np.pi/36) * omega_ww3[-1]
ftwl_ww3 = (9.81 / 6.0) / omega_ww3[-1] * (2*np.pi/36) * omega_ww3[-1]

# Compute wavenumber and group velocity
wn_ww3 = np.zeros(30)
cg_ww3 = np.zeros(30)
for ik in range(30):
    wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

# Create spectrum object
spectrum = WaveSpectrum(
    action=action,
    depth=depth,
    omega=omega_ww3,
    dintegral=dden_ww3,
    fte=fte_ww3,
    fttr=fttr_ww3,
    ftwl=ftwl_ww3,
    wavenumber=wn_ww3,
    group_velocity=cg_ww3
)

# Compute wave parameters
params = spectrum.compute_parameters()
print(f"HS: {params['hs']:.2f} m")
print(f"T01: {params['t01']:.2f} s")
print(f"Mean Direction: {np.degrees(params['thm']):.1f}°")
```

### Extracting Action Density from WW3

```python
# WW3 NetCDF output
import netCDF4 as nc

# Load WW3 output
ww3_file = nc.Dataset('ww3_output.nc')

# Extract action density (varies by WW3 version)
action = ww3_file.variables['efth'][0, :, :, :]  # [time, lat, lon, freq, dir]
depth = ww3_file.variables['dpt'][0, :, :]

# Get grid parameters from WW3 initialization
omega_ww3 = 2*np.pi * ww3_file.variables['frequency'][:]
direction_ww3 = ww3_file.variables['direction'][:]
```

### Multi-Point Analysis

```python
# Analyze spectrum at multiple locations
nlat, nlon = action.shape[:2]

for i in range(nlat):
    for j in range(nlon):
        spectrum = WaveSpectrum(
            action=action[i, j, :, :],
            depth=depth[i, j],
            omega=omega_ww3,
            dintegral=dden_ww3,
            fte=fte_ww3,
            fttr=fttr_ww3,
            ftwl=ftwl_ww3,
            wavenumber=wn_ww3,
            group_velocity=cg_ww3
        )
        params = spectrum.compute_parameters()

        # Store results
        hs_grid[i, j] = params['hs']
        t01_grid[i, j] = params['t01']
        thm_grid[i, j] = params['thm']
```

## Wave Parameter Definitions

All parameters are computed following the spectral moment method used in WW3:

### Heights and Periods
- **HS (Significant Wave Height)**: HS = 4√m₀
  - Average of highest 1/3 of wave heights
  - Most important parameter for operational wave forecasting

- **T01 (Mean Period)**: T01 = 2π(m₀/m₁)
  - Average period of all waves in spectrum

- **T02 (Zero-crossing Period)**: T02 = 2π√(m₀/m₂)
  - Average period of wave zero-crossings (wind-wave scale)

- **T0M1 (Energy Period)**: T0M1 = 2π(m₋₁/m₀)
  - Period associated with wave energy transport

- **TP (Peak Period)**: TP = 2π/ω_peak
  - Period of spectral peak (dominant wave period)

### Directional Properties
- **THM (Mean Direction)**: THM = atan2(mom_y, mom_x)
  - Direction of primary wave propagation [radians]
  - Convert to degrees: THM_deg = THM × 180/π

- **THS (Directional Spread)**: THS = √(2(1 - √((mom_x² + mom_y²)/m₀²)))
  - Width of directional distribution [radians]
  - 0 = all waves from one direction, π = spread equally

### Wavelength
- **WLM (Mean Wavelength)**: WLM = 2π(mom_wn/m₀)
  - Mean wavelength of all waves in spectrum [m]

## Physical Constants

- **GRAV**: 9.81 m/s² (gravitational acceleration)
- **TPI**: 2π (used for period calculations)
- **SMALL**: 1e-20 (floor value to avoid division by zero)

## Integration with WW3

### Parameter Extraction from WW3 Code

The required parameters can be extracted from WW3 source code:

#### 1. Angular Frequency Array (omega ← SIG)
**File:** `w3gridmd.F90:3320-3380`
```fortran
! Logarithmic frequency grid
SIG(IK) = 2*PI * FR1 * XFR**(IK-1)
```

#### 2. Integration Factor (dintegral ← DDEN)
**File:** `w3gridmd.F90:3420-3450`
```fortran
! Integration factor: DDEN = DTH × DSII × SIG
DDEN(IK) = DTH * DSII(IK) * SIG(IK)
```

#### 3. Tail Factors (fte, fttr, ftwl)
**File:** `w3gridmd.F90:3444-3448`
```fortran
FTE = 0.25 * SIG(NK) * DTH * SIG(NK)
FTTR = 0.20 * DTH * SIG(NK)
FTWL = (G/6.0) / SIG(NK) * DTH * SIG(NK)
```

#### 4. Wavenumber & Group Velocity (from WAVNU1)
**File:** `w3initmd.F90:1370-1399`
```fortran
! Compute once at initialization
WN(IK,ISEA) = k  ! from dispersion relation
CG(IK,ISEA) = ∂σ/∂k  ! group velocity
```

## Performance Notes

- **Memory**: Spectrum with (36, 30) grid = ~12 KB
- **Computation**: Parameter computation = O(n_freq × n_dir)
- **Typical Time**: ~1-5 ms per spectrum on modern CPU

## References

- WW3 Source: w3iogomd.F90, w3gridmd.F90, w3initmd.F90
- Hasselmann et al. (1973): Wave spectra and energy balance
- Pierson & Moskowitz (1964): Tail frequency extension
- IAHR (1989): Wave analysis standards
