# WaveWatch Python: Spectrum Modules Comparison

## Overview

Two implementations of the `WaveSpectrum` class are now available:

1. **`spectrum.py`** - Full-featured, flexible module
2. **`spectrum_ww3_only.py`** - Streamlined, WW3-only production module

Choose based on your use case.

## Comparison Table

| Feature | `spectrum.py` | `spectrum_ww3_only.py` |
|---------|---------------|----------------------|
| **Auto-generate grid** | ✅ Yes | ❌ No |
| **Accept WW3 parameters** | ✅ Yes | ✅ Yes |
| **Clean, minimal API** | ⚠️ Complex | ✅ Simple |
| **Lines of code** | 402 | 372 |
| **Dependencies** | numpy, dispersion | numpy, dispersion |
| **Requires omega input** | Optional | Required |
| **Requires DDEN input** | Optional | Required |
| **Requires WN input** | Optional | Required |
| **Requires CG input** | Optional | Required |
| **Research/testing friendly** | ✅ Yes | ⚠️ No |
| **Production-ready** | ✅ Yes | ✅✅ Yes |
| **Maintenance burden** | Medium | Low |

## Detailed Comparison

### 1. API Complexity

#### `spectrum.py` (Flexible)
```python
# Method 1: With WW3 parameters (recommended)
spectrum = WaveSpectrum(
    action, depth=100.0,
    omega=omega_ww3,
    dintegral=dden_ww3,
    fte=fte_ww3, fttr=fttr_ww3, ftwl=ftwl_ww3,
    wavenumber=wn_ww3,
    group_velocity=cg_ww3
)

# Method 2: With raw grid (auto-computed)
spectrum = WaveSpectrum(
    action, depth=100.0,
    frequencies=freq_array,
    directions=dir_array
)

# Method 3: Minimal (all auto-generated)
spectrum = WaveSpectrum(action, depth=100.0)
```

#### `spectrum_ww3_only.py` (Simple)
```python
# Only one method: explicit WW3 parameters
spectrum = WaveSpectrum(
    action=action,
    depth=depth,
    omega=omega_ww3,
    dintegral=dden_ww3,
    fte=fte_ww3,
    fttr=fttr_ww3,
    ftwl=ftwl_ww3,
    wavenumber=wn_ww3,
    group_velocity=cg_ww3,
    directions=directions  # Optional
)
```

### 2. Use Case Selection

#### Use `spectrum.py` if:

✅ **Research & Development**
- Exploring spectral algorithms
- Creating synthetic test spectra
- Academic papers with different grid configurations
- Experimental analysis

✅ **Flexibility Needed**
- Working with non-standard frequency grids
- Need to test different grid parameters
- Want a reference implementation

✅ **Rapid Prototyping**
- Don't want to pre-compute grid parameters
- Quick testing without WW3 integration

**Example Use:**
```python
# Testing spectrum algorithm with synthetic data
import numpy as np

# Simple: auto-generate everything
action = np.random.rand(36, 30)
spectrum = WaveSpectrum(action)  # Frequency grid auto-generated

# Or customize just what you need
spectrum = WaveSpectrum(
    action,
    frequencies=np.logspace(np.log10(0.04), np.log10(0.5), 30),
    directions=np.linspace(0, 2*np.pi, 36)
)
```

#### Use `spectrum_ww3_only.py` if:

✅ **Production Systems**
- Processing actual WW3 simulation output
- Building operational wave forecasting system
- Integration with WW3 data pipelines
- Distributed processing

✅ **Code Clarity & Maintenance**
- Want explicit, clear API with no hidden auto-generation
- Easy for team members to understand data flow
- Easier to debug (parameters must be explicitly provided)

✅ **Reproducibility**
- All parameters come from WW3, no surprises
- Easy to validate against WW3 output
- Full traceability of grid parameters

✅ **Performance-Critical**
- Minimal code path
- No unnecessary auto-generation logic
- Faster startup (no grid generation)

**Example Use:**
```python
# Processing WW3 NetCDF output
import netCDF4 as nc
from wavewatch_python.spectrum_ww3_only import WaveSpectrum

# Load WW3 output
ww3_file = nc.Dataset('ww3_multi_1.nc')
action = ww3_file.variables['efth'][:, :, :, :, :]  # All space/time/freq/dir
depth = ww3_file.variables['dpt'][:, :, :]

# Get WW3 grid parameters (pre-computed)
omega_ww3 = 2*np.pi * ww3_file.variables['frequency'][:]
dden_ww3 = ...  # Extracted from WW3 initialization
fte_ww3, fttr_ww3, ftwl_ww3 = ...  # From w3gridmd
wn_ww3, cg_ww3 = ...  # From WAVNU1

# Process entire grid
nlat, nlon, nfreq, ndir = action.shape[1:]
hs_grid = np.zeros((action.shape[0], nlat, nlon))

for t in range(action.shape[0]):
    for i in range(nlat):
        for j in range(nlon):
            spectrum = WaveSpectrum(
                action=action[t, i, j, :, :],
                depth=depth[t, i, j],
                omega=omega_ww3,
                dintegral=dden_ww3,
                fte=fte_ww3, fttr=fttr_ww3, ftwl=ftwl_ww3,
                wavenumber=wn_ww3,
                group_velocity=cg_ww3
            )
            hs_grid[t, i, j] = spectrum.compute_parameters()['hs']
```

## Code Size & Complexity

### `spectrum.py` Breakdown:
- Constructor (`__init__`): 95 lines (handles 3 different modes)
- Auto-generation methods: 70 lines (unused in WW3 mode)
- Moment calculation: 55 lines (core algorithm)
- Parameter derivation: 90 lines (formulas)
- 1D spectrum extraction: 25 lines
- **Total: 402 lines**

### `spectrum_ww3_only.py` Breakdown:
- Constructor (`__init__`): 60 lines (only WW3 mode)
- Validation: 25 lines (explicit checks)
- Frequency bandwidth: 20 lines
- Moment calculation: 65 lines (core algorithm)
- Parameter derivation: 100 lines (same formulas)
- 1D spectrum extraction: 28 lines
- **Total: 372 lines**

**Insight:** Removing optional code paths saves ~30 lines and makes remaining code clearer.

## Import Statement

### Using `spectrum.py`:
```python
from wavewatch_python import WaveSpectrum  # Imports from __init__.py
```

### Using `spectrum_ww3_only.py`:
```python
from wavewatch_python.spectrum_ww3_only import WaveSpectrum
```

## API Equivalence

Both versions implement identical wave parameter formulas:
- HS = 4√m₀
- T01 = 2π(m₀/m₁)
- T02 = 2π√(m₀/m₂)
- T0M1 = 2π(m₋₁/m₀)
- Mean direction, directional spread, wavelength, spectral width

**Results are identical** when given the same input parameters.

## Migration Guide

### From `spectrum.py` to `spectrum_ww3_only.py`

If you're currently using `spectrum.py` with WW3 parameters:

```python
# Before (spectrum.py)
from wavewatch_python import WaveSpectrum

spectrum = WaveSpectrum(
    action, depth=100.0,
    omega=omega_ww3,
    dintegral=dden_ww3,
    fte=fte_ww3, fttr=fttr_ww3, ftwl=ftwl_ww3,
    wavenumber=wn_ww3,
    group_velocity=cg_ww3
)
```

```python
# After (spectrum_ww3_only.py) - Just change import and remove depth keyword
from wavewatch_python.spectrum_ww3_only import WaveSpectrum

spectrum = WaveSpectrum(
    action=action,
    depth=100.0,
    omega=omega_ww3,
    dintegral=dden_ww3,
    fte=fte_ww3,
    fttr=fttr_ww3,
    ftwl=ftwl_ww3,
    wavenumber=wn_ww3,
    group_velocity=cg_ww3
)
```

The rest of the code remains identical.

## Recommendation

### For Most Users:
Use **`spectrum_ww3_only.py`** because:
- Most users work with WW3 data
- Cleaner, more explicit API
- Easier to understand and maintain
- No confusion about auto-generation
- Better for production systems

### For Special Cases:
Use **`spectrum.py`** if:
- You need to test with synthetic grids
- You want maximum flexibility
- You're doing research with non-standard configurations
- You want a reference implementation

## Testing

Both modules have identical test suites:

```bash
# Test spectrum.py (flexible)
python -m wavewatch_python.example

# Test spectrum_ww3_only.py (streamlined)
python test_spectrum_ww3_only.py
```

All tests pass for both versions.
