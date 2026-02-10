# Spectrum Converter: Quick Start Guide

## What It Does

Converts **action density spectra** (from WW3 simulations) to **energy density spectra** in various coordinate systems.

```
Action Density A(f,θ)
    ↓ (×DDEN/CG)
Energy Density E(f,θ)
    ├→ 2D energy spectrum
    ├→ 1D frequency spectrum
    └→ 1D directional spectrum
```

## Installation

The converter is included in the `wavewatch_python` package:

```python
from wavewatch_python import action_to_energy_2d
```

Or import individual functions:

```python
from wavewatch_python.spectrum_converter import (
    action_to_energy_2d,
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d,
    get_peak_frequency,
    get_peak_direction
)
```

## Basic Usage

### Convert 2D Action to Energy

```python
import numpy as np
from wavewatch_python import action_to_energy_2d
from wavewatch_python.dispersion import solve_dispersion

# WW3 parameters
omega_ww3 = 2*np.pi * 0.04 * 1.1**np.arange(30)  # Frequency grid
depth = 100.0

# Compute group velocity
cg_ww3 = np.zeros(30)
for i in range(30):
    _, cg_ww3[i] = solve_dispersion(omega_ww3[i], depth)

# Load action density from WW3
action = np.random.rand(36, 30) * 0.001  # Your WW3 data

# Convert to energy
energy_2d, freq, dirs, dden = action_to_energy_2d(
    action, omega_ww3, cg_ww3
)

print(f"Energy shape: {energy_2d.shape}")
print(f"Total energy: {np.sum(energy_2d):.2f} m²")
```

### Get 1D Frequency Spectrum

```python
from wavewatch_python import action_to_frequency_spectrum_1d

freq, e_freq = action_to_frequency_spectrum_1d(
    action, omega_ww3, cg_ww3
)

# Plot or analyze
peak_freq = freq[np.argmax(e_freq)]
print(f"Peak frequency: {peak_freq:.4f} Hz")
print(f"Peak period: {1/peak_freq:.2f} s")
```

### Get 1D Directional Spectrum

```python
from wavewatch_python import action_to_directional_spectrum_1d

dirs, e_dir = action_to_directional_spectrum_1d(
    action, omega_ww3, cg_ww3
)

# Find mean direction
mean_dir = np.arctan2(np.sum(e_dir * np.sin(dirs)),
                      np.sum(e_dir * np.cos(dirs)))
print(f"Mean direction: {np.degrees(mean_dir):.1f}°")
```

## Common Tasks

### Extract Peak Properties

```python
from wavewatch_python import get_peak_frequency, get_peak_direction

peak_freq = get_peak_frequency(action, freq)
peak_dir = get_peak_direction(action, dirs)

print(f"Peak: {peak_freq:.4f} Hz @ {np.degrees(peak_dir):.1f}°")
```

### Normalize Spectrum

```python
from wavewatch_python import normalize_spectrum, action_to_energy_2d

energy_2d, _, _, _ = action_to_energy_2d(action, omega_ww3, cg_ww3)

# Scale to specific energy
normalized = normalize_spectrum(energy_2d, target_energy=5.0)
print(f"Total energy: {np.sum(normalized):.2f} m²")
```

### Compare with WaveSpectrum Class

```python
# For quick conversion only
energy, freq, dirs, _ = action_to_energy_2d(action, omega_ww3, cg_ww3)

# For full analysis with wave parameters
from wavewatch_python import WaveSpectrum
spectrum = WaveSpectrum(action, depth, omega_ww3, dden_ww3,
                        fte_ww3, fttr_ww3, ftwl_ww3, wn_ww3, cg_ww3)
params = spectrum.compute_parameters()  # HS, T01, mean dir, etc.
```

## Function Summary

| Function | Input | Output | Use Case |
|----------|-------|--------|----------|
| `action_to_energy_2d` | 2D action | 2D energy | Full spectrum analysis |
| `action_to_frequency_spectrum_1d` | 2D action | 1D freq | Peak frequency, spectral shape |
| `action_to_directional_spectrum_1d` | 2D action | 1D dir | Mean direction, spread |
| `get_peak_frequency` | 2D action | float | Quick peak detection |
| `get_peak_direction` | 2D action | float | Quick peak detection |
| `normalize_spectrum` | 2D energy | 2D energy | Scale to target energy |

## Parameters Needed

From WW3:
- **action**: 2D array (ndir, nfreq) - your spectral data
- **omega**: Angular frequencies from SIG array
- **group_velocity**: From CG computation (or use `solve_dispersion`)
- **depth**: Water depth for dispersion relation

Optional:
- **dintegral**: Pre-computed DDEN factors (auto-computed if not provided)
- **directions**: Direction grid (auto-generated if not provided)

## Physical Meaning

The conversion uses:
```
E(f,θ) = A(f,θ) × (DDEN / CG)
```

Where:
- **E**: Energy density [m²/Hz/rad]
- **A**: Action density [m²·s·rad⁻¹]
- **DDEN**: Integration factor = DTH × DSII × ω
- **CG**: Group velocity [m/s]

This accounts for:
1. Directional binning (DTH)
2. Frequency binning (DSII)
3. Frequency scaling (ω)
4. Wave frame vs. lab frame (CG)

## Full Documentation

For detailed API reference, examples, and mathematical background, see:
- `SPECTRUM_CONVERTER_GUIDE.md` - Complete reference
- `test_spectrum_converter.py` - Test examples

## Performance

- Single spectrum: < 1 ms
- Batch (1000 spectra): ~1 second
- Memory: ~8.6 KB per (36,30) spectrum

## Common Errors

**Error**: `ValueError: omega length != nfreq`
- **Fix**: Ensure omega array matches action frequency bins

**Error**: `ValueError: group_velocity length != nfreq`
- **Fix**: Ensure CG array matches action frequency bins

**Warning**: Results don't match WaveSpectrum
- **Note**: Converter excludes tail extension by default. Use WaveSpectrum for full analysis.

## Next Steps

1. Try the test: `python test_spectrum_converter.py`
2. Read full guide: `SPECTRUM_CONVERTER_GUIDE.md`
3. Check examples in test suite or WW3_INTEGRATION.md
