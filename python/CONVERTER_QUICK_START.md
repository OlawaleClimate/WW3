# Spectrum Converter: Quick Start Guide

## What It Does

Converts **action density spectra** (from WW3 simulations) to **energy density spectra** in various coordinate systems.

```
Action Density A(k,θ)    [wavenumber-direction space]
    ↓ (coordinate transform & Jacobian)
Energy Density E(f,θ)    [frequency-direction space]
    ├→ 2D energy spectrum: E(f,θ)
    ├→ 1D frequency spectrum: E(f)
    └→ 1D directional spectrum: E(θ)
```

## Installation

The converter is included in the `wavewatch_python` package:

```python
from wavewatch_python import action_to_energy_2d, build_ww3_grid
```

Or import individual functions:

```python
from wavewatch_python.spectrum_converter import (
    build_ww3_grid,
    action_to_energy_2d,
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d,
    get_peak_frequency,
    get_peak_direction
)
```

## Basic Usage

### Build WW3 Grid

```python
from wavewatch_python import build_ww3_grid

# Basic grid (without depth)
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)

# Access grid parameters
freq = grid['freq']          # Frequencies [Hz]
sigma = grid['sigma']        # Angular frequencies [rad/s]
dsii = grid['dsii']          # Frequency bandwidths [rad/s]
dden = grid['dden']          # DDEN conversion factors

# NEW: Also auto-compute wavenumber and group velocity by providing depth!
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36, depth=100.0)

# Now grid includes:
wn = grid['wavenumber']      # Wavenumber [1/m] (auto-computed!)
cg = grid['group_velocity']  # Group velocity [m/s] (auto-computed!)
dth = grid['dth']            # Directional bin width [rad]
fte = grid['fte']            # Tail energy factor
```

### Convert 2D Action to Energy

```python
import numpy as np
from wavewatch_python import build_ww3_grid, action_to_energy_2d
from wavewatch_python.dispersion import solve_dispersion

# Build WW3 grid
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)
omega_ww3 = grid['sigma']
depth = 100.0

# Compute group velocity
cg_ww3 = np.zeros(30)
for i in range(30):
    _, cg_ww3[i] = solve_dispersion(omega_ww3[i], depth)

# Load action density from WW3
action = np.random.rand(36, 30) * 0.001  # Your WW3 data

# Convert to energy
energy_2d, freq, dirs, grid_params = action_to_energy_2d(
    action, omega_ww3, cg_ww3
)

print(f"Energy shape: {energy_2d.shape}")
print(f"Total energy: {np.sum(energy_2d):.2f} m²")
print(f"Frequency range: {freq[0]:.4f} - {freq[-1]:.4f} Hz")
```

### Get 1D Frequency Spectrum

```python
from wavewatch_python import action_to_frequency_spectrum_1d

freq, e_freq = action_to_frequency_spectrum_1d(
    action, omega_ww3, cg_ww3
)

# Find peak
peak_freq = freq[np.argmax(e_freq)]
peak_period = 1.0 / peak_freq if peak_freq > 0 else 0.0
print(f"Peak frequency: {peak_freq:.4f} Hz ({peak_period:.2f} s)")
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

peak_freq = get_peak_frequency(action, omega_ww3, cg_ww3)
peak_dir = get_peak_direction(action, omega_ww3, cg_ww3)

print(f"Peak: {peak_freq:.4f} Hz @ {np.degrees(peak_dir):.1f}°")
```

### Normalize Spectrum

```python
from wavewatch_python import normalize_spectrum, action_to_energy_2d

energy_2d, _, _, _ = action_to_energy_2d(action, omega_ww3, cg_ww3)

# Scale to specific energy
normalized = normalize_spectrum(energy_2d, target_energy=5.0)
print(f"Total energy after normalization: {np.sum(normalized):.2f} m²")
```

### Compare with WaveSpectrum Class

```python
# For quick conversion only
energy_2d, freq, dirs, _ = action_to_energy_2d(action, omega_ww3, cg_ww3)

# For full analysis with wave parameters
from wavewatch_python import WaveSpectrum
spectrum = WaveSpectrum(action, depth, omega_ww3, dden_ww3,
                        fte_ww3, fttr_ww3, ftwl_ww3, wn_ww3, cg_ww3)
params = spectrum.compute_parameters()  # HS, T01, mean dir, etc.
```

## Function Summary

| Function | Input | Output | Use Case |
|----------|-------|--------|----------|
| `build_ww3_grid` | fr1, xfr, nk, nth | Grid dict | Build frequency/directional grids |
| `action_to_energy_2d` | 2D action | 2D energy + grids | Full spectrum analysis |
| `action_to_frequency_spectrum_1d` | 2D action | 1D freq spectrum | Peak frequency, spectral shape |
| `action_to_directional_spectrum_1d` | 2D action | 1D directional spectrum | Mean direction, spread |
| `get_peak_frequency` | 2D action | float | Quick peak frequency detection |
| `get_peak_direction` | 2D action | float | Quick peak direction detection |
| `normalize_spectrum` | 2D energy | 2D energy | Scale to target energy |

## Parameters Needed

From WW3:
- **action**: 2D array (ndir, nfreq) - your action density data
- **omega**: Angular frequencies from SIG array [rad/s]
- **group_velocity**: From CG computation (or use `solve_dispersion`)
- **depth**: Water depth for dispersion relation

Optional:
- **directions**: Direction grid in radians (auto-generated if not provided)
- **ddir**: Directional bin width (auto-computed if not provided)
- **dsii**: Frequency bandwidths (auto-computed if not provided)

## Physical Meaning

The conversion transforms action density from wavenumber-direction space to energy density in frequency-direction space:

```
E(f,θ) = A(k,θ) × σ × (2π / CG)
       = A(k,θ) × SIG × (2π / CG)
```

**Key Relationships:**
- **E(f,θ)**: Energy density in frequency-direction space [m²/Hz/rad]
- **A(k,θ)**: Action density in wavenumber-direction space [m²·s·rad⁻¹]
- **Action-Energy**: A(k,θ) = F(k,θ) / σ (action = energy / intrinsic frequency)
- **σ = SIG**: Intrinsic (angular) frequency [rad/s] from WW3
- **∂k/∂f = 2π/CG**: Jacobian of the (k,θ) → (f,θ) transformation
- **CG**: Group velocity [m/s]

**Integration to 1D Spectra:**
- **1D Frequency**: E(f) = ∑_θ E(f,θ) × DTH (includes directional bin width)
- **1D Directional**: E(θ) = ∑_f E(f,θ) × DSII[f] (includes frequency bandwidths)

## Related Documentation

For detailed information, see:
- `SPECTRUM_CONVERTER_UPDATE.md` - Latest implementation changes
- `python/wavewatch_python/SPECTRUM_CONVERTER_GUIDE.md` - Complete API reference
- `test_spectrum_converter.py` - Working examples and test suite

## Performance

- Single spectrum: < 1 ms
- Batch (1000 spectra): ~1 second
- Memory: ~8.6 KB per (36,30) spectrum

## Common Errors

**Error**: `ValueError: omega length != nfreq`
- **Fix**: Ensure omega array matches action frequency bins

**Error**: `ValueError: group_velocity length != nfreq`
- **Fix**: Ensure CG array matches action frequency bins

**Error**: Peak detection returns wrong values
- **Fix**: Ensure omega and group_velocity are correctly computed for your depth

**Warning**: Results don't match WaveSpectrum
- **Note**: Converter excludes tail extension by default. Use WaveSpectrum for full analysis with tail.

## Next Steps

1. Run the test suite: `python test_spectrum_converter.py`
2. Read full implementation details: `SPECTRUM_CONVERTER_UPDATE.md`
3. Check API reference: `python/wavewatch_python/SPECTRUM_CONVERTER_GUIDE.md`
4. Explore test examples: `test_spectrum_converter.py`
