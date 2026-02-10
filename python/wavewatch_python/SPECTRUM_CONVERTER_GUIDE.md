# Wave Spectrum Converter: Action Density to Energy Spectrum

## Overview

The spectrum converter module provides utility functions to transform **action density spectra** (from WW3 simulations) into **energy density spectra** in various coordinate systems (frequency-direction, 1D frequency, 1D directional).

This is essential for post-processing WW3 output and analyzing wave energy distributions.

## Key Concepts

### Action Density vs. Energy Density

| Property | Action Density | Energy Density |
|----------|---|---|
| **Definition** | A(f,θ) = E(f,θ) / CG | E(f,θ) in energy domain |
| **Units** | m²·s·rad⁻¹ | m²/Hz/rad |
| **Reference frame** | Moving with waves | Fixed (lab) frame |
| **Reason for difference** | Group velocity causes Doppler shift | Standard wave measurements |
| **Conversion** | E = A × (DDEN / CG) | A = E / (DDEN / CG) |

### Transformation Process

```
Action Density [m²·s·rad⁻¹]
         ↓
         │ Apply conversion factor: DDEN / CG
         ↓
Energy Density [m²/Hz/rad]
         ↓
    ┌────┴─────┐
    ↓          ↓
1D Freq      1D Dir
Spectrum     Spectrum
```

## Function Reference

### Main Conversion Functions

#### `action_to_energy_2d()`

Converts 2D action density to 2D energy density in frequency-direction space.

```python
energy_2d, freq, dirs, dintegral_used = action_to_energy_2d(
    action,           # (ndir, nfreq) action spectrum
    omega,            # Angular frequencies [rad/s]
    group_velocity,   # Group velocity [m/s]
    dintegral=None,   # Optional: DDEN factors
    directions=None,  # Optional: direction grid
    ddir=None         # Optional: directional spacing
)
```

**Parameters:**
- `action` (ndarray): Action density (ndir, nfreq) [m²·s·rad⁻¹]
- `omega` (ndarray): Angular frequencies from SIG [rad/s]
- `group_velocity` (ndarray): Group velocity from CG [m/s]
- `dintegral` (ndarray, optional): DDEN = DTH × DSII × SIG
- `directions` (ndarray, optional): Direction grid [radians]
- `ddir` (float, optional): Directional spacing [radians]

**Returns:**
- `energy_2d`: Energy density (ndir, nfreq) [m²/Hz/rad]
- `freq`: Frequency array [Hz]
- `dirs`: Direction array [radians]
- `dintegral_used`: Integration factors actually used

**Example:**
```python
from wavewatch_python.spectrum_converter import action_to_energy_2d

# Convert WW3 action density to energy
energy, freq, dirs, dden = action_to_energy_2d(
    action=ww3_action,
    omega=omega_ww3,
    group_velocity=cg_ww3,
    dintegral=dden_ww3  # Pre-computed DDEN from WW3
)

print(f"Energy shape: {energy.shape}")
print(f"Total energy: {np.sum(energy):.2f} m²")
```

---

#### `action_to_frequency_spectrum_1d()`

Converts 2D action density to 1D frequency spectrum by integrating over all directions.

```python
freq, e_freq = action_to_frequency_spectrum_1d(
    action,           # (ndir, nfreq) action spectrum
    omega,            # Angular frequencies [rad/s]
    group_velocity,   # Group velocity [m/s]
    dintegral=None,   # Optional: DDEN factors
    directions=None,  # Optional: direction grid
    ddir=None         # Optional: directional spacing
)
```

**Returns:**
- `freq`: Frequency array [Hz]
- `e_freq`: 1D energy spectrum [m²/Hz]

**Physical Meaning:**
- Sums all directional components at each frequency
- Result is traditional "energy vs. frequency" spectrum
- Used for: peak frequency, mean period, spectral shape

**Example:**
```python
freq, e_freq = action_to_frequency_spectrum_1d(
    action_ww3, omega_ww3, cg_ww3
)

# Find peak frequency
peak_freq = freq[np.argmax(e_freq)]
print(f"Peak frequency: {peak_freq:.4f} Hz")
print(f"Peak period: {1/peak_freq:.2f} s")
```

---

#### `action_to_directional_spectrum_1d()`

Converts 2D action density to 1D directional spectrum by integrating over all frequencies.

```python
dirs, e_dir = action_to_directional_spectrum_1d(
    action,           # (ndir, nfreq) action spectrum
    omega,            # Angular frequencies [rad/s]
    group_velocity,   # Group velocity [m/s]
    dintegral=None,   # Optional: DDEN factors
    directions=None,  # Optional: direction grid
    ddir=None         # Optional: directional spacing
)
```

**Returns:**
- `dirs`: Direction array [radians]
- `e_dir`: 1D directional spectrum [m²/rad]

**Physical Meaning:**
- Sums all frequency components at each direction
- Result is "energy vs. direction" distribution
- Used for: mean direction, directional spread

**Example:**
```python
dirs, e_dir = action_to_directional_spectrum_1d(
    action_ww3, omega_ww3, cg_ww3
)

# Mean direction
mean_dir = np.sum(e_dir * np.cos(dirs)) / np.sum(e_dir)
print(f"Mean direction: {np.degrees(mean_dir):.1f}°")
```

---

### Utility Functions

#### `get_peak_frequency()`

Extract peak frequency from action spectrum.

```python
peak_freq = get_peak_frequency(action_2d, frequencies)
```

#### `get_peak_direction()`

Extract peak direction from action spectrum.

```python
peak_dir = get_peak_direction(action_2d, directions)
```

#### `normalize_spectrum()`

Normalize 2D energy spectrum to specific total energy.

```python
normalized = normalize_spectrum(energy_2d, target_energy=10.0)
```

## Usage Examples

### Example 1: Basic Conversion

```python
import numpy as np
from wavewatch_python.spectrum_converter import action_to_energy_2d
from wavewatch_python.dispersion import solve_dispersion

# WW3 parameters
nfreq, ndir = 30, 36
depth = 100.0
omega_ww3 = 2*np.pi * 0.04 * 1.1**np.arange(nfreq)

# Compute group velocity
cg_ww3 = np.zeros(nfreq)
for ik in range(nfreq):
    _, cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

# Load action density from WW3
action_ww3 = load_from_ww3_output(...)

# Convert to energy spectrum
energy_2d, freq, dirs, _ = action_to_energy_2d(
    action_ww3, omega_ww3, cg_ww3
)

print(f"Energy spectrum shape: {energy_2d.shape}")
print(f"Total energy: {np.sum(energy_2d):.2f} m²")
```

### Example 2: Extract and Analyze 1D Spectra

```python
from wavewatch_python.spectrum_converter import (
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d,
    get_peak_frequency,
    get_peak_direction
)

# Get 1D frequency spectrum
freq, e_freq = action_to_frequency_spectrum_1d(
    action_ww3, omega_ww3, cg_ww3
)

# Get 1D directional spectrum
dirs, e_dir = action_to_directional_spectrum_1d(
    action_ww3, omega_ww3, cg_ww3
)

# Extract peaks
peak_freq = get_peak_frequency(action_ww3, freq)
peak_dir = get_peak_direction(action_ww3, dirs)

print(f"Peak frequency: {peak_freq:.4f} Hz")
print(f"Peak direction: {np.degrees(peak_dir):.1f}°")
```

### Example 3: Multi-Point Processing

```python
# Process WW3 grid output (time, lat, lon, freq, dir)
action_data = ww3_output['efth']  # Shape: (nt, nlat, nlon, nfreq, ndir)

# Reshape to (ndir, nfreq)
nt, nlat, nlon = action_data.shape[:3]
hs_grid = np.zeros((nt, nlat, nlon))
peak_freq_grid = np.zeros((nt, nlat, nlon))

for t in range(nt):
    for i in range(nlat):
        for j in range(nlon):
            # Convert action to energy
            energy_2d, _, _, _ = action_to_energy_2d(
                action_data[t, i, j, :, :],  # (nfreq, ndir)
                omega_ww3, cg_ww3
            )

            # Get 1D frequency spectrum
            freq, e_freq = action_to_frequency_spectrum_1d(
                action_data[t, i, j, :, :],
                omega_ww3, cg_ww3
            )

            # Extract properties
            hs_grid[t, i, j] = 4 * np.sqrt(np.sum(energy_2d))
            peak_freq_grid[t, i, j] = freq[np.argmax(e_freq)]
```

### Example 4: Spectral Analysis

```python
# Compare spectra at different times
for t in [0, 6, 12]:
    freq, e_freq = action_to_frequency_spectrum_1d(
        action_ww3[t], omega_ww3, cg_ww3
    )

    # Analyze spectrum
    total_energy = np.sum(e_freq) * (freq[1] - freq[0])
    peak_freq = freq[np.argmax(e_freq)]
    peak_energy = e_freq.max()

    print(f"Time {t}:")
    print(f"  Total energy: {total_energy:.2f} m²")
    print(f"  Peak frequency: {peak_freq:.4f} Hz")
    print(f"  Peak energy: {peak_energy:.6e} m²/Hz")
```

## Mathematical Background

### Conversion Formula

The fundamental relationship between action and energy is:

```
E(f,θ) = A(f,θ) × (DDEN / CG)
```

Where:
- **E(f,θ)**: Energy density [m²/Hz/rad]
- **A(f,θ)**: Action density [m²·s·rad⁻¹]
- **DDEN**: Integration factor = DTH × DSII × ω
- **CG**: Group velocity [m/s]

### Why This Conversion?

1. **WW3 uses action density** because it's conserved during wave propagation across depth gradients
2. **Users need energy density** for standard wave statistics (HS, peak period, etc.)
3. **Group velocity compensation** accounts for the wave frame vs. lab frame difference

### Integration Factors

**DDEN = DTH × DSII × ω**

Where:
- **DTH**: Directional bin width = 2π/NTH [radians]
- **DSII**: Frequency bin width (variable for logarithmic grid)
- **ω**: Angular frequency [rad/s]

These are pre-computed in WW3's w3gridmd.F90 and provided as the DDEN array.

## Performance Notes

- **Memory**: (36, 30) spectrum ≈ 8.6 KB (float64)
- **Speed**: Single spectrum conversion < 1 ms on modern CPU
- **Batch processing**: ~1000 spectra/second with optimized code

## Integration with WaveSpectrum Class

The converter is standalone but compatible with `WaveSpectrum`:

```python
from wavewatch_python.spectrum_converter import action_to_energy_2d
from wavewatch_python.spectrum_ww3_only import WaveSpectrum

# Method 1: Use converter directly (faster for just energy)
energy, freq, dirs, _ = action_to_energy_2d(
    action, omega_ww3, cg_ww3
)

# Method 2: Use WaveSpectrum for full analysis (includes tail)
spectrum = WaveSpectrum(action, depth, omega_ww3, dden_ww3,
                        fte_ww3, fttr_ww3, ftwl_ww3,
                        wn_ww3, cg_ww3)
params = spectrum.compute_parameters()
```

**Key Differences:**
- **Converter**: Fast, pure energy conversion, no tail extension
- **WaveSpectrum**: Full analysis, includes Pierson-Moskowitz tail, computes all parameters

## References

- WW3 Source: w3iogomd.F90, w3gridmd.F90
- Action density concept: Hasselmann et al. (1973)
- Group velocity: Whitham (1974)
- Wave statistics: IAHR standards
