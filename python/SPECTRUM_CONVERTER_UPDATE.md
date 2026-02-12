# Spectrum Converter: Update to Reference Implementation

## Summary

The `spectrum_converter.py` module has been updated to align with the reference implementation from the `claude/wave-spectrum-analysis-TgpfA` branch. The updates ensure exact compliance with WW3 Fortran algorithms (w3iogomd.F90 and w3gridmd.F90) and correct handling of spectral transformations.

## Key Changes

### 1. Added `build_ww3_grid()` Function

**Purpose**: Construct WW3 spectral grid arrays exactly as in w3gridmd.F90

**Features**:
- Builds frequency and directional grids with proper geometric spacing
- Computes frequency bandwidths (DSII) with special cases for first and last bins:
  ```python
  dsii[0] = 0.5 * sigma[0] * (xfr - 1.0)           # first bin
  dsii[-1] = 0.5 * sigma[-1] * (xfr - 1.0) / xfr   # last bin
  ```
- Calculates DDEN factor: `DDEN = DTH × DSII × SIG`
- Computes tail energy factor: `FTE = 0.25 × SIG[-1] × DTH × SIG[-1]`

**Returns** a dict with: `freq`, `sigma`, `theta`, `dth`, `dsii`, `dden`, `fte`

### 2. Fixed `action_to_energy_2d()` Function

**Changes**:
- Now returns grid parameters (dth, dsii, dden) for use in 1D integration
- Clarified documentation that this is a **pure spectral conversion** with NO bin-width factors
- Bin widths are applied only when integrating to 1D spectra

**Physics**:
```
E(f,θ) = A(k,θ) × σ × (2π / CG)
```

**Key Point**: The 2D energy conversion itself does not include DTH or DSII factors.

### 3. Updated `action_to_frequency_spectrum_1d()` Function

**Major Fix**: Now properly includes DTH (directional bin width)

**Formula**:
```python
E(f) = ∑_θ E(f,θ) × DTH = ∑_θ [A(θ,f) × SIG × (2π/CG)] × DTH
```

**Why**: When integrating over directions, we must multiply by the directional bin width DTH to account for discrete integration over the directional domain.

### 4. Updated `action_to_directional_spectrum_1d()` Function

**Major Fix**: Now properly includes DSII (frequency bandwidth)

**Formula**:
```python
E(θ) = ∑_f [A(θ,f) × SIG × (2π/CG) × DSII]
```

**Why**: When integrating over frequencies, we must multiply by the frequency bandwidth DSII (which varies for each frequency bin) to account for discrete integration over the frequency domain.

### 5. Updated Peak Detection Functions

**Changes to signatures**:

Old:
```python
peak_freq = get_peak_frequency(action_2d, frequencies)
peak_dir = get_peak_direction(action_2d, directions)
```

New:
```python
peak_freq = get_peak_frequency(action_2d, omega, group_velocity,
                                directions=None, ddir=None)
peak_dir = get_peak_direction(action_2d, omega, group_velocity,
                               directions=None, ddir=None)
```

**Reason**: These functions now use the corrected 1D spectrum computation methods internally, requiring access to omega and group_velocity for proper conversion.

## Physics Explanation

### Action vs. Energy Coordinates

The key physics insight is that WW3 stores **action density** in **wavenumber-direction space** A(k,θ), but users need **energy density** in **frequency-direction space** E(f,θ).

```
Input:  A(k,θ)  [m²·s·rad⁻¹] - Action in wavenumber space
         ↓
         Coordinate transformation using Jacobian
         ↓
Output: E(f,θ)  [m²/Hz/rad] - Energy in frequency space
```

### Conversion Formula

The complete transformation combines:

1. **Action-to-Energy relationship**: `F(k,θ) = A(k,θ) × σ`
   - σ = SIG = intrinsic (angular) frequency [rad/s]

2. **Coordinate Jacobian**: `∂k/∂f = 2π / CG`
   - From dispersion relation: σ² = g·k·tanh(k·h)
   - CG = dσ/dk = group velocity [m/s]

3. **Result**: `E(f,θ) = A(k,θ) × SIG × (2π / CG)`

### Integration to 1D Spectra

When integrating the 2D energy to get 1D spectra, bin widths must be included:

**Frequency spectrum** (sum over directions):
```
E(f) = ∑_θ E(f,θ) × DTH
```
- Must include DTH because we're summing discrete directional samples
- DTH = 2π / NTH [radians]

**Directional spectrum** (sum over frequencies):
```
E(θ) = ∑_f E(f,θ) × DSII[f]
```
- Must include DSII because frequency bins have different widths (log-spaced grid)
- DSII varies per frequency bin

## Test Results

All tests pass successfully:

```
✓ Test 1: Action to Energy 2D Conversion
✓ Test 2: Action to 1D Frequency Spectrum
✓ Test 3: Action to 1D Directional Spectrum
✓ Test 4: Peak Frequency and Direction Detection
✓ Test 5: Spectrum Normalization
✓ Test 6: Verification Against WaveSpectrum Class
```

## Example Usage

### Basic Conversion

```python
from wavewatch_python import build_ww3_grid, action_to_energy_2d

# Build grid
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)

# Convert 2D action to energy
energy_2d, freq, dirs, grid_params = action_to_energy_2d(
    action, grid['sigma'], group_velocity
)
```

### Get 1D Spectra

```python
from wavewatch_python import (
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d
)

# Frequency spectrum
freq, e_freq = action_to_frequency_spectrum_1d(
    action, omega, group_velocity
)

# Directional spectrum
dirs, e_dir = action_to_directional_spectrum_1d(
    action, omega, group_velocity
)
```

### Detect Peaks

```python
from wavewatch_python import get_peak_frequency, get_peak_direction

peak_freq = get_peak_frequency(action, omega, group_velocity)
peak_dir = get_peak_direction(action, omega, group_velocity)
```

## References

- WW3 Fortran Source: w3iogomd.F90 (lines 1484-2113), w3gridmd.F90
- Physics: Hasselmann et al. (1973), Whitham (1974), IAHR standards
- Reference Implementation: `claude/wave-spectrum-analysis-TgpfA:model/tools/ww3_spectral_analysis.py`

## Backward Compatibility

This update is **NOT backward compatible** with previous versions due to:
1. Changed return format of `action_to_energy_2d()` (now returns grid_params dict)
2. Changed function signatures for `get_peak_frequency()` and `get_peak_direction()`

Users should migrate to the new API as shown in the examples above.

## Verification

The implementation has been verified to:
- ✓ Match the physics of the reference WW3 implementation
- ✓ Follow w3iogomd.F90 and w3gridmd.F90 algorithms exactly
- ✓ Handle DSII computation with proper special cases
- ✓ Properly include DTH and DSII factors in 1D integrations
- ✓ Pass all test cases
- ✓ Align with the reference implementation from the other branch
