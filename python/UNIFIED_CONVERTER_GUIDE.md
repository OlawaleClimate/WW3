# Unified Spectrum Converter - Single API for All Use Cases

## The Problem You Identified

Why have multiple converters for WW3 and non-WW3 data when the core computation is the same?

**Key Insight**:
- `omega` (angular frequency) can be computed from grid parameters (fr1, xfr, nk)
- `group_velocity` (cg) can be computed from omega and depth using the dispersion relation
- Everything else follows from these two parameters

## The Solution: `SpectrumConverter`

A **single, unified class** that intelligently handles:
- Generic data (auto-computes what's needed)
- WW3 data (uses pre-computed parameters for efficiency)
- Mixed scenarios (accepts any combination)

## Three Usage Patterns

All three give **identical results** but with different input/efficiency tradeoffs:

### Pattern 1: Generic Data (Auto-compute everything)

```python
from wavewatch_python import SpectrumConverter

converter = SpectrumConverter(
    action,
    depth=100.0,
    fr1=0.04,      # First frequency [Hz]
    xfr=1.1,       # Frequency ratio
    nk=30,         # Number of frequencies
    nth=36         # Number of directions
)

# All parameters auto-computed from grid specs + depth
energy_2d, freq, dirs = converter.to_energy_2d()
hs = converter.compute_wave_parameters()['hs']
```

**Use when**: You have any data source and want simplicity

---

### Pattern 2: WW3 Data (Pre-computed omega, cg)

```python
from wavewatch_python import SpectrumConverter, build_ww3_grid
from wavewatch_python.dispersion import solve_dispersion
import numpy as np

# Build grid and compute parameters (from WW3 simulation)
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)
omega = grid['sigma']
cg = np.array([solve_dispersion(omega[i], 100.0)[1] for i in range(len(omega))])

converter = SpectrumConverter(
    action,
    omega=omega,
    group_velocity=cg,
    depth=100.0
)

# Uses pre-computed omega and cg (faster, cleaner)
energy_2d, freq, dirs = converter.to_energy_2d()
hs = converter.compute_wave_parameters()['hs']
```

**Use when**: You have WW3 omega and cg pre-computed (most common)

---

### Pattern 3: WW3 Optimized (All parameters pre-computed via explicit parameters)

```python
from wavewatch_python import SpectrumConverter, build_ww3_grid
from wavewatch_python.dispersion import solve_dispersion
import numpy as np

# Pre-compute ALL parameters (from WW3)
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)
omega = grid['sigma']
wn = np.array([solve_dispersion(omega[i], 100.0)[0] for i in range(len(omega))])
cg = np.array([solve_dispersion(omega[i], 100.0)[1] for i in range(len(omega))])

converter = SpectrumConverter(
    action,
    omega=omega,
    group_velocity=cg,
    dintegral=grid['dden'],  # Pre-computed conversion factors
    wavenumber=wn,           # Pre-computed wavenumber
    depth=100.0
)

# Most efficient - uses all pre-computed parameters
energy_2d, freq, dirs = converter.to_energy_2d()
hs = converter.compute_wave_parameters(fte=grid['fte'])['hs']
```

**Use when**: You want maximum efficiency (skip all redundant computations)

---

### Pattern 4: WW3 Optimized (All parameters via dictionary - Cleaner API)

```python
from wavewatch_python import SpectrumConverter, build_ww3_grid
from wavewatch_python.dispersion import solve_dispersion
import numpy as np

# Pre-compute ALL parameters
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)
omega = grid['sigma']
wn = np.array([solve_dispersion(omega[i], 100.0)[0] for i in range(len(omega))])
cg = np.array([solve_dispersion(omega[i], 100.0)[1] for i in range(len(omega))])

# Your action density from WW3
action = load_action_spectrum('data.nc')

# Pass all pre-computed WW3 parameters via dict (cleaner!)
converter = SpectrumConverter(
    action,
    omega=omega,
    group_velocity=cg,
    depth=100.0,
    ww3_params={
        'dden': grid['dden'],          # Conversion factors
        'wn': wn,                       # Wavenumber
        'fte': grid['fte'],             # Tail factors
        'fttr': grid.get('fttr'),       # (optional)
        'ftwl': grid.get('ftwl')        # (optional)
    }
)

# Most efficient - uses all pre-computed parameters
energy_2d, freq, dirs = converter.to_energy_2d()
hs = converter.compute_wave_parameters()['hs']
```

**Use when**: You want maximum efficiency with cleaner, more readable code

**ww3_params Keys** (all optional):
- `'dden'` or `'dintegral'`: Pre-computed DDEN factors (DTH × DSII × SIG)
- `'wn'` or `'wavenumber'`: Pre-computed wavenumber [1/m]
- `'fte'`: Energy tail factor
- `'fttr'`: Period tail factor (defaults to fte × 0.20 if not provided)
- `'ftwl'`: Wavelength tail factor

---

## API Reference

### Initialization

```python
converter = SpectrumConverter(
    action,                    # Required: (ndir, nfreq) action spectrum
    depth=None,               # Optional: water depth [m]

    # Option A: Specify grid parameters (auto-compute omega, cg)
    fr1=None,                # First frequency [Hz]
    xfr=None,                # Frequency ratio
    nk=None,                 # Number of frequencies
    nth=None,                # Number of directions

    # Option B: Provide pre-computed parameters
    omega=None,              # Angular frequencies [rad/s]
    group_velocity=None,     # Group velocity [m/s]
    dintegral=None,          # DDEN factors (optional)
    wavenumber=None,         # Wavenumber [1/m] (optional)
    directions=None,         # Direction grid [radians]

    # Option C: Pass WW3 parameters via dict (cleaner API)
    ww3_params=None          # Dict with keys: 'dden'/'dintegral', 'wn'/'wavenumber',
                             # 'fte', 'fttr', 'ftwl'
)
```

### Core Methods

#### Convert to Energy Spectra

```python
# 2D energy spectrum: E(f,θ)
energy_2d, frequencies, directions = converter.to_energy_2d()

# 1D frequency spectrum: E(f)
frequencies, e_freq = converter.to_frequency_spectrum_1d()

# 1D directional spectrum: E(θ)
directions, e_dir = converter.to_directional_spectrum_1d()
```

#### Peak Detection

```python
fp = converter.get_peak_frequency()        # Peak frequency [Hz]
theta_peak = converter.get_peak_direction() # Peak direction [rad]
```

#### Wave Parameters

```python
params = converter.compute_wave_parameters()

# Returns dict with:
# - hs: Significant wave height [m]
# - tp: Peak period [s]
# - fp: Peak frequency [Hz]
# - t01, t02, t0m1: Various mean periods [s]
# - thm: Mean direction [rad]
# - ths: Directional spread [rad]
# - wlm: Mean wavelength [m]
# - width: Spectral width
# - m0, m1, m2, m_1: Spectral moments
```

#### Normalization

```python
normalized = converter.normalize(target_energy=5.0)  # Scale to 5 m²
```

---

## Physics

All approaches use the same physics:

```
Action Density: A(k,θ)
    ↓
Coordinate Transform: E(f,θ) = A(k,θ) × SIG × (2π/CG)
    ↓
1D Integration:
  - Frequency:   E(f) = ∑_θ [E(f,θ) × DTH]
  - Directional: E(θ) = ∑_f [E(f,θ) × DSII]
```

Where:
- **SIG = omega**: Angular frequency [rad/s]
- **CG**: Group velocity [m/s]
- **DTH**: Directional bin width [rad]
- **DSII**: Frequency bin width [rad/s]

---

## When to Use Each Pattern

| Pattern | Pros | Cons | Best For |
|---------|------|------|----------|
| **Generic** | Simple, clean, no WW3 knowledge needed | Slightly slower (computes everything) | Non-WW3 data, learning |
| **Pre-ω,cg** | Good balance of simplicity + efficiency | Still computes DSII, wavenumber | Most WW3 workflows |
| **Optimized (explicit)** | Maximum efficiency, skips all redundant computation | Verbose parameter list | High-performance batch processing |
| **Optimized (dict)** | Maximum efficiency + clean code with dict API | Requires pre-computing everything | Production-quality code |

---

## Convenience: Auto-compute Wavenumber and Group Velocity

WW3 parameters like **wavenumber** can be precomputed from grid parameters, but require **water depth**.

`build_ww3_grid()` now supports this automatically:

```python
from wavewatch_python import build_ww3_grid

# Without depth: just basic grid
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)
# Returns: freq, sigma, dsii, dden, fte (no wavenumber/group_velocity)

# With depth: auto-compute wavenumber and group velocity!
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36, depth=100.0)
# Returns: freq, sigma, dsii, dden, fte, wavenumber, group_velocity
```

This eliminates the need to manually compute dispersion for every pattern:

**Before** (manual):
```python
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)
omega = grid['sigma']
wn = np.array([solve_dispersion(omega[i], 100.0)[0] for i in range(len(omega))])
cg = np.array([solve_dispersion(omega[i], 100.0)[1] for i in range(len(omega))])
```

**After** (cleaner):
```python
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36, depth=100.0)
wn = grid['wavenumber']
cg = grid['group_velocity']
```

---

## Examples

### Example 1: Simple Analysis

```python
from wavewatch_python import SpectrumConverter
import numpy as np

# Load your action spectrum
action = np.load('wave_action.npy')  # Shape: (36, 30)

# Simple one-liner: auto-compute everything
converter = SpectrumConverter(action, depth=50.0, fr1=0.04, xfr=1.1, nk=30, nth=36)

# Get results
energy_2d, freq, dirs = converter.to_energy_2d()
print(f"Peak frequency: {converter.get_peak_frequency():.2f} Hz")
print(f"Peak direction: {np.degrees(converter.get_peak_direction()):.1f}°")
```

### Example 2: WW3 Workflow (Using auto-computed grid)

```python
from wavewatch_python import SpectrumConverter, build_ww3_grid
import numpy as np

# From WW3 simulation (NOW with auto-computed wavenumber & group velocity!)
depth = 50.0
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36, depth=depth)
omega = grid['sigma']
cg = grid['group_velocity']  # Auto-computed!

# Your action density from WW3
action = load_from_ww3_netcdf('wave.nc')

# Create converter with pre-computed parameters
converter = SpectrumConverter(action, omega=omega, group_velocity=cg, depth=depth)

# Get wave parameters
params = converter.compute_wave_parameters()
print(f"Hs: {params['hs']:.2f} m")
print(f"Tp: {params['tp']:.2f} s")
print(f"Direction: {np.degrees(params['thm']):.1f}°")
```

### Example 3: Batch Processing (High Performance)

```python
from wavewatch_python import SpectrumConverter, build_ww3_grid
import numpy as np

# Pre-compute ALL parameters once (using convenient depth parameter!)
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36, depth=100.0)
omega = grid['sigma']
wn = grid['wavenumber']  # Auto-computed!
cg = grid['group_velocity']  # Auto-computed!

# Process many spectra efficiently
for i in range(1000):
    action = load_action_spectrum(i)

    converter = SpectrumConverter(
        action,
        omega=omega,
        group_velocity=cg,
        dintegral=grid['dden'],
        wavenumber=wn,
        depth=100.0
    )

    params = converter.compute_wave_parameters(fte=grid['fte'])
    # ... process results
```

### Example 4: Batch Processing (ww3_params Dictionary - Cleaner)

```python
from wavewatch_python import SpectrumConverter, build_ww3_grid
import numpy as np

# Pre-compute ALL parameters once (using convenient depth parameter!)
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36, depth=100.0)
omega = grid['sigma']
wn = grid['wavenumber']  # Auto-computed!
cg = grid['group_velocity']  # Auto-computed!

# Use ww3_params dict for cleaner API
ww3_params = {
    'dden': grid['dden'],
    'wn': wn,
    'fte': grid['fte']
}

# Process many spectra efficiently (cleaner code!)
results = []
for i in range(1000):
    action = load_action_spectrum(i)

    converter = SpectrumConverter(
        action,
        omega=omega,
        group_velocity=cg,
        depth=100.0,
        ww3_params=ww3_params
    )

    params = converter.compute_wave_parameters()
    results.append(params)
```

---

## Comparison: Old vs New

### Before (Fragmented API - Multiple Modules)

```python
# For WW3 data - had to use spectrum.py WaveSpectrum
from wavewatch_python import WaveSpectrum
spectrum = WaveSpectrum(action, depth, omega=omega, dintegral=dden,
                       wavenumber=wn, group_velocity=cg)
params = spectrum.compute_parameters()

# For generic data - had to use spectrum_converter.py functions
from wavewatch_python import action_to_energy_2d, action_to_frequency_spectrum_1d
energy_2d, freq, dirs, _ = action_to_energy_2d(action, omega, cg)

# Different APIs, required learning multiple modules!
```

### After (Unified API)

```python
# ANY data source - single, simple API
from wavewatch_python import SpectrumConverter

converter = SpectrumConverter(action, depth=100, fr1=0.04, xfr=1.1, nk=30, nth=36)
# ... or with pre-computed parameters
converter = SpectrumConverter(action, omega=omega, group_velocity=cg)

# Everything works the same way
energy_2d, freq, dirs = converter.to_energy_2d()
params = converter.compute_wave_parameters()

# Cleaner, simpler, more intuitive!
```

---

## Testing

Run comprehensive tests showing all three approaches give identical results:

```bash
cd /home/user/WW3/python
python test_unified_converter.py
```

Expected output:
```
✓ Generic data works (auto-computes everything)
✓ WW3 data works (uses pre-computed omega, cg)
✓ WW3 optimized works (all parameters pre-computed)
✓ All approaches give consistent results
```

---

## Key Features

✅ **Single API** - No more switching between modules
✅ **Flexible** - Works with any data source
✅ **Efficient** - Uses pre-computed values when available
✅ **Consistent** - All approaches give identical results
✅ **Well-tested** - Comprehensive test coverage
✅ **Physics-correct** - Follows w3iogomd.F90 exactly
✅ **Well-documented** - Clear examples for all use cases

---

## Migration from Old Code

If you're using the old approach, migration is simple:

```python
# Old (spectrum_converter.py functions)
energy_2d, freq, dirs, _ = action_to_energy_2d(action, omega, cg)
freq, e_freq = action_to_frequency_spectrum_1d(action, omega, cg)

# New (unified converter)
converter = SpectrumConverter(action, omega=omega, group_velocity=cg)
energy_2d, freq, dirs = converter.to_energy_2d()
freq, e_freq = converter.to_frequency_spectrum_1d()
```

The old functions still work, but the new `SpectrumConverter` is recommended!

---

## Why This Design?

1. **Physics is Universal** - Coordinate transformation works for any action spectrum
2. **Parameters are Computable** - omega and cg aren't unique to WW3
3. **Single Source of Truth** - One converter, one physics implementation
4. **Flexible** - Accepts any combination of inputs
5. **Efficient** - Only computes what's needed
6. **Simple** - Intuitive API, no confusion about which module to use

## Questions?

See the test file for working examples of all three approaches:
- `test_unified_converter.py` - All use cases demonstrated
