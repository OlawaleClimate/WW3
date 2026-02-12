# WaveWatch III Python Tools

Complete Python implementation for converting WaveWatch III (WW3) action density spectra to wave energy spectra and computing wave parameters.

## 📋 Overview

This package provides:
- **Spectrum Conversion**: Transform action density A(k,θ) in wavenumber-direction space to energy density E(f,θ) in frequency-direction space
- **Wave Parameters**: Compute significant wave height (Hs), peak period (Tp), mean direction, directional spread, and other bulk parameters
- **Dispersion Relations**: Compute wavenumber and group velocity from frequency and water depth
- **Grid Management**: Build and manage WW3 spectral grids with proper frequency bandwidth calculations

## 📁 Directory Structure

```
python/
├── wavewatch_python/              # Main package
│   ├── __init__.py               # Package exports
│   ├── spectrum_converter.py      # Core conversion functions (NEW - UPDATED)
│   ├── spectrum.py               # Full-featured WaveSpectrum class
│   ├── spectrum_ww3_only.py      # Streamlined WW3-only WaveSpectrum
│   ├── dispersion.py             # Dispersion relation solver
│   ├── constants.py              # Physical constants
│   ├── README.md                 # Package documentation
│   ├── SPECTRUM_CONVERTER_GUIDE.md
│   ├── SPECTRUM_WW3_ONLY_GUIDE.md
│   └── WW3_INTEGRATION.md
├── test_spectrum_converter.py    # Test suite for spectrum converter
├── test_spectrum_ww3_only.py     # Test suite for WaveSpectrum
├── CONVERTER_QUICK_START.md      # Quick start guide (UPDATED)
├── SPECTRUM_CONVERTER_UPDATE.md  # Latest implementation changes (NEW)
├── SPECTRUM_MODULES_COMPARISON.md # Comparison of spectrum modules
└── WW3_WAVE_PARAMETERS_COMPLETE_GUIDE.md # Detailed physics explanation
```

## 🚀 Quick Start

### Installation

```bash
# Add to your Python path
import sys
sys.path.insert(0, '/home/user/WW3/python')

# Or install the package
from wavewatch_python import build_ww3_grid, action_to_energy_2d
```

### Basic Usage

```python
import numpy as np
from wavewatch_python import build_ww3_grid, action_to_energy_2d
from wavewatch_python.dispersion import solve_dispersion

# Step 1: Build WW3 grid
grid = build_ww3_grid(fr1=0.04, xfr=1.1, nk=30, nth=36)
omega = grid['sigma']
depth = 100.0

# Step 2: Compute group velocity
cg = np.array([solve_dispersion(omega[i], depth)[1] for i in range(len(omega))])

# Step 3: Load your WW3 action density
action = np.random.rand(36, 30) * 0.001

# Step 4: Convert to energy
energy_2d, freq, dirs, grid_params = action_to_energy_2d(
    action, omega, cg
)

print(f"Peak frequency: {freq[np.argmax(np.sum(energy_2d, axis=0))]:.4f} Hz")
```

## 📚 Documentation

### For Users Starting Out
- **[CONVERTER_QUICK_START.md](CONVERTER_QUICK_START.md)** - Quick reference with common tasks and examples
- **[SPECTRUM_MODULES_COMPARISON.md](SPECTRUM_MODULES_COMPARISON.md)** - Choose between `spectrum.py` and `spectrum_ww3_only.py`

### For Understanding Implementation Details
- **[SPECTRUM_CONVERTER_UPDATE.md](SPECTRUM_CONVERTER_UPDATE.md)** - Latest changes and physics corrections
- **[WW3_WAVE_PARAMETERS_COMPLETE_GUIDE.md](WW3_WAVE_PARAMETERS_COMPLETE_GUIDE.md)** - Complete physics explanation with Fortran references

### For API Reference
- **[wavewatch_python/SPECTRUM_CONVERTER_GUIDE.md](wavewatch_python/SPECTRUM_CONVERTER_GUIDE.md)** - Full spectrum_converter.py API
- **[wavewatch_python/SPECTRUM_WW3_ONLY_GUIDE.md](wavewatch_python/SPECTRUM_WW3_ONLY_GUIDE.md)** - Full WaveSpectrum API
- **[wavewatch_python/README.md](wavewatch_python/README.md)** - Package overview

### For Learning with Examples
- **[test_spectrum_converter.py](test_spectrum_converter.py)** - 6 comprehensive test cases
- **[wavewatch_python/WW3_INTEGRATION.md](wavewatch_python/WW3_INTEGRATION.md)** - Integration examples

## 🔧 Core Functions

### `spectrum_converter.py` - Coordinate Transformation
Converts action density spectra between coordinate systems with proper physics.

```python
from wavewatch_python import (
    build_ww3_grid,
    action_to_energy_2d,
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d,
    get_peak_frequency,
    get_peak_direction
)
```

**Key Functions:**
- `build_ww3_grid(fr1, xfr, nk, nth)` - Build frequency/directional grids
- `action_to_energy_2d(action, omega, cg, ...)` - 2D energy conversion
- `action_to_frequency_spectrum_1d(action, omega, cg, ...)` - 1D frequency spectrum
- `action_to_directional_spectrum_1d(action, omega, cg, ...)` - 1D directional spectrum
- `get_peak_frequency(action, omega, cg)` - Find peak frequency
- `get_peak_direction(action, omega, cg)` - Find peak direction

### `spectrum.py` - Full-Featured Analysis
Flexible class supporting multiple input modes for research and development.

```python
from wavewatch_python import WaveSpectrum

# Method 1: With WW3 parameters (recommended)
spectrum = WaveSpectrum(
    action, depth=100.0,
    omega=omega_ww3, dintegral=dden_ww3,
    wavenumber=wn_ww3, group_velocity=cg_ww3
)

# Method 2: Auto-generate with minimal inputs
spectrum = WaveSpectrum(action, depth=100.0)

params = spectrum.compute_parameters()
```

### `spectrum_ww3_only.py` - Production Ready
Streamlined class requiring explicit WW3 parameters for clarity.

```python
from wavewatch_python.spectrum_ww3_only import WaveSpectrum

spectrum = WaveSpectrum(
    action=action, depth=depth,
    omega=omega_ww3, dintegral=dden_ww3,
    fte=fte_ww3, fttr=fttr_ww3, ftwl=ftwl_ww3,
    wavenumber=wn_ww3, group_velocity=cg_ww3
)

params = spectrum.compute_parameters()
```

### `dispersion.py` - Wave Physics
Compute wavenumber and group velocity from dispersion relations.

```python
from wavewatch_python import solve_dispersion

wn, cg = solve_dispersion(omega, depth)
```

## 📊 Physics: Coordinate Transformation

The core transformation converts action density from wavenumber-direction space to energy density in frequency-direction space:

```
Action Density: A(k,θ) [wavenumber-direction, conserved]
                    ↓
         Jacobian Transformation: ∂k/∂f = 2π/CG
                    ↓
Energy Density: E(f,θ) = A(k,θ) × SIG × (2π/CG)
                [frequency-direction, derived]
```

**Key Physics:**
- **Action-Energy Relationship**: A = F / σ
- **Intrinsic Frequency**: σ = SIG [rad/s]
- **Group Velocity**: CG = dσ/dk [m/s] from dispersion relation
- **Jacobian Factor**: 2π/CG from coordinate transformation
- **Integration Factors**: DTH (directional), DSII (frequency) for 1D spectra

## ✅ Testing

Run the comprehensive test suite:

```bash
python test_spectrum_converter.py
```

This runs 6 tests:
1. 2D action to energy conversion
2. 1D frequency spectrum extraction
3. 1D directional spectrum extraction
4. Peak frequency and direction detection
5. Spectrum normalization
6. Verification against WaveSpectrum class

All tests verify correctness of physics and proper handling of bin widths (DTH, DSII).

## 🔄 Latest Updates

**Version 2.0 - Aligned with Reference Implementation (Feb 2026)**

Key improvements:
- ✅ Added `build_ww3_grid()` function mirroring w3gridmd.F90
- ✅ Fixed DSII computation with special cases for first/last bins
- ✅ Corrected DTH integration in 1D frequency spectrum
- ✅ Corrected DSII integration in 1D directional spectrum
- ✅ Updated peak detection functions with proper API
- ✅ All tests passing with physics-correct implementation
- ✅ Full alignment with reference WW3 implementation

See **[SPECTRUM_CONVERTER_UPDATE.md](SPECTRUM_CONVERTER_UPDATE.md)** for details.

## 🎯 Use Cases

| Task | Module | Example |
|------|--------|---------|
| Quick conversion | `spectrum_converter.py` | `action_to_energy_2d()` |
| Peak detection | `spectrum_converter.py` | `get_peak_frequency()` |
| Wave parameters | `spectrum.py` or `spectrum_ww3_only.py` | `compute_parameters()` |
| Research/testing | `spectrum.py` | Multiple input modes |
| Production WW3 | `spectrum_ww3_only.py` | Explicit parameters |
| Batch processing | `spectrum_converter.py` | Vectorized operations |

## 🌊 Wave Parameters Computed

When using `WaveSpectrum` class:

| Parameter | Symbol | Units | Meaning |
|-----------|--------|-------|---------|
| Significant Wave Height | Hs | m | 4 × √(m₀) |
| Peak Period | Tp | s | 1/fp |
| Peak Frequency | fp | Hz | Frequency of max energy |
| Mean Period 0,1 | Tm01 | s | m₀/m₁ × 2π |
| Mean Period 0,2 | Tm02 | s | 2π × √(m₀/m₂) |
| Mean Period -1,0 | Tm-10 | s | m₋₁/m₀ × 2π |
| Mean Direction | θm | deg | Mean wave direction |
| Directional Spread | σθ | deg | Directional spreading |
| Mean Wavelength | WLM | m | Mean wavelength |
| Peakedness Qp | Qp | - | Spectral peakedness |
| Energy Flux | CGE | W/m | Wave power |

## 📖 References

- **WW3 Source Code**: w3iogomd.F90, w3gridmd.F90, w3dispmd.F90
- **Reference Implementation**: `claude/wave-spectrum-analysis-TgpfA` branch
- **Physics**: Hasselmann et al. (1973), Whitham (1974), IAHR Wave Spectra
- **WW3 Documentation**: [Official WW3 Manual](https://github.com/NOAA-EMC/WW3)

## 🤝 Contributing

When making changes:
1. Update corresponding test cases
2. Verify against reference implementation (w3iogomd.F90)
3. Update documentation
4. Run full test suite: `python test_spectrum_converter.py`
5. Commit with clear message referencing physics or algorithms

## 📝 File Descriptions

### Core Implementation Files
- `spectrum_converter.py` - Pure spectral transformation (UPDATED v2.0)
- `spectrum.py` - Full-featured WaveSpectrum class
- `spectrum_ww3_only.py` - Streamlined WaveSpectrum
- `dispersion.py` - Dispersion relation solvers
- `constants.py` - Physical constants (GRAV, TPI, etc.)

### Documentation Files
- `CONVERTER_QUICK_START.md` - For users new to spectrum conversion (UPDATED)
- `SPECTRUM_CONVERTER_UPDATE.md` - Implementation details of v2.0 (NEW)
- `WW3_WAVE_PARAMETERS_COMPLETE_GUIDE.md` - Physics deep-dive
- `SPECTRUM_MODULES_COMPARISON.md` - When to use which module
- `README.md` (this file) - Master overview

### Test Files
- `test_spectrum_converter.py` - Tests for spectrum_converter functions
- `test_spectrum_ww3_only.py` - Tests for WaveSpectrum class

## ⚠️ Backward Compatibility

Version 2.0 introduces API changes:
- `action_to_energy_2d()` now returns grid_params dict instead of group_velocity
- `get_peak_frequency()` and `get_peak_direction()` require omega and group_velocity

See **[SPECTRUM_CONVERTER_UPDATE.md](SPECTRUM_CONVERTER_UPDATE.md)** for migration guide.

## 🔗 Quick Links

- Start here: **[CONVERTER_QUICK_START.md](CONVERTER_QUICK_START.md)**
- API docs: **[wavewatch_python/SPECTRUM_CONVERTER_GUIDE.md](wavewatch_python/SPECTRUM_CONVERTER_GUIDE.md)**
- Physics: **[WW3_WAVE_PARAMETERS_COMPLETE_GUIDE.md](WW3_WAVE_PARAMETERS_COMPLETE_GUIDE.md)**
- Tests: **[test_spectrum_converter.py](test_spectrum_converter.py)**
