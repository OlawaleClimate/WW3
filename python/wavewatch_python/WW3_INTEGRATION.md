# Using WW3 Pre-computed Grid Parameters with Python Package

## Overview

The WaveWatch Python package can now accept **pre-computed grid parameters from a WW3 simulation** instead of auto-generating them. This ensures:

✓ **Exact consistency** with WW3 model calculations
✓ **No duplication** of grid computations
✓ **Direct compatibility** with WW3 output data
✓ **Maximum efficiency** by reusing WW3's pre-computed values

---

## WW3 Grid Parameters to Provide

When you have a WW3 simulation, extract these pre-computed values at grid initialization:

### 1. Angular Frequency Array (SIG)
```fortran
! From w3gridmd.F90:
omega(1:NK) = SIG(1:NK)  [rad/s]
```
**In Python:**
```python
omega_ww3 = SIG_from_ww3  # numpy array of size NK
```

### 2. Integration Factor (DDEN)
```fortran
! From w3gridmd.F90:
DDEN(IK) = DTH × DSII(IK) × SIG(IK)
```
**In Python:**
```python
dintegral_ww3 = DDEN_from_ww3  # numpy array of size NK
```

### 3. Tail Factors (Pre-computed in w3gridmd.F90)
```fortran
FTE   = 0.25 * SIG(NK) * DTH * SIG(NK)
FTTR  = 0.20 * DTH * SIG(NK)
FTWL  = (GRAV/6) / SIG(NK) * DTH * SIG(NK)
```
**In Python:**
```python
fte_ww3 = 0.25 * omega_ww3[-1] * dth * omega_ww3[-1]
fttr_ww3 = 0.20 * dth * omega_ww3[-1]
ftwl_ww3 = (9.81 / 6.0) / omega_ww3[-1] * dth * omega_ww3[-1]
```

### 4. Wavenumber Array (WN)
```fortran
! From w3initmd.F90:1370-1399
! For each sea point and frequency:
CALL WAVNU1(SIG(IK), DEPTH, WN(IK,ISEA), CG(IK,ISEA))
```
**In Python:**
```python
wn_ww3 = WN_from_ww3[:, isea]  # numpy array of size NK
```

### 5. Group Velocity Array (CG)
```fortran
! From WAVNU1 subroutine (w3dispmd.F90)
CG(IK,ISEA) = group velocity at frequency IK
```
**In Python:**
```python
cg_ww3 = CG_from_ww3[:, isea]  # numpy array of size NK
```

---

## Complete Usage Example

### Method 1: Using WW3 Pre-computed Parameters (RECOMMENDED)

```python
import numpy as np
from wavewatch import WaveSpectrum
from wavewatch.dispersion import solve_dispersion

# Get action density spectrum from WW3 output
action_ww3 = load_ww3_spectrum()  # shape (36, 30) = (NTH, NK)

# Get pre-computed WW3 grid parameters
omega_ww3 = load_ww3_omega()        # shape (30,)
dintegral_ww3 = load_ww3_dden()     # shape (30,)
wn_ww3 = load_ww3_wavenumber()      # shape (30,)
cg_ww3 = load_ww3_group_velocity()  # shape (30,)

# Get or compute tail factors
fte_ww3 = load_or_compute_fte()
fttr_ww3 = load_or_compute_fttr()
ftwl_ww3 = load_or_compute_ftwl()

# Get water depth for this location
depth = 100.0

# Initialize spectrum with WW3 parameters
spectrum = WaveSpectrum(
    action_ww3,
    depth=depth,
    # All WW3 pre-computed parameters
    omega=omega_ww3,
    dintegral=dintegral_ww3,
    fte=fte_ww3,
    fttr=fttr_ww3,
    ftwl=ftwl_ww3,
    wavenumber=wn_ww3,
    group_velocity=cg_ww3
)

# Compute wave parameters
params = spectrum.compute_parameters()

print(f"HS = {params['hs']:.3f} m")
print(f"T01 = {params['t01']:.2f} s")
print(f"T02 = {params['t02']:.2f} s")
```

### Method 2: Using Raw Grids (Auto-computed)

```python
# Alternative: If you only have frequencies and directions
spectrum = WaveSpectrum(
    action,
    depth=100.0,
    frequencies=freq_array,
    directions=dir_array
)
```

---

## Extracting WW3 Parameters

### From WW3 Model Code

**In w3gridmd.F90** (at grid initialization):
```fortran
! Extract these arrays after W3GRID initialization
SIG(1:NK)           ! Angular frequencies
DSII(1:NK)          ! Frequency bandwidths
DDEN(1:NK)          ! Integration factors
DTH                 ! Directional step
```

**In w3initmd.F90** (at model initialization):
```fortran
! After calling WAVNU1:
WN(1:NK, 1:NSEA)    ! Wavenumber array
CG(1:NK, 1:NSEA)    ! Group velocity array
```

### From WW3 Output/Restart Files

Many WW3 configurations output these as part of model state. Check:
- NetCDF output files (dimensions and variables)
- Restart files (.ww3 format)
- Grid definition files

### Computing on-the-fly

If not available, compute them:

```python
import numpy as np
from wavewatch.dispersion import solve_dispersion

# Create omega array
fr1 = 0.04
xfr = 1.1
nfreq = 30
omega = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)

# Compute tail factors
dth = 2.0 * np.pi / 36  # NTH = 36
fte = 0.25 * omega[-1] * dth * omega[-1]
fttr = 0.20 * dth * omega[-1]
ftwl = (9.81 / 6.0) / omega[-1] * dth * omega[-1]

# Compute dintegral
dsii = np.zeros(nfreq)
dsii[0] = (omega[1] - omega[0]) / 2.0
for i in range(1, nfreq-1):
    dsii[i] = (omega[i+1] - omega[i-1]) / 2.0
dsii[-1] = (omega[-1] - omega[-2]) / 2.0
dintegral = dth * dsii * omega

# Compute WN and CG
depth = 100.0
wn = np.zeros(nfreq)
cg = np.zeros(nfreq)
for ik in range(nfreq):
    wn[ik], cg[ik] = solve_dispersion(omega[ik], depth)

# Now pass to WaveSpectrum
spectrum = WaveSpectrum(
    action,
    depth=depth,
    omega=omega,
    dintegral=dintegral,
    fte=fte,
    fttr=fttr,
    ftwl=ftwl,
    wavenumber=wn,
    group_velocity=cg
)
```

---

## WW3 File Reference

### Key WW3 Source Files

| File | Lines | Purpose |
|------|-------|---------|
| w3gridmd.F90 | 3444-3448 | Compute DDEN, FTE, FTTR, FTWL |
| w3initmd.F90 | 1370-1399 | Compute WN and CG via WAVNU1 |
| w3dispmd.F90 | Various | WAVNU1 subroutine (dispersion solver) |
| w3iogomd.F90 | 1552 | Use DDEN and CG in calculation |

### Reading WW3 Grid Data

```python
# From NetCDF file (if WW3 output includes it)
import xarray as xr

ds = xr.open_dataset('ww3_grid.nc')
omega = ds['SIG'].values  # Assume variable named SIG
wn = ds['WN'].values
cg = ds['CG'].values

# From Fortran unformatted file (WW3 native)
import f90

with open('ww3.grid', 'rb') as f:
    # Read Fortran record structure
    # (Implementation depends on WW3 version)
    pass
```

---

## Parameter Validation

Before using WW3 parameters, validate them:

```python
import numpy as np

def validate_ww3_parameters(omega, dintegral, fte, fttr, ftwl, wn, cg):
    """Validate WW3 parameters."""

    checks = []

    # Check shapes
    nfreq = len(omega)
    assert len(dintegral) == nfreq, "Shape mismatch: dintegral"
    assert len(wn) == nfreq, "Shape mismatch: wavenumber"
    assert len(cg) == nfreq, "Shape mismatch: group_velocity"

    # Check omega is increasing
    assert np.all(np.diff(omega) > 0), "omega must be monotonically increasing"
    checks.append("✓ omega monotonically increasing")

    # Check dintegral is positive
    assert np.all(dintegral > 0), "dintegral must be positive"
    checks.append("✓ dintegral all positive")

    # Check wavenumber is positive and increasing
    assert np.all(wn > 0), "wavenumber must be positive"
    assert np.all(np.diff(wn) > 0), "wavenumber must be increasing"
    checks.append("✓ wavenumber positive and increasing")

    # Check group velocity is positive
    assert np.all(cg > 0), "group_velocity must be positive"
    checks.append("✓ group_velocity all positive")

    # Check tail factors are reasonable
    assert fte > 0 and fte < 1.0, "FTE should be between 0 and 1"
    assert fttr > 0 and fttr < 1.0, "FTTR should be between 0 and 1"
    assert ftwl > 0 and ftwl < 10.0, "FTWL should be reasonable"
    checks.append("✓ tail factors reasonable")

    # Check Airy relation (deep water approximation)
    omega_check = np.sqrt(9.81 * wn)
    error = np.max(np.abs(omega - omega_check)) / np.max(omega)
    if error > 0.1:  # 10% tolerance for shallow water
        checks.append(f"⚠ Airy relation error {error*100:.1f}% (OK for shallow water)")
    else:
        checks.append("✓ Airy relation satisfied")

    for check in checks:
        print(check)

    return True

# Use it
validate_ww3_parameters(omega, dintegral, fte, fttr, ftwl, wn, cg)
```

---

## Integration with WW3 Workflow

### Full Workflow

```
WW3 Simulation Run
    ↓
Output A(θ,f) spectrum
    ↓
Extract SIG, DDEN, FTE, etc.
    ↓
Extract WN, CG from WAVNU1
    ↓
Pass all to Python WaveSpectrum
    ↓
Compute HS, T01, T02, etc.
    ↓
Validate/Compare with WW3 output
```

### Processing WW3 Simulation Results

```python
import numpy as np
from wavewatch import WaveSpectrum

def process_ww3_output(ww3_spectrum, ww3_params, depth):
    """
    Process WW3 simulation output.

    Args:
        ww3_spectrum: Action density A(θ,f) from WW3
        ww3_params: Dict with 'omega', 'dintegral', 'fte', 'fttr', 'ftwl', 'wn', 'cg'
        depth: Water depth

    Returns:
        Wave parameters from Python calculation
    """

    spectrum = WaveSpectrum(
        ww3_spectrum,
        depth=depth,
        **ww3_params  # Unpack all parameters
    )

    params = spectrum.compute_parameters()

    return params

# Example usage
params = process_ww3_output(
    action_from_ww3,
    {
        'omega': omega_ww3,
        'dintegral': dintegral_ww3,
        'fte': fte_ww3,
        'fttr': fttr_ww3,
        'ftwl': ftwl_ww3,
        'wavenumber': wn_ww3,
        'group_velocity': cg_ww3
    },
    depth=100.0
)
```

---

## Summary

**Use WW3 pre-computed parameters by passing:**

| Parameter | Source | WW3 Variable |
|-----------|--------|--------------|
| `omega` | w3gridmd | SIG |
| `dintegral` | w3gridmd | DDEN = DTH × DSII × SIG |
| `fte` | w3gridmd (line 3444) | FTE |
| `fttr` | w3gridmd (line 3447) | FTTR |
| `ftwl` | w3gridmd (line 3448) | FTWL |
| `wavenumber` | w3initmd + WAVNU1 | WN |
| `group_velocity` | w3initmd + WAVNU1 | CG |

**All passed to WaveSpectrum constructor:**

```python
spectrum = WaveSpectrum(
    action,
    depth=depth,
    omega=omega_ww3,
    dintegral=dintegral_ww3,
    fte=fte_ww3,
    fttr=fttr_ww3,
    ftwl=ftwl_ww3,
    wavenumber=wn_ww3,
    group_velocity=cg_ww3
)
```

This ensures **100% consistency** with WW3's calculations!
