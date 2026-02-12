# WW3 Wave Parameters: Complete Conversion Process

## Overview

This guide explains how WaveWatch III (WW3) converts from **action density in wavenumber-direction space** to **mean wave parameters** used in operational forecasting.

---

## **Phase 1: Fundamental Quantities in WW3**

### Action Density Spectrum A(k,θ)

**Storage Format:**
- Wavenumber-direction space: (k,θ)
- Units: [m²·s·rad⁻¹]
- Conserved quantity during wave propagation

**Relationship to Energy:**
```
A(k,θ) = F(k,θ) / σ(k)

Where:
- F(k,θ) = Energy spectrum in (k,θ) space
- σ(k) = SIG = Intrinsic (angular) frequency [rad/s]
```

### WW3 Grid Parameters (Computed Once at Initialization)

**From w3gridmd.F90:**

```fortran
! Frequency grid (logarithmic spacing)
FR1 = 0.04                                    ! First frequency [Hz]
XFR = 1.1                                     ! Frequency ratio
NK = 30                                       ! Number of frequencies
NTH = 36                                      ! Number of directions

! Angular frequency array (this is σ in formulas)
SIG(IK) = 2π × FR1 × XFR^(IK-1)              ! [rad/s]

! Directional spacing
DTH = 2π / NTH = 0.1745 rad                  ! [radians]

! Frequency bandwidth (from logarithmic grid)
DSII(IK) = (SIG(IK+1) - SIG(IK-1)) / 2       ! [rad/s]
```

**From w3initmd.F90 (WAVNU1 subroutine):**

For each frequency and water depth location:
```fortran
INPUT:  SIG(IK), DEPTH(ISEA)
OUTPUT: WN(IK,ISEA), CG(IK,ISEA)

! These satisfy the dispersion relation:
! ω² = g·k·tanh(k·h)
```

---

## **Phase 2: Coordinate Transformation (k,θ) → (f,θ)**

### The Challenge

- WW3 stores action in **wavenumber-direction space (k,θ)**
- Users need energy in **frequency-direction space (f,θ)**
- These are DIFFERENT coordinate systems!

### The Solution: Jacobian Transformation

**Dispersion Relation Jacobian:**
```
∂k/∂f = 2π / CG
```

Where:
- From the dispersion relation: ω = 2πf = σ(k)
- Group velocity: CG = dω/dk = dσ/dk

**Energy Transformation:**
```
E(f,θ) df dθ = F(k,θ) dk dθ

Therefore:
E(f,θ) = F(k,θ) × |∂k/∂f| = F(k,θ) × (2π/CG)
```

### Complete Conversion Formula

```
E(f,θ) = F(k,θ) × (2π/CG)
        = A(k,θ) × σ × (2π/CG)
        = A(k,θ) × SIG × (2π/CG)
```

**This combines:**
1. **Action-to-energy relationship:** × σ (intrinsic frequency)
2. **Coordinate transformation Jacobian:** × (2π/CG)

---

## **Phase 3: Spectral Moment Accumulation**

### Integration Process

From the 2D energy spectrum E(f,θ), compute moments:

```python
m_n = ∫∫ E(f,θ) × f^n df dθ

# In discrete form (WW3 implementation):
for each frequency bin ik:
    # 1. Sum action over all directions
    action_band = ∑_θ A(θ,ik)

    # 2. Convert to energy (apply the conversion formula)
    energy_band = action_band × SIG(ik) × (2π/CG(ik))

    # 3. Accumulate spectral moments
    m0 += energy_band
    m1 += energy_band × SIG(ik)
    m2 += energy_band × SIG(ik)²
    m_1 += energy_band / SIG(ik)

    # 4. Directional moments
    mom_x += ∑_θ [A(θ,ik) × cos(θ)] × SIG × (2π/CG)
    mom_y += ∑_θ [A(θ,ik) × sin(θ)] × SIG × (2π/CG)

    # 5. Wavelength moment (using wavenumber)
    mom_wn += energy_band / WN(ik)

# 6. Tail extension (Pierson-Moskowitz f⁻⁵ tail)
# Beyond cutoff frequency, use pre-computed factors:
m0 += FTE × energy_tail
m1 += FTE × 0.20 × energy_tail
m_1 += FTTR × energy_tail
mom_wn += FTWL × energy_tail
```

### Spectral Moments

**Definition:**
```
m₀ = Total energy [m²]
m₁ = First moment (frequency-weighted) [m²·s]
m₂ = Second moment (frequency²-weighted) [m²·s²]
m₋₁ = Inverse moment [m²·s]
mom_x = E-W directional moment
mom_y = N-S directional moment
mom_wn = Wavelength moment
```

---

## **Phase 4: Derive Mean Wave Parameters**

### From Spectral Moments to Wave Parameters

```python
# WAVE HEIGHTS
HS = 4 × √(m₀)                          # Significant wave height [m]

# WAVE PERIODS
T01 = 2π × (m₀ / m₁)                   # Mean period [s]
T02 = 2π × √(m₀ / m₂)                  # Zero-crossing period [s]
T0M1 = 2π × (m₋₁ / m₀)                 # Energy period [s]
TP = 2π / ω_peak = 2π / SIG(peak_ik)   # Peak period [s]

# DIRECTIONAL PROPERTIES
THM = atan2(mom_y, mom_x)               # Mean direction [radians]
                                         # Convert to degrees: THM × 180/π
THS = √(2 × (1 - √((mom_x² + mom_y²) / m₀²)))  # Directional spread [radians]

# WAVELENGTH
WLM = 2π × (mom_wn / m₀)                # Mean wavelength [m]

# SPECTRAL WIDTH
WIDTH = √(m₂×m₀/m₁² - 1)                # Spectral narrowness (dimensionless)
```

---

## **Phase 5: Complete Processing Chain**

```
┌──────────────────────────────────────────────┐
│ WW3 Action Density: A(k,θ)                   │
│ • Storage: Wavenumber-direction space        │
│ • Units: [m²·s·rad⁻¹]                       │
│ • Conserved quantity                         │
└──────────┬───────────────────────────────────┘
           │
           ↓
┌──────────────────────────────────────────────────────┐
│ Energy in (k,θ) Space: F(k,θ) = A(k,θ) × SIG       │
│ • Computed from action × intrinsic frequency        │
│ • Units: [m²]                                       │
└──────────┬───────────────────────────────────────────┘
           │
           ↓ Coordinate Transformation
           │ Using Jacobian: 2π/CG
           │ (From dispersion relation)
           │
┌──────────────────────────────────────────────────────┐
│ Energy in (f,θ) Space:                               │
│ E(f,θ) = A(k,θ) × SIG × (2π/CG)                     │
│ • Output in frequency-direction space               │
│ • Units: [m²/Hz/rad]                                │
└──────────┬───────────────────────────────────────────┘
           │
           ↓ Integrate Over Frequencies
           │
┌──────────────────────────────────────────────────────┐
│ Spectral Moments: m₀, m₁, m₂, m₋₁                   │
│ Directional Moments: mom_x, mom_y, mom_wn           │
│ • Accumulated from 2D energy spectrum               │
│ • Include tail extension (f⁻⁵ tail)                 │
└──────────┬───────────────────────────────────────────┘
           │
           ↓ Apply Spectral Formulas
           │
┌──────────────────────────────────────────────────────┐
│ MEAN WAVE PARAMETERS:                                │
│ • HS, T01, T02, T0M1, TP (Heights & Periods)       │
│ • THM, THS (Directional Properties)                 │
│ • WLM (Mean Wavelength)                             │
│ • WIDTH (Spectral Narrowness)                       │
└──────────────────────────────────────────────────────┘
```

---

## **Important Constants & Variables**

### Physical Constants
```python
GRAV = 9.81 m/s²               # Gravitational acceleration
TPI = 2π ≈ 6.2832              # Radian conversion
TPIINV = 1/(2π) ≈ 0.15915      # Inverse of 2π
SMALL = 1e-20                   # Floor value for division
```

### Grid Parameters (Fixed During Simulation)
```python
FR1 = 0.04          # First frequency [Hz]
XFR = 1.1           # Frequency ratio
NK = 30             # Number of frequencies
NTH = 36            # Number of directions
SIG = ω [rad/s]     # Angular frequency array
DTH [rad]           # Directional spacing = 2π/NTH
DSII [rad/s]        # Frequency bandwidth
DDEN [rad²/s]       # Full integration factor = DTH × DSII × SIG
```

### Dispersion Variables (Fixed During Simulation)
```python
WN(ik) [1/m]        # Wavenumber at each frequency
CG(ik) [m/s]        # Group velocity at each frequency
                     # Both computed from: ω² = g·k·tanh(k·h)
```

### Tail Factors (Fixed During Simulation)
```python
FTE = 0.25 × SIG(NK) × DTH × SIG(NK)      # Energy tail
FTTR = 0.20 × DTH × SIG(NK)               # Period tail
FTWL = (GRAV/6) / SIG(NK) × DTH × SIG(NK) # Wavelength tail
```

---

## **Python Implementation**

```python
from wavewatch_python import (
    action_to_energy_2d,
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d,
    GRAV, TPI
)
from wavewatch_python.dispersion import solve_dispersion

# WW3 pre-computed parameters
omega = 2*np.pi * 0.04 * 1.1**np.arange(30)  # SIG
depth = 100.0

# Compute dispersion
wn = np.zeros(30)
cg = np.zeros(30)
for ik in range(30):
    wn[ik], cg[ik] = solve_dispersion(omega[ik], depth)

# Your action density from WW3
action = np.array(...)  # Shape: (36, 30)

# Step 1: Convert to 2D energy
energy_2d, freq, dirs, _ = action_to_energy_2d(
    action, omega, cg
)

# Step 2: Compute moments
m0 = np.sum(energy_2d)
m1 = np.sum(energy_2d[:, :] * omega[np.newaxis, :])

# Step 3: Compute wave parameters
hs = 4 * np.sqrt(m0)
t01 = TPI * (m0 / m1)

print(f"HS: {hs:.2f} m")
print(f"T01: {t01:.2f} s")
```

---

## **Key Physics Insights**

1. **Action is Conserved**
   - Doesn't change during wave propagation
   - Fundamental quantity WW3 solves for

2. **Energy is Not Conserved**
   - Changes due to depth, current, dissipation
   - Derived quantity for user output

3. **SIG (σ) is Essential**
   - Intrinsic frequency in the wave frame
   - Must be included in conversion formula

4. **Group Velocity is Critical**
   - Controls energy transport speed
   - Depth-dependent
   - Essential for Jacobian calculation

5. **Coordinate Transformation is Fundamental**
   - Converting from (k,θ) to (f,θ) requires Jacobian
   - Cannot be ignored or approximated away

---

## **References**

- WW3 Source Code: w3gridmd.F90, w3initmd.F90, w3iogomd.F90
- Hasselmann et al. (1973): Wave spectral energy balance
- Pierson & Moskowitz (1964): Tail frequency extension
- IAHR (1989): Wave analysis standards
- Whitham (1974): Linear and nonlinear waves
