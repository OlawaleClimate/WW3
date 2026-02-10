# Variables Required for SWH Computation in WaveWatch III
## Beyond the Action Density Spectrum

---

## Executive Summary

While the **action density spectrum** `A(θ,f)` is the primary input to SWH calculation, WW3 requires numerous additional variables for accurate wave parameter computation. These variables fall into three categories:

1. **Grid Parameters** - Spectral discretization information
2. **Physical Properties** - Frequency, wavenumber, and velocity relationships
3. **Environment Variables** - Water depth and bathymetry information

---

## 1. Grid and Discretization Parameters

### 1.1 Frequency Grid Parameters

| Variable | Type | Units | Purpose | From File |
|----------|------|-------|---------|-----------|
| **NK** | Integer | — | Number of frequency bins (typically 30-35) | w3gridmd.F90 |
| **NTH** | Integer | — | Number of directional bins (typically 36-72) | w3gridmd.F90 |
| **SIG(IK)** | Real Array | rad/s | Angular frequency at bin IK | w3gridmd.F90 |
| **DSII(IK)** | Real Array | rad/s | Frequency bandwidth for bin IK | w3gridmd.F90 |
| **DDEN(IK)** | Real Array | — | Full discretization factor (DTH × DSII × SIG) | w3gridmd.F90 |
| **DTH** | Real | radians | Directional step size (2π/NTH) | w3gridmd.F90 |
| **FR1** | Real | Hz | First (minimum) frequency | Grid file input |
| **XFR** | Real | — | Frequency ratio (typically 1.1) | Grid file input |

### 1.2 Integration Factors (Pre-computed at Grid Init)

| Variable | Value | Formula | Purpose |
|----------|-------|---------|---------|
| **FTE** | Factor | 0.25 × SIG(NK) × DTH × SIG(NK) | Energy tail contribution |
| **FTF** | Factor | 0.20 × DTH × SIG(NK) | Frequency tail factor |
| **FTTR** | Factor | FTF | Period tail factor |
| **FTWL** | Factor | (GRAV/6) × DTH / SIG(NK) | Wavelength tail factor |
| **FTWN** | Factor | 0.20 × √GRAV × DTH × SIG(NK) | Wavenumber tail factor |

---

## 2. Physical and Spectral Parameters

### 2.1 Wavenumber and Velocity Arrays

| Variable | Type | Units | Purpose | Computed By |
|----------|------|-------|---------|-------------|
| **WN(IK, ISEA)** | Real 2D Array | 1/m | Wavenumber at frequency IK and location ISEA | WAVNU1/WAVNU2 (w3dispmd.F90) |
| **CG(IK, ISEA)** | Real 2D Array | m/s | Group velocity at frequency IK and location ISEA | WAVNU1/WAVNU2 (w3dispmd.F90) |

#### Dispersion Relation
```
Deep Water:     σ² = g·k
Shallow Water:  σ² = g·k·tanh(k·h)

Group Velocity (deep water): CG = σ/(2k) = √(g/k)

WAVNU1 iteratively solves for k given σ using Newton-Raphson
```

### 2.2 Spectral Moment Arrays (Accumulated During Loop)

| Variable | Type | Units | Formula | Meaning |
|----------|------|-------|---------|---------|
| **ET(JSEA)** | Real Array | m² | ∑∫E(f)df | Total energy (M₀) |
| **ET1(JSEA)** | Real Array | m²·s | ∑∫f·E(f)df | First moment (M₁) |
| **ET02(JSEA)** | Real Array | m²·s² | ∑∫f²·E(f)df | Second moment (M₂) |
| **ETR(JSEA)** | Real Array | m²·s | ∑∫E(f)/f·df | Inverse moment (M₋₁) |
| **EWN(JSEA)** | Real Array | m³ | ∑∫E(f)/k·df | Wavelength moment |

### 2.3 Directional Moment Arrays

| Variable | Type | Formula | Purpose |
|----------|------|---------|---------|
| **ETX(JSEA)** | Real Array | ∑∫E(f,θ)·cos(θ)·df·dθ | East-West directional moment |
| **ETY(JSEA)** | Real Array | ∑∫E(f,θ)·sin(θ)·df·dθ | North-South directional moment |
| **EET1(JSEA)** | Real Array | ∑∫E(f)²·f·df / DSII(IK) | For peakedness (QP) calculation |

---

## 3. Environmental and Location Parameters

### 3.1 Water Depth and Bathymetry

| Variable | Type | Units | Source | Impact |
|----------|------|-------|--------|--------|
| **DW(ISEA)** | Real Array | m | w3adatmd.F90 | Water depth at each point |
| **DEPTH** | Real | m | Local variable | Passed to WAVNU1 for dispersion |
| **h** | Real | m | In dispersion relation | Affects wavenumber and group velocity |

**Effect on Calculation:**
- Determines which dispersion relation to use:
  - Deep water (h > λ/2): σ = √(gk)
  - Intermediate (λ/20 < h < λ/2): σ = √(gk·tanh(kh))
  - Shallow water (h < λ/20): σ ≈ √(gh)·k
- Changes group velocity: CG = ∂σ/∂k (depth-dependent)
- Different depths → different wavenumbers at same frequency

### 3.2 Location Parameters

| Variable | Type | Purpose |
|----------|------|---------|
| **JSEA** | Integer | Sea point index in local array |
| **ISEA** | Integer | Global sea point index |
| **IX, IY** | Integer | Grid indices (x, y position) |
| **MAPSF(ISEA, 1/2)** | Integer Array | Maps sea points to grid (x/y) |

---

## 4. Output Array Parameters

### 4.1 Wave Parameters Computed from Spectral Moments

| Parameter | Variable | Formula | Units | Code Line |
|-----------|----------|---------|-------|-----------|
| **Significant Wave Height** | HS | 4√ET | m | 1997 |
| **Zero-crossing Period** | T02 | 2π√(ET/ET02) | s | 2039 |
| **Mean Period** | T01 | 2π(ET/ET1) | s | 2040 |
| **Energy Period** | T0M1 | 2π(ETR/ET) | s | 2006 |
| **Mean Wavelength** | WLM | 2π(EWN/ET) | m | 2005 |
| **Mean Direction** | THM | atan2(ETY, ETX) | rad | 2019 |
| **Directional Spread** | THS | Longuet-Higgins formula | rad | 2007 |
| **Peakedness** | QP | 2·EET1/ET² | — | 2004 |

### 4.2 Intermediate Computed Values

| Variable | Type | Purpose | Code Location |
|----------|------|---------|----------------|
| **EBD(IK, JSEA)** | Real 2D Array | Energy in frequency band IK | w3iogomd:1553 |
| **AB(JSEA)** | Real Array | Action integrated over directions | w3iogomd:1513 |
| **ABX(JSEA), ABY(JSEA)** | Real Arrays | Directional components of action | w3iogomd:1514-1515 |
| **FACTOR** | Real | Conversion factor (DDEN/CG) | w3iogomd:1552 |

---

## 5. Constants and Conversion Factors

### 5.1 Physical Constants

```fortran
! From constants.F90

REAL, PARAMETER :: GRAV = 9.81         ! Gravitational acceleration (m/s²)
REAL, PARAMETER :: DWAT = 1025.0       ! Water density (kg/m³)
REAL, PARAMETER :: DAIR = 1.225        ! Air density (kg/m³)
REAL, PARAMETER :: TPI = 2.0*PI        ! 2π (6.283185...)
REAL, PARAMETER :: TPIINV = 1.0/TPI    ! 1/(2π) (0.159155...)
REAL, PARAMETER :: RADE = 180.0/PI     ! Conversion radians to degrees
```

### 5.2 Conversion Factors

| Factor | Expression | Purpose |
|--------|-----------|---------|
| **Radian to Hz** | σ/(2π) | Convert angular frequency to frequency |
| **Frequency to Period** | 1/f or 2π/σ | Convert frequency to period |
| **Wavenumber to Wavelength** | 2π/k | Convert wavenumber to wavelength |

---

## 6. Detailed Computation Flow with All Variables

### Step-by-Step Variable Usage

```fortran
! ============================================================
! SUBROUTINE W3OUTG - Main Calculation
! ============================================================

SUBROUTINE W3OUTG (A, FLPART, FLOUTG, FLOUTG2)

    ! INPUT VARIABLES
    REAL :: A(NTH, NK, 0:NSEAL)              ! Action density spectrum

    ! GRID PARAMETERS
    USE W3GDATMD                             ! Access to:
        ! - NK, NTH (dimensions)
        ! - SIG(IK), DSII(IK), DDEN(IK), DTH (frequency grid)
        ! - FTE, FTTR, FTWL, FTWN (tail factors)

    ! PHYSICAL ARRAYS
    USE W3ADATMD                             ! Access to:
        ! - CG(IK, ISEA) (group velocity)
        ! - WN(IK, ISEA) (wavenumber)
        ! - DW(ISEA) (water depth)
        ! - Output arrays: HS, T01, T02, T0M1, etc.

    ! LOOP OVER FREQUENCY BINS
    DO IK = 1, NK

        ! Get parameters for this frequency
        SIGMA = SIG(IK)                      ! Angular frequency
        FREQ_BANDWIDTH = DSII(IK)            ! Frequency bandwidth
        INTEGRATION_FACTOR = DDEN(IK)        ! Full integration factor

        ! LOOP OVER DIRECTIONS
        DO ITH = 1, NTH

            DO JSEA = 1, NSEAL

                ! Get depth for dispersion relation
                ISEA = MAP_JSEA_TO_ISEA(JSEA)
                DEPTH = DW(ISEA)

                ! Get action density at (f, θ)
                ACTION = A(ITH, IK, JSEA)

                ! Accumulate directional integral
                AB(JSEA) = AB(JSEA) + ACTION
                ABX(JSEA) = ABX(JSEA) + ACTION*COS(THETA(ITH))
                ABY(JSEA) = ABY(JSEA) + ACTION*SIN(THETA(ITH))

            END DO
        END DO

        ! ENERGY CONVERSION
        DO JSEA = 1, NSEAL

            ISEA = MAP_JSEA_TO_ISEA(JSEA)

            ! Convert action to energy using group velocity
            CG_VALUE = CG(IK, ISEA)           ! Group velocity (depth-dependent)
            FACTOR = DDEN(IK) / CG_VALUE      ! Conversion factor
            EBD(IK, JSEA) = AB(JSEA) * FACTOR ! Energy in band

            ! Accumulate spectral moments
            ET(JSEA) = ET(JSEA) + EBD(IK, JSEA)
            ET1(JSEA) = ET1(JSEA) + EBD(IK, JSEA)*SIG(IK)
            ET02(JSEA) = ET02(JSEA) + EBD(IK, JSEA)*SIG(IK)**2
            ETR(JSEA) = ETR(JSEA) + EBD(IK, JSEA)/SIG(IK)
            EWN(JSEA) = EWN(JSEA) + EBD(IK, JSEA)/WN(IK, ISEA)
            ETX(JSEA) = ETX(JSEA) + ABX(JSEA)*FACTOR
            ETY(JSEA) = ETY(JSEA) + ABY(JSEA)*FACTOR

        END DO
    END DO

    ! ADD TAIL CONTRIBUTION
    DO JSEA = 1, NSEAL

        ISEA = MAP_JSEA_TO_ISEA(JSEA)
        EBAND = AB(JSEA) / CG(NK, ISEA)

        ! Add Pierson-Moskowitz tail to moments
        ET(JSEA) = ET(JSEA) + FTE*EBAND
        ET1(JSEA) = ET1(JSEA) + FT1*EBAND
        ET02(JSEA) = ET02(JSEA) + 0.5*SIG(NK)**4*DTH*EBAND
        ETR(JSEA) = ETR(JSEA) + FTTR*EBAND

    END DO

    ! COMPUTE WAVE PARAMETERS FROM MOMENTS
    DO JSEA = 1, NSEAL

        ! Significant wave height from M₀
        HS(JSEA) = 4.0 * SQRT(ET(JSEA))

        ! Periods from spectral moments
        IF (ET02(JSEA) > 1.E-7 .AND. ET(JSEA) > 0) THEN
            T02(JSEA) = TPI * SQRT(ET(JSEA) / ET02(JSEA))
            T01(JSEA) = TPI * ET(JSEA) / ET1(JSEA)
        ELSE
            T02(JSEA) = TPI / SIG(NK)
            T01(JSEA) = T02(JSEA)
        END IF

        ! Energy period
        IF (ET(JSEA) > 1.E-7) THEN
            T0M1(JSEA) = ETR(JSEA) / ET(JSEA) * TPI
        ELSE
            T0M1(JSEA) = TPI / SIG(NK)
        END IF

        ! Mean direction
        IF (ABS(ETX(JSEA)) + ABS(ETY(JSEA)) > 1.E-7) THEN
            THM(JSEA) = ATAN2(ETY(JSEA), ETX(JSEA))
        ELSE
            THM(JSEA) = 0.
        END IF

        ! Mean wavelength
        IF (ET(JSEA) > 1.E-7) THEN
            WLM(JSEA) = EWN(JSEA) / ET(JSEA) * TPI
        ELSE
            WLM(JSEA) = 0.
        END IF

    END DO

END SUBROUTINE
```

---

## 7. Summary Table: All Variables Required

| Category | Variable | Type | Source | Used For |
|----------|----------|------|--------|----------|
| **Input** | A(NTH,NK,JSEA) | Real 3D | Model state | Action density spectrum |
| **Grid** | NK, NTH | Integer | Grid file | Array dimensions |
| **Grid** | SIG(IK) | Real Array | w3gridmd | Angular frequency |
| **Grid** | DSII(IK), DDEN(IK) | Real Array | w3gridmd | Integration factors |
| **Grid** | DTH | Real | w3gridmd | Directional step |
| **Grid** | FTE, FTTR, FTWL | Real | w3gridmd | Tail factors |
| **Physics** | WN(IK, ISEA) | Real 2D | WAVNU1 | Wavenumber (from dispersion) |
| **Physics** | CG(IK, ISEA) | Real 2D | WAVNU1 | Group velocity (depth-dependent) |
| **Environment** | DW(ISEA) | Real Array | w3adatmd | Water depth |
| **Constants** | GRAV, TPI, RADE | Real | constants.F90 | Conversions |
| **Intermediate** | ET, ET1, ET02, ETR | Real Array | Computed | Spectral moments |
| **Intermediate** | ETX, ETY, EWN | Real Array | Computed | Directional/wavelength moments |
| **Output** | HS, T01, T02, T0M1 | Real Array | Computed | Wave parameters |

---

## 8. Critical Dependencies

### 8.1 Dispersion Relation Dependency

```
Water Depth (DW) → WAVNU1 subroutine
    ↓
Computes Wavenumber (WN) and Group Velocity (CG)
    ↓
Used in energy conversion (FACTOR = DDEN/CG)
    ↓
Affects spectral moment accumulation
    ↓
Changes final HS, T01, T02 values
```

**Impact**: Same spectrum in different water depths yields different wave parameters!

### 8.2 Frequency Grid Dependency

```
FR1, XFR → Computed in w3gridmd
    ↓
SIG(IK), DSII(IK), DDEN(IK) arrays
    ↓
Grid spacing and integration factors
    ↓
Affects resolution and accuracy
    ↓
Coarser grid → less accurate SWH
```

### 8.3 Tail Factor Dependency

```
SIG(NK), DTH → FTE, FTTR, FTWL computed at init
    ↓
Used to extend spectrum to high frequencies
    ↓
Contributes to M₀, M₁, M₂ moments
    ↓
Critical for accurate HS (5-15% contribution)
```

---

## 9. Physical Meaning of Each Variable Category

### 9.1 Why We Need Grid Parameters
- **Spectral resolution**: More bins = better accuracy but more computation
- **Frequency range**: FR1 and XFR define which waves are resolved
- **Integration factors**: Account for non-uniform log-spacing in frequency

### 9.2 Why We Need Wavenumber and Group Velocity
- **Action to Energy conversion**: Energy = Action × group velocity effects
- **Depth-dependent physics**: Different depths have different wave speeds
- **Jacobian transformation**: Change of variables from k-space to f-space

### 9.3 Why We Need Water Depth
- **Dispersion relation**: Determines how frequency relates to wavenumber
- **Group velocity**: Affects wave propagation and energy conversion
- **Shallow water effects**: Waves behave differently in shallow water

### 9.4 Why We Need Tail Factors
- **High-frequency extension**: Real spectra extend beyond NK bins
- **Energy conservation**: Without tail, 5-15% of energy is missing
- **Extreme event statistics**: Tail affects maximum wave height calculations

---

## 10. Practical Example: Variable Values

```
Typical Configuration:
═════════════════════════════════════════════════════

NK = 30 frequency bins
NTH = 36 directions
FR1 = 0.04 Hz
XFR = 1.1

Grid Parameters:
  SIG(1) = 2π × 0.04 = 0.251 rad/s
  SIG(30) = 2π × 0.04 × 1.1^29 ≈ 3.14 rad/s (f ≈ 0.5 Hz)
  DTH = 2π/36 = 0.1745 rad
  DSII(15) ≈ 0.301 rad/s
  DDEN(15) ≈ 0.1658

For a point at depth h = 50 m:
  WN(15, ISEA) = 0.0625 m^-1 (from WAVNU1)
  CG(15, ISEA) = 1.25 m/s (group velocity)
  FACTOR = 0.1658 / 1.25 = 0.133

If action at (f15, θ10) = 0.002 m²:
  EBD(15, JSEA) contribution = 0.002 × 0.133 = 0.000266 m²

After summing all IK and ITH:
  ET(JSEA) = 2.5 m²
  HS = 4 × √2.5 = 6.32 m
```

---

## 11. Dependency Diagram

```
WW3 SWH Computation Dependencies:
═════════════════════════════════════════════════════════

INPUT SPECTRUM A(θ,f)
    │
    ├─→ WATER DEPTH DW ─→ WAVNU1 ──→ WN(f), CG(f)
    │                                  │
    │   GRID PARAMS                    │
    │   ├─ NK, NTH                    │
    │   ├─ SIG(f)                     │
    │   ├─ DSII(f)                    │
    │   ├─ DDEN(f)                    │
    │   ├─ FTE, FTTR, FTWL            │
    │   └─ DTH                        │
    │                                  │
    └────→ Directional Loop ←──────────┘
             ∑(θ) A(θ,f) → AB(f)

                ↓

           Energy Conversion
           FACTOR = DDEN(f) / CG(f)
           EBD = AB × FACTOR

                ↓

           Spectral Moments
           ET = ∑(f) EBD
           ET1 = ∑(f) EBD × σ
           ET02 = ∑(f) EBD × σ²
           ETR = ∑(f) EBD / σ

                ↓

           Add Tail
           ET += FTE × EBAND
           (etc.)

                ↓

           WAVE PARAMETERS
           ├─ HS = 4√ET
           ├─ T01 = 2π(ET/ET1)
           ├─ T02 = 2π√(ET/ET02)
           ├─ T0M1 = 2π(ETR/ET)
           └─ THM = atan2(ETY, ETX)
```

---

## Summary

To compute SWH from action density spectrum, WW3 requires:

### Must-Have Variables:
1. **A(θ,f)** - Action density spectrum
2. **DW** - Water depth (for dispersion relation)
3. **Grid parameters** - NK, NTH, FR1, XFR, SIG(f), DSII(f)
4. **Physical constants** - GRAV, TPI

### Computed Variables:
1. **WN(f)** - Wavenumber from dispersion relation
2. **CG(f)** - Group velocity (depth-dependent)
3. **Spectral moments** - ET, ET1, ET02, ETR from integration

### Key Insight:
**The same action spectrum produces different SWH values in different water depths!** This is because water depth affects the dispersion relation, which changes the group velocity used in the energy conversion step.

---

**File**: w3iogomd.F90 (W3OUTG subroutine, lines 1197-2379)
**Related**: w3gridmd.F90, w3dispmd.F90, w3adatmd.F90, constants.F90
