# Constants and Variables Classification for SWH Computation
## What Changes During Simulation and What Does Not

---

## **WHAT DOES NOT CHANGE DURING SIMULATION**

### Physical Constants (Never Change)
```
GRAV = 9.81 m/s²           (Gravitational acceleration)
TPI = 2π ≈ 6.283           (2π factor)
TPIINV = 1/(2π) ≈ 0.159    (Inverse 2π)
RADE = 180/π               (Radians to degrees conversion)
```

### Grid Parameters (Fixed at Initialization)
```
Frequency Grid Setup:
├─ NK = 30                              (Number of frequency bins)
├─ NTH = 36                             (Number of directions)
├─ FR1 = 0.04 Hz                        (First frequency)
├─ XFR = 1.1                            (Frequency ratio)

Computed Once at Initialization:
├─ SIG(IK) = 2π×FR1×XFR^(IK-1)         (Angular frequencies - 1D array)
├─ DSII(IK)                             (Frequency bandwidths - 1D array)
├─ DDEN(IK) = DTH×DSII×SIG            (Integration factors - 1D array)
└─ DTH = 2π/NTH                         (Directional step size)

Tail Factors (Computed Once):
├─ FTE = 0.25×SIG(NK)×DTH×SIG(NK)     (Energy tail factor)
├─ FTF = 0.20×DTH×SIG(NK)             (Frequency tail factor)
├─ FTTR = FTF                           (Period tail factor)
└─ FTWL = (GRAV/6)×DTH/SIG(NK)        (Wavelength tail factor)
```

### Environmental Data (Static Bathymetry)
```
Water Depth:
└─ DW(ISEA)                             (Water depth at each location)
   ├─ Read from bathymetry data at init
   ├─ Never modified during simulation
   └─ Static for entire model run

Computed from Water Depth (Once at Init):
├─ WN(IK, ISEA) = Wavenumber array    (2D: NK × NSEA)
│  ├─ Computed via dispersion relation
│  ├─ Depends on: SIG(IK) and DW(ISEA)
│  ├─ Computed in w3initmd.F90:1370-1399
│  └─ Stored and never recalculated
│
└─ CG(IK, ISEA) = Group velocity array (2D: NK × NSEA)
   ├─ Computed via WAVNU1 subroutine
   ├─ Depends on: SIG(IK) and DW(ISEA)
   ├─ Computed in w3initmd.F90:1370-1399
   └─ Stored and never recalculated
```

---

## **WHAT CHANGES DURING SIMULATION**

### Primary Input (Changes Every Time Step)
```
Action Density Spectrum:
└─ A(ITH, IK, JSEA, t)                 (3D array + time dimension)
   ├─ Varies with direction ITH (1 to NTH)
   ├─ Varies with frequency IK (1 to NK)
   ├─ Varies with location JSEA (1 to NSEAL)
   ├─ Varies with time t
   └─ Updated at each model time step
```

### Spectral Moments (Computed Every Time Step)
```
└─ All computed from A(θ,f) and static parameters:

   ET(JSEA, t) = ∑∫ E(f) df              (Total energy, M₀)
   ET1(JSEA, t) = ∑∫ f·E(f) df           (First moment, M₁)
   ET02(JSEA, t) = ∑∫ f²·E(f) df         (Second moment, M₂)
   ETR(JSEA, t) = ∑∫ E(f)/f df           (Inverse moment, M₋₁)
   EWN(JSEA, t) = ∑∫ E(f)/k df           (Wavelength moment)
   ETX(JSEA, t) = ∑∫ E(f,θ)·cos(θ) dθ df (E-W directional moment)
   ETY(JSEA, t) = ∑∫ E(f,θ)·sin(θ) dθ df (N-S directional moment)
```

### Output Wave Parameters (Computed Every Time Step)
```
└─ All computed from spectral moments:

   HS(JSEA, t) = 4√ET                    (Significant Wave Height)
   T01(JSEA, t) = 2π(ET/ET1)             (Mean period)
   T02(JSEA, t) = 2π√(ET/ET02)           (Zero-crossing period)
   T0M1(JSEA, t) = 2π(ETR/ET)            (Energy period)
   THM(JSEA, t) = atan2(ETY, ETX)        (Mean direction)
   WLM(JSEA, t) = 2π(EWN/ET)             (Mean wavelength)
   THS(JSEA, t) = Longuet-Higgins formula (Directional spread)
   QP(JSEA, t) = 2·EET1/ET²              (Peakedness)
```

---

## **CODE VERIFICATION: WN and CG Do Not Change**

### Evidence from Code Analysis

**Location of Computation:**
```
File: w3initmd.F90
Lines: 1370-1399
Section: "5.b Fill wavenumber and group velocity arrays"
Subroutine: W3INIT (initialization phase only)
Timing: ONCE at model startup
```

**Initialization Code (w3initmd.F90:1370-1399):**
```fortran
! 5.b Fill wavenumber and group velocity arrays.
!
DO IS=0, NSEA
  DEPTH = MAX(DMIN, DW(IS))    ← Read water depth from DW array

  DO IK=0, NK+1
    ! Calculate wavenumbers and group velocities.
#ifdef W3_PDLIB
    CALL WAVNU3(SIG(IK), DEPTH, WN(IK,IS), CG(IK,IS))
#else
    CALL WAVNU1(SIG(IK), DEPTH, WN(IK,IS), CG(IK,IS))
#endif
    ! Results stored in WN and CG arrays
  END DO
END DO

! After this, never called again during simulation
```

**Search Results:**

1. **Are WN and CG recomputed during simulation?**
   - Search in w3wavemd.F90 (main simulation loop)
   - Result: ❌ NO `CALL WAVNU` statements found
   - Conclusion: Never recalculated during simulation

2. **Are there assignments to WN or CG during simulation?**
   - Search for: `WN(.*) =` in w3wavemd.F90
   - Result: ❌ NO assignments found
   - Search for: `CG(.*) =` in w3wavemd.F90
   - Result: ❌ NO assignments found
   - Conclusion: Never modified during simulation

3. **How are WN and CG used during simulation?**
   - File: w3iogomd.F90 (output calculations, called every time step)
   - Usage examples:
     * Line 1529: `CG(IK,ISEA) / SIG(IK) * WN(IK,ISEA)` ← READ only
     * Line 1552: `FACTOR = DDEN(IK) / CG(IK,ISEA)` ← READ only
     * Line 1558: `EBD(IK,JSEA) * CG(IK,ISEA)` ← READ only
     * Line 1598: `TPI*SIG(IK) / CG(IK,ISEA)` ← READ only
   - Operation: READING VALUES ONLY, never modified

4. **Is water depth DW static?**
   - Search in w3wavemd.F90 for assignments: `DW(ISEA) =`
   - Result: ❌ NO assignments found during simulation
   - How DW is used:
     * Line 972: `DEPTH = MAX(DMIN, DW(ISEA))` ← READ only
     * Line 1763: `DEPTH = MAX(DMIN, DW(ISEA))` ← READ only
     * Line 2075: `DEPTH = MAX(DMIN, DW(ISEA))` ← READ only
   - Conclusion: DW is READ-ONLY during simulation

---

## **Dependency Chain**

### At Initialization (One Time Only)
```
DW(ISEA) ──────────→ WAVNU1 ──────────→ WN(IK,ISEA), CG(IK,ISEA)
SIG(IK)  ──────────→ Subroutine ───────→ Stored in 2D arrays
                     (line 1391)        (Never recalculated)
```

### At Each Time Step (Repeatedly)
```
A(θ,f,t) ──────────→ Integration Loop
                     ├─ Uses fixed WN(IK,ISEA)
                     ├─ Uses fixed CG(IK,ISEA)
                     ├─ Uses fixed DDEN(IK)
                     ├─ Uses fixed SIG(IK)
                     ↓
                     Computes: ET, ET1, ET02, ETR
                     ↓
                     Computes: HS = 4√ET, T01, T02, etc.
                     ↓
                     Output stored for this time step
```

---

## **Complete Classification Table**

| Variable | Type | Varies With | Dimension | Status | When Set |
|----------|------|-------------|-----------|--------|----------|
| **GRAV** | Constant | — | Scalar | Static | Compile-time |
| **TPI** | Constant | — | Scalar | Static | Compile-time |
| **NK** | Int | — | Scalar | Static | Grid init |
| **NTH** | Int | — | Scalar | Static | Grid init |
| **SIG(IK)** | Array | Frequency bin | 1D (NK) | Static | Grid init |
| **DSII(IK)** | Array | Frequency bin | 1D (NK) | Static | Grid init |
| **DDEN(IK)** | Array | Frequency bin | 1D (NK) | Static | Grid init |
| **DTH** | Scalar | — | Scalar | Static | Grid init |
| **FTE, FTTR, FTWL** | Scalar | — | Scalar | Static | Grid init |
| **DW(ISEA)** | Array | Location | 1D (NSEA) | Static | Bathymetry init |
| **WN(IK,ISEA)** | Array | Freq + Loc | 2D (NK×NSEA) | Static | Model init |
| **CG(IK,ISEA)** | Array | Freq + Loc | 2D (NK×NSEA) | Static | Model init |
| **A(θ,f,JSEA)** | Array | All+Time | 3D+T | **Dynamic** | Every timestep |
| **ET(JSEA)** | Array | Loc+Time | 1D+T | **Dynamic** | Every timestep |
| **ET1(JSEA)** | Array | Loc+Time | 1D+T | **Dynamic** | Every timestep |
| **ET02(JSEA)** | Array | Loc+Time | 1D+T | **Dynamic** | Every timestep |
| **HS(JSEA)** | Array | Loc+Time | 1D+T | **Dynamic** | Every timestep |
| **T01(JSEA)** | Array | Loc+Time | 1D+T | **Dynamic** | Every timestep |
| **T02(JSEA)** | Array | Loc+Time | 1D+T | **Dynamic** | Every timestep |

---

## **Summary**

### ✅ Static During Simulation
- Physical constants (GRAV, TPI, etc.)
- Grid parameters (NK, NTH, SIG, DSII, DDEN)
- Tail factors (FTE, FTTR, FTWL)
- Water depth (DW)
- **Wavenumber (WN) - computed once from DW**
- **Group velocity (CG) - computed once from DW**

### 🔄 Dynamic During Simulation
- Action density spectrum (A)
- Spectral moments (ET, ET1, ET02, ETR)
- Wave parameters (HS, T01, T02, T0M1)
- Directional moments (ETX, ETY)
- All other derived outputs

### Key Insight
**WN and CG are computed ONCE at model initialization from the static water depth DW and frequency grid SIG, then stored in memory and never recalculated. They are used as lookup tables during the simulation.**

---

## **Files Referenced**

| File | Purpose | Line Range |
|------|---------|-----------|
| w3initmd.F90 | WN/CG initialization | 1370-1399 |
| w3dispmd.F90 | WAVNU1 subroutine | Dispersion solver |
| w3wavemd.F90 | Main simulation loop | Reads but never modifies |
| w3iogomd.F90 | Output calculation | Uses WN/CG as read-only |
| w3gdatmd.F90 | Data structures | WN/CG definitions |
| constants.F90 | Physical constants | GRAV, TPI, etc. |
