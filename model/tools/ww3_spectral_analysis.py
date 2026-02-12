"""
WW3 Spectral Analysis: Action Density to Wave Parameters
=========================================================

Converts WW3 action density N(k, theta) to:
  - Wave energy spectrum E(f, theta)    [m^2 / Hz / rad]
  - Frequency spectrum E(f)             [m^2 / Hz]
  - Directional spectrum E(theta)       [m^2 / rad]
  - Significant wave height Hs          [m]
  - Mean wave parameters (Tm01, Tm02, Tm-10, fp, theta_m, sigma_theta, ...)

All transformations mirror the Fortran code in w3iogomd.F90 (W3OUTG subroutine)
and w3gridmd.F90 (grid setup).

Assumes the input comes from a WW3 simulation where:
  - A(NTH, NK)  is the action density N(k, theta) stored on a (k, theta) grid
  - SIG(NK)     is the angular frequency array  sigma = 2*pi*f  [rad/s]
  - WN(NK)      is the wavenumber array  k  [rad/m]
  - CG(NK)      is the group velocity array  [m/s]
  - NTH         is the number of directional bins
  - NK          is the number of frequency bins
  - XFR         is the frequency increment factor  (f_{k+1} / f_k)
  - FR1         is the first frequency  [Hz]
  - DTH         is the directional bin width  = 2*pi / NTH  [rad]
"""

import numpy as np


# ---------------------------------------------------------------------------
# 1. Build WW3 spectral grid (mirrors w3gridmd.F90 lines 1297-1316)
# ---------------------------------------------------------------------------

def build_ww3_grid(fr1, xfr, nk, nth):
    """
    Construct the WW3 spectral grid arrays exactly as in w3gridmd.F90.

    Parameters
    ----------
    fr1 : float
        First discrete frequency [Hz].
    xfr : float
        Frequency increment factor (geometric ratio).
    nk : int
        Number of frequency bins.
    nth : int
        Number of directional bins.

    Returns
    -------
    grid : dict with keys
        'freq'   : (NK,)   frequencies f [Hz]
        'sigma'  : (NK,)   angular frequencies sigma = 2*pi*f [rad/s]
        'theta'  : (NTH,)  direction centres [rad]
        'dth'    : float   directional bin width [rad]
        'dsii'   : (NK,)   sigma bandwidths [rad/s]
        'dden'   : (NK,)   DTH * DSII * SIG  (WW3 conversion factor)
        'dsip'   : (NK,)   sigma bandwidth (interior formula)
        'fte'    : float   tail energy factor
    """
    tpi = 2.0 * np.pi
    dth = tpi / nth

    # --- Directions (evenly spaced, centred on 0) ---
    theta = np.array([(ith - 0.5) * dth for ith in range(1, nth + 1)])

    # --- Frequencies (geometric series, mirrors Fortran IK=1..NK) ---
    # In WW3 the loop starts from IK=0 with SIGMA = FR1*TPI/XFR**2
    # then multiplies by XFR each step.  IK=1 gives SIGMA = FR1*TPI/XFR,
    # but the standard convention gives SIG(1) = FR1*TPI.
    # Following the Fortran exactly:
    sigma_start = fr1 * tpi / xfr**2
    sigma = np.zeros(nk + 2)  # indices 0..NK+1 in Fortran
    for ik in range(nk + 2):
        sigma_start *= xfr
        sigma[ik] = sigma_start

    # Trim to IK=1..NK  (Fortran 1-based → Python 0-based)
    sig = sigma[1:nk + 1]    # SIG(1:NK) in Fortran
    freq = sig / tpi

    # --- Sigma bandwidths (DSIP, DSII) ---
    sxfr = 0.5 * (xfr - 1.0 / xfr)
    dsip = sig * sxfr

    dsii = np.copy(dsip)
    dsii[0] = 0.5 * sig[0] * (xfr - 1.0)           # first bin
    dsii[-1] = 0.5 * sig[-1] * (xfr - 1.0) / xfr   # last bin

    # --- DDEN: combined conversion factor ---
    dden = dth * dsii * sig

    # --- Tail energy factor (w3gridmd.F90 line 3444) ---
    fte = 0.25 * sig[-1] * dth * sig[-1]

    return {
        'freq': freq,
        'sigma': sig,
        'theta': theta,
        'dth': dth,
        'dsii': dsii,
        'dsip': dsip,
        'dden': dden,
        'fte': fte,
    }


# ---------------------------------------------------------------------------
# 2. Action density → Energy spectra  (mirrors w3iogomd.F90 lines 1484-1663)
# ---------------------------------------------------------------------------

def action_to_energy(action, cg, grid):
    """
    Convert WW3 action density N(k, theta) to energy spectrum E(f, theta).

    This is the core Jacobian transformation:
        E(f, theta) = N(k, theta) * sigma / Cg * (dsigma/df)
    with dsigma/df = 2*pi.

    Parameters
    ----------
    action : ndarray, shape (NTH, NK)
        Action density A(ITH, IK) as stored by WW3 (N(k,theta) on the
        native (k, theta) grid).
    cg : ndarray, shape (NK,)
        Group velocity [m/s] at each frequency bin (depth-dependent).
    grid : dict
        Output of build_ww3_grid().

    Returns
    -------
    result : dict with keys
        'Ef_theta' : (NTH, NK) E(f, theta) in [m^2 / Hz / rad]
        'Ef'       : (NK,)     E(f)         in [m^2 / Hz]
        'Etheta'   : (NTH,)    E(theta)     in [m^2 / rad]
        'EBD'      : (NK,)     E(f)*df      in [m^2]  (energy per freq band)
    """
    tpi = 2.0 * np.pi
    nth, nk = action.shape
    sig = grid['sigma']
    dth = grid['dth']
    dsii = grid['dsii']
    dden = grid['dden']

    # ---- Step 1:  Direction-integrated action per frequency band ----
    #   AB(IK) = sum_theta A(ITH, IK)         [line 1513]
    AB = np.sum(action, axis=0)  # shape (NK,)

    # ---- Step 2:  Energy per frequency band  EBD = E(f)*df  ----
    #   FACTOR = DDEN(IK) / CG(IK)            [line 1552]
    #   EBD(IK) = AB(IK) * FACTOR             [line 1553]
    #   where DDEN = DTH * DSII * SIG
    factor = dden / cg                         # shape (NK,)
    EBD = AB * factor                          # E(f)*df  [m^2]

    # ---- Step 3:  E(f)  spectral density  ----
    #   EBD(IK) = EBD(IK) / DSII(IK)   →  E(sigma)    [line 1660]
    #   EF(IK)  = EBD(IK) * TPI        →  E(f)         [line 1663]
    E_sigma = EBD / dsii                       # E(sigma) [m^2 s/rad]
    Ef = E_sigma * tpi                         # E(f)     [m^2/Hz]

    # ---- Step 4:  Full 2-D energy spectrum E(f, theta)  ----
    # Per-bin, per-direction:
    #   E(sigma, theta) = A(ITH, IK) * SIG(IK) / CG(IK)      (omit DTH,DSII)
    #   E(f, theta)     = E(sigma, theta) * 2*pi
    Ef_theta = np.zeros_like(action)
    for ik in range(nk):
        Ef_theta[:, ik] = action[:, ik] * sig[ik] / cg[ik] * tpi

    # ---- Step 5:  E(theta)  direction spectrum ----
    # Integrate E(f, theta) over frequency:
    #   E(theta_j) = sum_IK  A(ITH_j, IK) * DSII(IK) * SIG(IK) / CG(IK)
    # which is sum of E(sigma, theta)*dsigma.  Convert to per-radian:
    Etheta = np.zeros(nth)
    for ith in range(nth):
        for ik in range(nk):
            Etheta[ith] += action[ith, ik] * dsii[ik] * sig[ik] / cg[ik]

    return {
        'Ef_theta': Ef_theta,
        'Ef': Ef,
        'Etheta': Etheta,
        'EBD': EBD,
    }


# ---------------------------------------------------------------------------
# 3. Mean wave parameters  (mirrors w3iogomd.F90 lines 1550-2050)
# ---------------------------------------------------------------------------

def compute_mean_params(action, cg, wn, grid):
    """
    Compute significant wave height and bulk mean parameters from action
    density, following the WW3 w3iogomd.F90 W3OUTG subroutine.

    Parameters
    ----------
    action : ndarray, shape (NTH, NK)
        Action density on the WW3 (k, theta) grid.
    cg : ndarray, shape (NK,)
        Group velocity [m/s].
    wn : ndarray, shape (NK,)
        Wavenumber [rad/m].
    grid : dict
        Output of build_ww3_grid().

    Returns
    -------
    params : dict with all computed parameters.
    spectra : dict with spectral arrays (from action_to_energy).
    """
    tpi = 2.0 * np.pi
    grav = 9.81
    rade = 180.0 / np.pi

    nth, nk = action.shape
    sig = grid['sigma']
    dth = grid['dth']
    dsii = grid['dsii']
    dden = grid['dden']
    fte = grid['fte']
    theta = grid['theta']

    ecos = np.cos(theta)
    esin = np.sin(theta)

    # --- Energy spectra ---
    spectra = action_to_energy(action, cg, grid)
    EBD = spectra['EBD']  # E(f)*df per band

    # -----------------------------------------------------------------
    # Accumulate spectral moments over discrete bands (lines 1484-1685)
    # -----------------------------------------------------------------
    ET = 0.0     # m0   = integral E df
    ET1 = 0.0    # m1   = integral E * sigma df
    ET02 = 0.0   # m2   = integral E * sigma^2 df
    ETR = 0.0    # m-1  = integral E / sigma df
    EWN = 0.0    # integral E / k df
    ETF = 0.0    # integral Cg * E df  (energy flux)
    EET1 = 0.0   # for Qp = (2/m0^2) * integral (E(f)*df)^2 * sigma / dsigma

    ETX = 0.0    # integral E * cos(theta) df  (x-component of mean dir)
    ETY = 0.0    # integral E * sin(theta) df  (y-component of mean dir)

    for ik in range(nk):
        factor = dden[ik] / cg[ik]

        # Direction-integrated sums for this frequency band
        AB = np.sum(action[:, ik])
        ABX = np.sum(action[:, ik] * ecos)
        ABY = np.sum(action[:, ik] * esin)

        ebd_ik = AB * factor      # E(f)*df

        # Spectral moments  [lines 1554-1564]
        ET += ebd_ik
        ET1 += ebd_ik * sig[ik]
        ET02 += ebd_ik * sig[ik]**2
        ETR += ebd_ik / sig[ik]
        EWN += ebd_ik / wn[ik]
        ETF += ebd_ik * cg[ik]
        EET1 += ebd_ik**2 * sig[ik] / dsii[ik]

        # Directional moments  [lines 1565-1566]
        ETX += ABX * factor
        ETY += ABY * factor

    # -----------------------------------------------------------------
    # Add diagnostic tail  (lines 1956-1974)
    # -----------------------------------------------------------------
    # EBAND = AB_last / CG(NK)  where AB_last = sum_theta A(ITH, NK)
    AB_last = np.sum(action[:, -1])
    ABX_last = np.sum(action[:, -1] * ecos)
    ABY_last = np.sum(action[:, -1] * esin)
    EBAND = AB_last / cg[-1]

    ET += fte * EBAND
    ETR += fte * EBAND    # FTR = FTE for standard WW3
    ET1 += fte * EBAND    # approximate (FT1 ≈ FTE for standard tail)
    ET02 += EBAND * 0.5 * sig[-1]**4 * dth
    ETX += fte * ABX_last / cg[-1]
    ETY += fte * ABY_last / cg[-1]
    ETF += grav * fte * EBAND   # deep water: Cg ≈ g/(2*sigma)

    # -----------------------------------------------------------------
    # Compute bulk parameters  (lines 1989-2050)
    # -----------------------------------------------------------------

    # Significant wave height  [line 1997]
    Hs = 4.0 * np.sqrt(max(0.0, ET))

    # Mean periods  [lines 2006, 2039-2040]
    if ET > 1.0e-7:
        Tm_10 = ETR / ET * tpi            # T_{-1,0} = m_{-1}/m_0 * 2pi
        Tm01 = ET / ET1 * tpi             # T_{0,1}  = m_0/m_1 * 2pi
        Tm02 = tpi * np.sqrt(ET / ET02)   # T_{0,2}  = 2pi * sqrt(m0/m2)
    else:
        Tm_10 = 0.0
        Tm01 = 0.0
        Tm02 = 0.0

    # Mean wavelength  [line 2005]
    WLM = EWN / ET * tpi if ET > 1.0e-7 else 0.0

    # Mean direction  [line 2019]
    if abs(ETX) + abs(ETY) > 1.0e-7:
        theta_m = np.arctan2(ETY, ETX)
    else:
        theta_m = 0.0

    # Directional spread  [lines 2007-2009]
    if ET > 1.0e-7:
        m1 = np.sqrt(ETX**2 + ETY**2) / ET
        sigma_theta = np.sqrt(max(0.0, 2.0 * (1.0 - m1)))   # [rad]
    else:
        sigma_theta = 0.0

    # Peakedness parameter Qp  [line 2004]
    Qp = (2.0 / ET**2) * EET1 if ET > 1.0e-7 else 0.0

    # Energy flux  [line 2037]
    CGE = 1025.0 * grav * ETF   # rho_w * g * integral(Cg*E)  [W/m]

    # -----------------------------------------------------------------
    # Peak frequency  (lines 2077-2113)
    # Find IK where E(f) = EBD/DSII is maximum, then fp = sig(IKp)/(2pi)
    # -----------------------------------------------------------------
    E_sigma = EBD / dsii
    ikp = np.argmax(E_sigma)
    fp = sig[ikp] / tpi
    Tp = 1.0 / fp if fp > 0 else 0.0

    # Peak direction: direction with max energy at peak frequency
    theta_p = theta[np.argmax(action[:, ikp])]

    params = {
        'Hs': Hs,
        'Tp': Tp,
        'fp': fp,
        'Tm01': Tm01,
        'Tm02': Tm02,
        'Tm_10': Tm_10,
        'theta_m': theta_m,
        'theta_m_deg': np.degrees(theta_m),
        'theta_p': theta_p,
        'theta_p_deg': np.degrees(theta_p),
        'sigma_theta': sigma_theta,
        'sigma_theta_deg': np.degrees(sigma_theta),
        'Qp': Qp,
        'WLM': WLM,
        'CGE': CGE,
        'm0': ET,
        'm1': ET1,
        'm2': ET02,
        'm_1': ETR,
    }

    return params, spectra


# ---------------------------------------------------------------------------
# 4. Demonstration with synthetic WW3-like data
# ---------------------------------------------------------------------------

def demo():
    """
    Demonstrate the full chain N(k,theta) → E(f,theta) → E(f) → E(theta) → Hs
    using synthetic data that mimics a typical WW3 output.
    """
    # --- WW3 grid parameters (typical operational settings) ---
    fr1 = 0.0373        # first frequency [Hz]
    xfr = 1.1           # geometric ratio
    nk = 32             # number of frequency bins
    nth = 36            # number of direction bins (10-degree resolution)

    grid = build_ww3_grid(fr1, xfr, nk, nth)
    freq = grid['freq']
    sigma = grid['sigma']
    theta = grid['theta']
    dth = grid['dth']

    # --- Deep-water dispersion: k = sigma^2 / g,  Cg = g / (2*sigma) ---
    grav = 9.81
    wn = sigma**2 / grav
    cg = grav / (2.0 * sigma)

    # --- Build a synthetic action density A(NTH, NK) ---
    # Create a JONSWAP-like shape with cos^2s directional spreading
    fp_true = 0.08                      # peak frequency [Hz]
    sigp = fp_true * 2.0 * np.pi
    alpha = 0.01
    gamma = 3.3
    theta_peak = np.radians(210.0)      # waves coming from SW

    # JONSWAP E(f) on the WW3 frequency grid
    Ef_jonswap = np.zeros(nk)
    for ik in range(nk):
        f = freq[ik]
        sig_j = 0.07 if f <= fp_true else 0.09
        Ef_jonswap[ik] = (alpha * grav**2 / (2.0 * np.pi)**4 / f**5
                          * np.exp(-1.25 * (fp_true / f)**4)
                          * gamma**np.exp(-0.5 * ((f / fp_true - 1.0) / sig_j)**2))

    # cos^2s directional distribution  D(theta)
    s = 10.0
    D = np.cos(0.5 * (theta - theta_peak))**( 2.0 * s)
    D /= (np.sum(D) * dth)  # normalise so integral D dtheta = 1

    # E(f, theta) = E(f) * D(theta)
    Ef_theta_true = np.outer(D, Ef_jonswap)   # (NTH, NK)

    # Convert E(f,theta) → N(k,theta) to create our "WW3 input"
    # E(f,theta) = A(ITH,IK) * sigma / Cg * 2*pi   →
    # A(ITH,IK) = E(f,theta) * Cg / (sigma * 2*pi)
    action = np.zeros((nth, nk))
    for ik in range(nk):
        action[:, ik] = Ef_theta_true[:, ik] * cg[ik] / (sigma[ik] * 2.0 * np.pi)

    # ===================================================================
    # Run the full conversion chain
    # ===================================================================
    params, spectra = compute_mean_params(action, cg, wn, grid)

    # --- Print results ---
    print("=" * 65)
    print("  WW3 Spectral Analysis: N(k,theta) → E(f,theta) → Wave Params")
    print("=" * 65)
    print()
    print("Grid setup:")
    print(f"  Frequencies : NK={nk}, f1={freq[0]:.4f} Hz, fN={freq[-1]:.4f} Hz")
    print(f"  Directions  : NTH={nth}, dtheta={np.degrees(dth):.1f} deg")
    print()
    print("--- Bulk Wave Parameters ---")
    print(f"  Hs            = {params['Hs']:.3f} m")
    print(f"  Tp            = {params['Tp']:.2f} s    (fp = {params['fp']:.4f} Hz)")
    print(f"  Tm01          = {params['Tm01']:.2f} s")
    print(f"  Tm02          = {params['Tm02']:.2f} s")
    print(f"  Tm-10         = {params['Tm_10']:.2f} s")
    print(f"  Mean dir      = {params['theta_m_deg']:.1f} deg")
    print(f"  Peak dir      = {params['theta_p_deg']:.1f} deg")
    print(f"  Dir spread    = {params['sigma_theta_deg']:.1f} deg")
    print(f"  Mean wvlen    = {params['WLM']:.1f} m")
    print(f"  Peakedness Qp = {params['Qp']:.2f}")
    print(f"  Energy flux   = {params['CGE']:.1f} W/m")
    print()
    print("--- Spectral Moments ---")
    print(f"  m0   = {params['m0']:.6f} m^2")
    print(f"  m-1  = {params['m_1']:.6f} m^2 s")
    print(f"  m1   = {params['m1']:.6f} m^2/s")
    print(f"  m2   = {params['m2']:.6f} m^2/s^2")
    print()

    # --- Verify: Hs = 4*sqrt(m0) ---
    Hs_check = 4.0 * np.sqrt(params['m0'])
    print(f"  Hs check: 4*sqrt(m0) = {Hs_check:.3f} m  (matches Hs = {params['Hs']:.3f} m)")
    print()

    # --- Verify: integral of E(f)*df ≈ m0 ---
    Ef = spectra['Ef']
    m0_from_Ef = np.sum(Ef * grid['dsii'] / (2.0 * np.pi))
    # Equivalently from EBD:
    m0_from_EBD = np.sum(spectra['EBD'])
    print(f"  m0 from sum(EBD)          = {m0_from_EBD:.6f} m^2")
    print(f"  m0 from sum(E(f)*df)      = {m0_from_Ef:.6f} m^2")
    print()

    # --- Verify: integral of E(theta)*dtheta ≈ m0 ---
    Etheta = spectra['Etheta']
    m0_from_Etheta = np.sum(Etheta) * dth
    print(f"  m0 from sum(E(theta)*dth) = {m0_from_Etheta:.6f} m^2")
    print()

    # --- Verify: integral of E(f,theta)*df*dtheta ≈ m0 ---
    Ef_theta = spectra['Ef_theta']
    m0_from_2d = 0.0
    for ik in range(nk):
        m0_from_2d += np.sum(Ef_theta[:, ik]) * dth * grid['dsii'][ik] / (2.0 * np.pi)
    print(f"  m0 from E(f,theta)*df*dth = {m0_from_2d:.6f} m^2")
    print()
    print("=" * 65)

    return params, spectra, grid


if __name__ == "__main__":
    demo()
