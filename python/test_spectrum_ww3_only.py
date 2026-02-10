#!/usr/bin/env python
"""
Test script for WW3-only spectrum module.

Demonstrates the streamlined WaveSpectrum class that only accepts
pre-computed WW3 grid parameters.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/user/WW3/python')

from wavewatch_python.spectrum_ww3_only import WaveSpectrum
from wavewatch_python.dispersion import solve_dispersion


def test_ww3_only_spectrum():
    """Test WaveSpectrum with WW3 pre-computed parameters."""

    print("="*70)
    print("Testing WW3-Only WaveSpectrum Class")
    print("="*70)

    # ===== WW3 Grid Parameters (Pre-computed) =====
    nfreq = 30
    ndir = 36
    depth = 100.0

    # Frequency grid from WW3 (FR1=0.04, XFR=1.1)
    fr1 = 0.04
    xfr = 1.1
    omega_ww3 = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)

    # Compute DDEN = DTH × DSII × SIG
    dth = 2.0 * np.pi / ndir
    dsii = np.zeros(nfreq)
    dsii[0] = (omega_ww3[1] - omega_ww3[0]) / 2.0
    for i in range(1, nfreq - 1):
        dsii[i] = (omega_ww3[i+1] - omega_ww3[i-1]) / 2.0
    dsii[-1] = (omega_ww3[-1] - omega_ww3[-2]) / 2.0
    dden_ww3 = dth * dsii * omega_ww3

    # Tail factors (from w3gridmd.F90:3444-3448)
    fte_ww3 = 0.25 * omega_ww3[-1] * dth * omega_ww3[-1]
    fttr_ww3 = 0.20 * dth * omega_ww3[-1]
    ftwl_ww3 = (9.81 / 6.0) / omega_ww3[-1] * dth * omega_ww3[-1]

    # Wavenumber and group velocity (from WAVNU1)
    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # ===== Create test action density spectrum =====
    # JONSWAP-like spectrum
    action_ww3 = np.zeros((ndir, nfreq))
    for ifreq in range(nfreq):
        # JONSWAP-like energy distribution
        f = omega_ww3[ifreq] / (2.0 * np.pi)
        fp = 0.1  # Peak frequency
        gam = 3.3
        alpha = 0.0081

        e_f = (alpha * 9.81**2 / (2*np.pi)**4 * f**(-5) *
               np.exp(-1.25 * (fp/f)**4) * gam**np.exp(-(f-fp)**2 / (2*fp**2/gam)))

        # Directional spreading
        s = 2.0
        for idir in range(ndir):
            dir_spread = np.cos(np.linspace(0, 2*np.pi, ndir)[idir])**(2*s)
            action_ww3[idir, ifreq] = e_f / cg_ww3[ifreq] * dir_spread

    # ===== Create WaveSpectrum with WW3 parameters =====
    print("\n1. Creating WaveSpectrum with WW3 pre-computed parameters...")
    spectrum = WaveSpectrum(
        action=action_ww3,
        depth=depth,
        omega=omega_ww3,
        dintegral=dden_ww3,
        fte=fte_ww3,
        fttr=fttr_ww3,
        ftwl=ftwl_ww3,
        wavenumber=wn_ww3,
        group_velocity=cg_ww3
    )

    print(f"   ✓ Spectrum created successfully")
    print(f"   - Grid dimensions: {ndir} directions × {nfreq} frequencies")
    print(f"   - Frequency range: {spectrum.frequencies[0]:.4f} - {spectrum.frequencies[-1]:.4f} Hz")
    print(f"   - Depth: {depth} m")

    # ===== Compute wave parameters =====
    print("\n2. Computing wave parameters...")
    params = spectrum.compute_parameters()

    print(f"   ✓ Parameters computed:")
    print(f"   - HS: {params['hs']:.3f} m")
    print(f"   - T01: {params['t01']:.2f} s")
    print(f"   - T02: {params['t02']:.2f} s")
    print(f"   - T0M1: {params['t0m1']:.2f} s")
    print(f"   - TP: {params['tp']:.2f} s")
    print(f"   - Mean Direction: {np.degrees(params['thm']):.1f}°")
    print(f"   - Directional Spread: {np.degrees(params['ths']):.1f}°")
    print(f"   - Mean Wavelength: {params['wlm']:.2f} m")
    print(f"   - Spectral Width: {params['width']:.4f}")

    # ===== Check moments =====
    print("\n3. Spectral moments:")
    moments = params['moments']
    print(f"   - M0: {moments['m0']:.6f} m²")
    print(f"   - M1: {moments['m1']:.6f} m²·s")
    print(f"   - M2: {moments['m2']:.6f} m²·s²")
    print(f"   - M-1: {moments['m_1']:.6f} m²·s")

    # ===== Extract 1D spectra =====
    print("\n4. Extracting 1D spectra...")
    freq, e_freq = spectrum.get_spectrum_1d()
    dirs, e_dir = spectrum.get_spectrum_directional_1d()

    print(f"   ✓ 1D spectra extracted:")
    print(f"   - Frequency spectrum shape: {e_freq.shape}")
    print(f"   - Directional spectrum shape: {e_dir.shape}")
    print(f"   - Freq range: {freq[0]:.4f} - {freq[-1]:.4f} Hz")
    print(f"   - Dir range: {np.degrees(dirs[0]):.1f}° - {np.degrees(dirs[-1]):.1f}°")

    # ===== Test with different depths =====
    print("\n5. Testing depth dependency...")
    depths = [5.0, 20.0, 100.0, 1000.0]

    print(f"   {'Depth (m)':>10} | {'HS (m)':>8} | {'T01 (s)':>8} | {'T02 (s)':>8} | {'CG[0] (m/s)':>10}")
    print("   " + "-"*60)

    for test_depth in depths:
        # Recompute wavenumber and group velocity for different depth
        wn_test = np.zeros(nfreq)
        cg_test = np.zeros(nfreq)
        for ik in range(nfreq):
            wn_test[ik], cg_test[ik] = solve_dispersion(omega_ww3[ik], test_depth)

        # Create new spectrum with different depth
        spectrum_test = WaveSpectrum(
            action=action_ww3,
            depth=test_depth,
            omega=omega_ww3,
            dintegral=dden_ww3,
            fte=fte_ww3,
            fttr=fttr_ww3,
            ftwl=ftwl_ww3,
            wavenumber=wn_test,
            group_velocity=cg_test
        )

        params_test = spectrum_test.compute_parameters()
        print(f"   {test_depth:>10.1f} | {params_test['hs']:>8.3f} | "
              f"{params_test['t01']:>8.2f} | {params_test['t02']:>8.2f} | {cg_test[0]:>10.3f}")

    print("\n" + "="*70)
    print("All tests completed successfully!")
    print("="*70)


if __name__ == "__main__":
    test_ww3_only_spectrum()
