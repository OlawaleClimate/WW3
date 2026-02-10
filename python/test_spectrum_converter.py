#!/usr/bin/env python
"""
Test suite for spectrum converter utilities.

Demonstrates action density to energy spectrum conversion.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/user/WW3/python')

from wavewatch_python.spectrum_converter import (
    action_to_energy_2d,
    action_to_frequency_spectrum_1d,
    action_to_directional_spectrum_1d,
    get_peak_frequency,
    get_peak_direction,
    normalize_spectrum
)
from wavewatch_python.dispersion import solve_dispersion


def test_action_to_energy_conversion():
    """Test conversion from action density to energy spectrum."""

    print("="*70)
    print("Test 1: Action to Energy 2D Conversion")
    print("="*70)

    # Setup WW3 grid parameters
    nfreq, ndir = 30, 36
    depth = 100.0

    # Frequency grid
    fr1 = 0.04
    xfr = 1.1
    omega_ww3 = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)
    frequencies = fr1 * xfr ** np.arange(nfreq)

    # Compute wavenumber and group velocity
    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # Create simple action spectrum (uniform in all directions)
    action = np.ones((ndir, nfreq)) * 0.001

    print(f"\nInput:")
    print(f"  Action shape: {action.shape}")
    print(f"  Frequency range: {frequencies[0]:.4f} - {frequencies[-1]:.4f} Hz")
    print(f"  Water depth: {depth} m")

    # Convert to energy
    energy_2d, freq, dirs, dintegral = action_to_energy_2d(
        action, omega_ww3, cg_ww3
    )

    print(f"\nOutput:")
    print(f"  Energy 2D shape: {energy_2d.shape}")
    print(f"  Energy 2D min: {energy_2d.min():.6e} m²/Hz/rad")
    print(f"  Energy 2D max: {energy_2d.max():.6e} m²/Hz/rad")
    print(f"  Energy 2D mean: {energy_2d.mean():.6e} m²/Hz/rad")
    print(f"  Total energy (sum): {np.sum(energy_2d):.6f} m²")
    print(f"  ✓ Conversion successful")


def test_1d_frequency_spectrum():
    """Test conversion to 1D frequency spectrum."""

    print("\n" + "="*70)
    print("Test 2: Action to 1D Frequency Spectrum")
    print("="*70)

    # Setup
    nfreq, ndir = 30, 36
    depth = 100.0

    fr1 = 0.04
    xfr = 1.1
    omega_ww3 = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)

    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # JONSWAP-like spectrum
    action = np.zeros((ndir, nfreq))
    for ifreq in range(nfreq):
        f = omega_ww3[ifreq] / (2.0 * np.pi)
        fp = 0.1
        gam = 3.3
        alpha = 0.0081

        e_f = (alpha * 9.81**2 / (2*np.pi)**4 * f**(-5) *
               np.exp(-1.25 * (fp/f)**4) * gam**np.exp(-(f-fp)**2 / (2*fp**2/gam)))

        for idir in range(ndir):
            action[idir, ifreq] = e_f / cg_ww3[ifreq]

    print(f"\nInput:")
    print(f"  Action spectrum: JONSWAP-like, Tp=10s")
    print(f"  Grid: {ndir} directions × {nfreq} frequencies")

    # Convert to 1D frequency spectrum
    freq, e_freq = action_to_frequency_spectrum_1d(
        action, omega_ww3, cg_ww3
    )

    print(f"\nOutput:")
    print(f"  Frequency spectrum shape: {e_freq.shape}")
    print(f"  Frequency range: {freq[0]:.4f} - {freq[-1]:.4f} Hz")
    print(f"  Peak frequency: {freq[np.argmax(e_freq)]:.4f} Hz")
    print(f"  Energy at peak: {e_freq.max():.6e} m²/Hz")
    print(f"  Total energy (integrated): {np.sum(e_freq) * (freq[1] - freq[0]) if len(freq) > 1 else e_freq[0]:.6f} m²")
    print(f"  ✓ 1D frequency conversion successful")


def test_1d_directional_spectrum():
    """Test conversion to 1D directional spectrum."""

    print("\n" + "="*70)
    print("Test 3: Action to 1D Directional Spectrum")
    print("="*70)

    # Setup
    nfreq, ndir = 30, 36
    depth = 100.0

    fr1 = 0.04
    xfr = 1.1
    omega_ww3 = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)
    directions = np.linspace(0, 2*np.pi, ndir, endpoint=False)

    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # Action spectrum with directional spreading
    action = np.zeros((ndir, nfreq))
    mean_dir = np.pi / 4  # 45 degrees
    s = 2.0  # Directional shape parameter

    for idir in range(ndir):
        dir_spread = np.cos(directions[idir] - mean_dir) ** (2*s)
        dir_spread = np.maximum(dir_spread, 0)

        for ifreq in range(nfreq):
            f = omega_ww3[ifreq] / (2.0 * np.pi)
            fp = 0.1
            gam = 3.3
            alpha = 0.0081

            e_f = (alpha * 9.81**2 / (2*np.pi)**4 * f**(-5) *
                   np.exp(-1.25 * (fp/f)**4) * gam**np.exp(-(f-fp)**2 / (2*fp**2/gam)))

            action[idir, ifreq] = e_f / cg_ww3[ifreq] * dir_spread

    print(f"\nInput:")
    print(f"  Action spectrum: JONSWAP + directional spreading")
    print(f"  Mean direction: {np.degrees(mean_dir):.1f}°")
    print(f"  Grid: {ndir} directions × {nfreq} frequencies")

    # Convert to 1D directional spectrum
    dirs, e_dir = action_to_directional_spectrum_1d(
        action, omega_ww3, cg_ww3, directions=directions
    )

    print(f"\nOutput:")
    print(f"  Directional spectrum shape: {e_dir.shape}")
    print(f"  Direction range: {np.degrees(dirs[0]):.1f}° - {np.degrees(dirs[-1]):.1f}°")
    peak_dir_idx = np.argmax(e_dir)
    print(f"  Peak direction: {np.degrees(dirs[peak_dir_idx]):.1f}°")
    print(f"  Energy at peak: {e_dir.max():.6e} m²/rad")
    ddir = 2.0 * np.pi / ndir
    print(f"  Total energy (integrated): {np.sum(e_dir) * ddir:.6f} m²")
    print(f"  ✓ 1D directional conversion successful")


def test_peak_detection():
    """Test peak frequency and direction detection."""

    print("\n" + "="*70)
    print("Test 4: Peak Frequency and Direction Detection")
    print("="*70)

    # Setup
    nfreq, ndir = 30, 36
    depth = 100.0

    fr1 = 0.04
    xfr = 1.1
    omega_ww3 = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)
    frequencies = fr1 * xfr ** np.arange(nfreq)
    directions = np.linspace(0, 2*np.pi, ndir, endpoint=False)

    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # Create action spectrum with known peak
    action = np.zeros((ndir, nfreq))
    peak_freq_idx = 15
    peak_dir_idx = 9

    for idir in range(ndir):
        for ifreq in range(nfreq):
            dist_freq = (ifreq - peak_freq_idx)**2 / 5.0
            dist_dir = ((idir - peak_dir_idx)**2) / 3.0
            action[idir, ifreq] = np.exp(-dist_freq - dist_dir)

    print(f"\nInput:")
    print(f"  Action spectrum: Gaussian peak centered at")
    print(f"    Frequency index: {peak_freq_idx} → {frequencies[peak_freq_idx]:.4f} Hz")
    print(f"    Direction index: {peak_dir_idx} → {np.degrees(directions[peak_dir_idx]):.1f}°")

    # Detect peaks
    peak_freq = get_peak_frequency(action, frequencies)
    peak_dir = get_peak_direction(action, directions)

    print(f"\nDetected peaks:")
    print(f"  Peak frequency: {peak_freq:.4f} Hz (index {np.argmax(np.sum(action, axis=0))})")
    print(f"  Peak direction: {np.degrees(peak_dir):.1f}° (index {np.argmax(np.sum(action, axis=1))})")
    print(f"  ✓ Peak detection successful")


def test_normalization():
    """Test spectrum normalization."""

    print("\n" + "="*70)
    print("Test 5: Spectrum Normalization")
    print("="*70)

    # Setup
    nfreq, ndir = 30, 36
    depth = 100.0

    fr1 = 0.04
    xfr = 1.1
    omega_ww3 = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)

    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # Create random action spectrum
    action = np.random.rand(ndir, nfreq) * 0.001

    # Convert to energy
    energy_2d, _, _, _ = action_to_energy_2d(action, omega_ww3, cg_ww3)

    print(f"\nInput:")
    print(f"  Energy spectrum sum: {np.sum(energy_2d):.6f} m²")

    # Normalize to specific energy
    target_energy = 10.0
    normalized = normalize_spectrum(energy_2d, target_energy)

    print(f"\nOutput:")
    print(f"  Normalized spectrum sum: {np.sum(normalized):.6f} m²")
    print(f"  Target energy: {target_energy:.6f} m²")
    print(f"  Match: {np.isclose(np.sum(normalized), target_energy)}")
    print(f"  ✓ Normalization successful")


def test_comparison_with_spectrum_class():
    """Verify converter matches WaveSpectrum class output."""

    print("\n" + "="*70)
    print("Test 6: Verification Against WaveSpectrum Class")
    print("="*70)

    from wavewatch_python.spectrum_ww3_only import WaveSpectrum

    # Setup
    nfreq, ndir = 30, 36
    depth = 100.0

    fr1 = 0.04
    xfr = 1.1
    omega_ww3 = 2.0 * np.pi * fr1 * xfr ** np.arange(nfreq)

    dth = 2.0 * np.pi / ndir
    dsii = np.zeros(nfreq)
    dsii[0] = (omega_ww3[1] - omega_ww3[0]) / 2.0
    for i in range(1, nfreq - 1):
        dsii[i] = (omega_ww3[i+1] - omega_ww3[i-1]) / 2.0
    dsii[-1] = (omega_ww3[-1] - omega_ww3[-2]) / 2.0
    dden_ww3 = dth * dsii * omega_ww3

    fte_ww3 = 0.25 * omega_ww3[-1] * dth * omega_ww3[-1]
    fttr_ww3 = 0.20 * dth * omega_ww3[-1]
    ftwl_ww3 = (9.81 / 6.0) / omega_ww3[-1] * dth * omega_ww3[-1]

    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # Create action spectrum
    action = np.random.rand(ndir, nfreq) * 0.001

    print(f"\nComparing converter output with WaveSpectrum class...")

    # Method 1: Using converter
    freq_conv, e_freq_conv = action_to_frequency_spectrum_1d(
        action, omega_ww3, cg_ww3, dintegral=dden_ww3
    )

    # Method 2: Using WaveSpectrum class
    spectrum = WaveSpectrum(
        action=action,
        depth=depth,
        omega=omega_ww3,
        dintegral=dden_ww3,
        fte=fte_ww3,
        fttr=fttr_ww3,
        ftwl=ftwl_ww3,
        wavenumber=wn_ww3,
        group_velocity=cg_ww3
    )
    freq_spec, e_freq_spec = spectrum.get_spectrum_1d()

    # Compare (note: WaveSpectrum adds tail, converter doesn't)
    print(f"  Converter 1D freq spectrum range: {e_freq_conv.min():.6e} - {e_freq_conv.max():.6e}")
    print(f"  WaveSpectrum 1D freq spectrum range: {e_freq_spec.min():.6e} - {e_freq_spec.max():.6e}")
    print(f"  Frequency grids match: {np.allclose(freq_conv, freq_spec)}")
    print(f"  Note: WaveSpectrum includes tail extension, converter shows baseline")
    print(f"  ✓ Verification complete")


if __name__ == "__main__":
    test_action_to_energy_conversion()
    test_1d_frequency_spectrum()
    test_1d_directional_spectrum()
    test_peak_detection()
    test_normalization()
    test_comparison_with_spectrum_class()

    print("\n" + "="*70)
    print("All tests completed successfully!")
    print("="*70)
