"""
Example usage of WaveWatch Python package.

Demonstrates how to:
  1. Create a wave spectrum from action density
  2. Compute spectral moments
  3. Calculate wave parameters
  4. Plot results
"""

import numpy as np
import matplotlib.pyplot as plt
from spectrum import WaveSpectrum


def create_jonswap_spectrum(nfreq=30, ndir=36, tp=8.0, hs=2.0, gamma=3.3):
    """
    Create a synthetic JONSWAP spectrum (2D action density).

    Arguments:
        nfreq (int): Number of frequency bins
        ndir (int): Number of directional bins
        tp (float): Peak period [s]
        hs (float): Significant wave height [m]
        gamma (float): JONSWAP peak parameter

    Returns:
        action (ndarray): Action density spectrum (ndir, nfreq)
    """

    # Generate frequency grid
    fr1 = 0.04
    xfr = 1.1
    freq = fr1 * xfr ** np.arange(nfreq)
    omega = 2.0 * np.pi * freq

    # Generate directional grid
    directions = np.linspace(0, 2*np.pi, ndir, endpoint=False)

    # JONSWAP spectrum: E(f) = alpha*g²/(2π)^4 * f^(-5) * exp(-5/4*(f/fp)^(-4)) * gamma^exp
    fp = 1.0 / tp
    alpha = hs**2 / (0.0081 * 9.81**2 / (2*np.pi)**4 * (1.0/fp)**4)

    # Frequency spectrum
    e_freq = np.zeros(nfreq)
    for i, f in enumerate(freq):
        f_ratio = f / fp
        sigma = 0.07 if f <= fp else 0.09
        exp_term = np.exp(-(f_ratio - 1.0)**2 / (2 * sigma**2))
        gamma_term = gamma**exp_term
        e_freq[i] = (alpha * 9.81**2 / (2*np.pi)**4 * f**(-5) *
                     np.exp(-1.25 * (fp/f)**4) * gamma_term)

    # Directional spreading (cosine-2s)
    s = 2
    dir_spreading = (2**s * np.math.factorial(s)**2 / np.math.pi /
                     np.math.factorial(2*s) * np.cos(directions)**2s)
    dir_spreading = dir_spreading / np.sum(dir_spreading) * ndir

    # 2D spectrum: E(f,θ) = E(f) × D(θ)
    spectrum_2d = np.outer(dir_spreading, e_freq)

    # Convert energy to action density: N(f,θ) = E(f,θ) / (g/(2π)f) = E(f,θ)/(g*T)
    # For simplicity, just use energy normalized
    action = spectrum_2d / (9.81 * omega[np.newaxis, :] + 1e-10)

    return action


def example_1_basic_spectrum():
    """Example 1: Basic spectrum with random action density."""

    print("=" * 70)
    print("EXAMPLE 1: Basic Spectrum with Random Action Density")
    print("=" * 70)

    # Create random action density
    ndir, nfreq = 36, 30
    action = np.random.rand(ndir, nfreq) * 0.001

    # Initialize spectrum
    spectrum = WaveSpectrum(action, depth=100.0)

    # Compute parameters
    params = spectrum.compute_parameters()

    # Display results
    print(f"\nWater Depth: {params['depth']:.1f} m")
    print(f"\nWave Parameters:")
    print(f"  Significant Wave Height (HS): {params['hs']:.3f} m")
    print(f"  Mean Period (T01): {params['t01']:.2f} s")
    print(f"  Zero-crossing Period (T02): {params['t02']:.2f} s")
    print(f"  Energy Period (T0M1): {params['t0m1']:.2f} s")
    print(f"  Peak Period (Tp): {params['tp']:.2f} s")
    print(f"  Mean Direction: {np.degrees(params['thm']):.1f}°")
    print(f"  Directional Spread: {np.degrees(params['ths']):.1f}°")
    print(f"  Mean Wavelength: {params['wlm']:.1f} m")
    print(f"  Spectral Width: {params['width']:.3f}")

    print(f"\nSpectral Moments:")
    print(f"  M0 (Energy): {params['moments']['m0']:.6f} m²")
    print(f"  M1 (Freq-weighted): {params['moments']['m1']:.6f} m²·s")
    print(f"  M2 (Freq²-weighted): {params['moments']['m2']:.6f} m²·s²")
    print(f"  M-1 (Inverse): {params['moments']['m_1']:.6f} m²·s")


def example_2_jonswap_spectrum():
    """Example 2: JONSWAP spectrum (realistic wind sea)."""

    print("\n" + "=" * 70)
    print("EXAMPLE 2: JONSWAP Spectrum (Realistic Wind Sea)")
    print("=" * 70)

    # Create JONSWAP spectrum
    action = create_jonswap_spectrum(nfreq=30, ndir=36, tp=8.0, hs=2.5)

    # Test different water depths
    depths = [5.0, 20.0, 100.0, 1000.0]

    print(f"\nJONSWAP Spectrum: Tp=8.0s, Hs=2.5m, γ=3.3")
    print(f"\nWave Parameters at Different Depths:")
    print(f"{'Depth (m)':>10} | {'HS (m)':>8} | {'T01 (s)':>8} | {'T02 (s)':>8} | "
          f"{'T0M1 (s)':>8} | {'λ (m)':>8}")
    print("-" * 70)

    for depth in depths:
        spectrum = WaveSpectrum(action, depth=depth)
        params = spectrum.compute_parameters()

        print(f"{depth:>10.1f} | {params['hs']:>8.3f} | {params['t01']:>8.2f} | "
              f"{params['t02']:>8.2f} | {params['t0m1']:>8.2f} | {params['wlm']:>8.1f}")

    # Detailed output for deep water
    print(f"\nDetailed Output (Deep Water, h=1000m):")
    spectrum = WaveSpectrum(action, depth=1000.0)
    params = spectrum.compute_parameters()

    print(f"  HS: {params['hs']:.3f} m")
    print(f"  T01: {params['t01']:.2f} s")
    print(f"  T02: {params['t02']:.2f} s")
    print(f"  T0M1: {params['t0m1']:.2f} s")
    print(f"  Mean Direction: {np.degrees(params['thm']):.1f}°")
    print(f"  Directional Spread: {np.degrees(params['ths']):.1f}°")


def example_3_depth_effect():
    """Example 3: Demonstrate how water depth affects wave parameters."""

    print("\n" + "=" * 70)
    print("EXAMPLE 3: Water Depth Effect on Wave Parameters")
    print("=" * 70)

    # Create a fixed spectrum
    action = create_jonswap_spectrum(nfreq=30, ndir=36, tp=10.0, hs=3.0)

    # Compute parameters at different depths
    depths = np.logspace(0, 3.5, 20)  # 1 to 3162 meters
    results = {
        'depth': [],
        'hs': [],
        't01': [],
        't02': [],
        't0m1': [],
    }

    for depth in depths:
        spectrum = WaveSpectrum(action, depth=depth)
        params = spectrum.compute_parameters()

        results['depth'].append(depth)
        results['hs'].append(params['hs'])
        results['t01'].append(params['t01'])
        results['t02'].append(params['t02'])
        results['t0m1'].append(params['t0m1'])

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Water Depth Effect on Wave Parameters', fontsize=14, fontweight='bold')

    # HS vs Depth
    ax = axes[0, 0]
    ax.semilogx(results['depth'], results['hs'], 'o-', linewidth=2)
    ax.set_xlabel('Water Depth [m]')
    ax.set_ylabel('Significant Wave Height [m]')
    ax.set_title('HS vs Depth')
    ax.grid(True, alpha=0.3)

    # T01 vs Depth
    ax = axes[0, 1]
    ax.semilogx(results['depth'], results['t01'], 'o-', linewidth=2, color='orange')
    ax.set_xlabel('Water Depth [m]')
    ax.set_ylabel('Mean Period [s]')
    ax.set_title('T01 vs Depth')
    ax.grid(True, alpha=0.3)

    # T02 vs Depth
    ax = axes[1, 0]
    ax.semilogx(results['depth'], results['t02'], 'o-', linewidth=2, color='green')
    ax.set_xlabel('Water Depth [m]')
    ax.set_ylabel('Zero-crossing Period [s]')
    ax.set_title('T02 vs Depth')
    ax.grid(True, alpha=0.3)

    # T0M1 vs Depth
    ax = axes[1, 1]
    ax.semilogx(results['depth'], results['t0m1'], 'o-', linewidth=2, color='red')
    ax.set_xlabel('Water Depth [m]')
    ax.set_ylabel('Energy Period [s]')
    ax.set_title('T0M1 vs Depth')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/wavewatch_python/depth_effect.png', dpi=100, bbox_inches='tight')
    print("\nPlot saved: /home/user/wavewatch_python/depth_effect.png")
    plt.close()


def example_4_1d_spectra():
    """Example 4: Compute and plot 1D frequency and directional spectra."""

    print("\n" + "=" * 70)
    print("EXAMPLE 4: 1D Frequency and Directional Spectra")
    print("=" * 70)

    # Create JONSWAP spectrum
    action = create_jonswap_spectrum(nfreq=30, ndir=36, tp=8.0, hs=2.5)

    # Initialize spectrum
    spectrum = WaveSpectrum(action, depth=100.0)

    # Get 1D spectra
    freq, e_freq = spectrum.get_spectrum_1d()
    dirs, e_dir = spectrum.get_spectrum_directional_1d()

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('1D Wave Spectra (h=100m)', fontsize=14, fontweight='bold')

    # Frequency spectrum
    ax = axes[0]
    ax.plot(freq, e_freq, 'b-', linewidth=2)
    ax.fill_between(freq, e_freq, alpha=0.3)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Energy Density [m²/Hz]')
    ax.set_title('Frequency Spectrum E(f)')
    ax.grid(True, alpha=0.3)

    # Directional spectrum
    ax = axes[1]
    dirs_deg = np.degrees(dirs)
    ax.plot(dirs_deg, e_dir, 'r-', linewidth=2)
    ax.fill_between(dirs_deg, e_dir, alpha=0.3)
    ax.set_xlabel('Direction [°]')
    ax.set_ylabel('Energy Density [m²/rad]')
    ax.set_title('Directional Spectrum D(θ)')
    ax.set_xlim(0, 360)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/wavewatch_python/1d_spectra.png', dpi=100, bbox_inches='tight')
    print("\nPlot saved: /home/user/wavewatch_python/1d_spectra.png")
    plt.close()


def example_0_ww3_parameters():
    """Example 0: Using WW3 pre-computed parameters (RECOMMENDED)."""

    print("=" * 70)
    print("EXAMPLE 0: Using WW3 Pre-computed Parameters (RECOMMENDED)")
    print("=" * 70)

    # Simulate WW3 grid initialization values
    # These would be read from a WW3 model run
    nfreq, ndir = 30, 36

    # Frequency grid from WW3 (SIG array in w3gridmd.F90)
    fr1_ww3 = 0.04
    xfr_ww3 = 1.1
    omega_ww3 = 2.0 * np.pi * fr1_ww3 * xfr_ww3 ** np.arange(nfreq)

    # Integration factors from WW3 (DDEN in w3gridmd.F90)
    dth = 2.0 * np.pi / ndir
    dsii = np.zeros(nfreq)
    dsii[0] = (omega_ww3[1] - omega_ww3[0]) / 2.0
    for i in range(1, nfreq - 1):
        dsii[i] = (omega_ww3[i+1] - omega_ww3[i-1]) / 2.0
    dsii[-1] = (omega_ww3[-1] - omega_ww3[-2]) / 2.0

    dden_ww3 = dth * dsii * omega_ww3

    # Tail factors from WW3 (computed in w3gridmd.F90)
    fte_ww3 = 0.25 * omega_ww3[-1] * dth * omega_ww3[-1]
    fttr_ww3 = 0.20 * dth * omega_ww3[-1]
    ftwl_ww3 = (9.81 / 6.0) / omega_ww3[-1] * dth * omega_ww3[-1]

    # Wavenumber and group velocity from WW3 (computed in w3initmd.F90)
    from dispersion import solve_dispersion
    depth = 100.0
    wn_ww3 = np.zeros(nfreq)
    cg_ww3 = np.zeros(nfreq)
    for ik in range(nfreq):
        wn_ww3[ik], cg_ww3[ik] = solve_dispersion(omega_ww3[ik], depth)

    # Action density spectrum (from WW3 simulation output)
    action = create_jonswap_spectrum(nfreq=nfreq, ndir=ndir, tp=8.0, hs=2.5)

    print(f"\nWW3 Grid Parameters (from w3gridmd.F90):")
    print(f"  NK={nfreq}, NTH={ndir}")
    print(f"  omega[0] = {omega_ww3[0]:.4f} rad/s (f={omega_ww3[0]/(2*np.pi):.4f} Hz)")
    print(f"  omega[-1] = {omega_ww3[-1]:.4f} rad/s (f={omega_ww3[-1]/(2*np.pi):.4f} Hz)")
    print(f"  FTE (tail energy) = {fte_ww3:.6f}")
    print(f"  FTTR (tail period) = {fttr_ww3:.6f}")
    print(f"  FTWL (tail wavelength) = {ftwl_ww3:.6f}")
    print(f"  Water depth = {depth:.1f} m")

    # Initialize spectrum with WW3 parameters (RECOMMENDED METHOD)
    spectrum = WaveSpectrum(
        action,
        depth=depth,
        # Pre-computed WW3 grid parameters
        omega=omega_ww3,
        dintegral=dden_ww3,
        fte=fte_ww3,
        fttr=fttr_ww3,
        ftwl=ftwl_ww3,
        wavenumber=wn_ww3,
        group_velocity=cg_ww3
    )

    # Compute parameters
    params = spectrum.compute_parameters()

    print(f"\nComputed Wave Parameters:")
    print(f"  HS: {params['hs']:.3f} m")
    print(f"  T01: {params['t01']:.2f} s")
    print(f"  T02: {params['t02']:.2f} s")
    print(f"  T0M1: {params['t0m1']:.2f} s")
    print(f"  Mean Direction: {np.degrees(params['thm']):.1f}°")
    print(f"  Spectral Moments:")
    print(f"    M0: {params['moments']['m0']:.6f} m²")
    print(f"    M1: {params['moments']['m1']:.6f} m²·s")


if __name__ == "__main__":
    # Run examples
    print("\n\n")
    example_0_ww3_parameters()
    example_1_basic_spectrum()
    example_2_jonswap_spectrum()
    example_3_depth_effect()
    example_4_1d_spectra()

    print("\n" + "=" * 70)
    print("All examples completed!")
    print("=" * 70)
