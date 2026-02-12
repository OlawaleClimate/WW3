"""
Wave spectrum processing for WW3 pre-computed grid parameters.

Streamlined, production-ready module for computing spectral moments from action
density spectrum using pre-computed WW3 grid parameters (SIG, DDEN, WN, CG,
tail factors). Updated in v2.0 to align with spectrum_converter.py physics
corrections, properly handling coordinate transformations from A(k,θ) to E(f,θ).

Reference:
  - w3iogomd.F90 (W3OUTG subroutine)
  - w3gridmd.F90 (grid initialization)
  - spectrum_converter.py (coordinate transformation formulas)
"""

import numpy as np
from .constants import GRAV, TPI, SMALL
from .dispersion import solve_dispersion


class WaveSpectrum:
    """
    Wave spectrum object for WW3 pre-computed grid parameters.

    Requires all grid parameters to be pre-computed from WW3 initialization:
    - omega (SIG): Angular frequency array from w3gridmd.F90
    - dintegral (DDEN): Integration factors = DTH × DSII × SIG
    - fte, fttr, ftwl: Tail factors from w3gridmd.F90:3444-3448
    - wavenumber (WN): Wavenumber array from WAVNU1 subroutine
    - group_velocity (CG): Group velocity array from WAVNU1 subroutine

    Attributes:
        action (ndarray): Action density spectrum A(θ,f) [m²·s·rad⁻¹]
                         Shape: (ndir, nfreq)
        omega (ndarray): Angular frequency [rad/s]
        frequencies (ndarray): Frequency grid [Hz]
        directions (ndarray): Direction grid [radians]
        depth (float): Water depth [m]
        freq_dim (int): Number of frequency bins
        dir_dim (int): Number of directional bins
        dintegral (ndarray): Integration factors DDEN = DTH × DSII × SIG
        dfreq (ndarray): Frequency bandwidth for each bin
        ddir (float): Directional spacing [radians]
        fte, fttr, ftwl (float): Tail factors for f⁻⁵ extension
        wavenumber (ndarray): Wavenumber at each frequency [1/m]
        group_velocity (ndarray): Group velocity at each frequency [m/s]
    """

    def __init__(self, action, depth, omega, dintegral, fte, fttr, ftwl,
                 wavenumber, group_velocity, directions=None):
        """
        Initialize wave spectrum object with WW3 pre-computed parameters.

        Arguments:
            action (ndarray): Action density spectrum (ndir, nfreq) [m²·s·rad⁻¹]
            depth (float): Water depth [m]
            omega (ndarray): Angular frequencies from w3gridmd SIG [rad/s]
            dintegral (ndarray): Integration factors DDEN = DTH × DSII × SIG
            fte (float): Energy tail factor
            fttr (float): Period tail factor
            ftwl (float): Wavelength tail factor
            wavenumber (ndarray): Wavenumber from WAVNU1 [1/m]
            group_velocity (ndarray): Group velocity from WAVNU1 [m/s]
            directions (ndarray, optional): Direction grid [radians]
                                           Default: uniform spacing

        Raises:
            ValueError: If array dimensions don't match
            TypeError: If required parameters are missing
        """

        # Validate and reshape action spectrum
        self.action = np.atleast_2d(np.asarray(action, dtype=float))

        # Ensure action is (ndir, nfreq)
        if self.action.shape[0] < self.action.shape[1]:
            self.action = self.action.T

        self.dir_dim, self.freq_dim = self.action.shape
        self.depth = float(depth)

        # Validate and store WW3 pre-computed parameters
        self.omega = np.atleast_1d(np.asarray(omega, dtype=float))
        self.dintegral = np.atleast_1d(np.asarray(dintegral, dtype=float))
        self.wavenumber = np.atleast_1d(np.asarray(wavenumber, dtype=float))
        self.group_velocity = np.atleast_1d(np.asarray(group_velocity, dtype=float))

        # Validate dimension consistency
        if len(self.omega) != self.freq_dim:
            raise ValueError(
                f"omega length ({len(self.omega)}) != freq_dim ({self.freq_dim})"
            )
        if len(self.dintegral) != self.freq_dim:
            raise ValueError(
                f"dintegral length ({len(self.dintegral)}) != freq_dim ({self.freq_dim})"
            )
        if len(self.wavenumber) != self.freq_dim:
            raise ValueError(
                f"wavenumber length ({len(self.wavenumber)}) != freq_dim ({self.freq_dim})"
            )
        if len(self.group_velocity) != self.freq_dim:
            raise ValueError(
                f"group_velocity length ({len(self.group_velocity)}) != freq_dim ({self.freq_dim})"
            )

        # Convert angular frequency to Hz
        self.frequencies = self.omega / (2.0 * np.pi)

        # Generate or use provided directional grid
        if directions is None:
            self.directions = np.linspace(0, 2*np.pi, self.dir_dim, endpoint=False)
        else:
            self.directions = np.atleast_1d(np.asarray(directions, dtype=float))

        # Directional spacing
        self.ddir = 2.0 * np.pi / self.dir_dim

        # Compute frequency bandwidth from omega
        self.dfreq = self._compute_freq_bandwidth()

        # Store tail factors
        self.fte = float(fte)
        self.fttr = float(fttr)
        self.ftwl = float(ftwl)

    def _compute_freq_bandwidth(self):
        """
        Compute frequency bandwidth for each bin from omega.

        Uses central differences between adjacent frequencies.
        """
        dfreq = np.zeros(self.freq_dim)

        if self.freq_dim == 1:
            dfreq[0] = 1.0  # Fallback for single frequency
        else:
            # First bin: half-distance to next
            dfreq[0] = (self.omega[1] - self.omega[0]) / 2.0

            # Interior bins: central differences
            for i in range(1, self.freq_dim - 1):
                dfreq[i] = (self.omega[i+1] - self.omega[i-1]) / 2.0

            # Last bin: half-distance from previous
            dfreq[-1] = (self.omega[-1] - self.omega[-2]) / 2.0

        return dfreq

    def integrate_moments(self):
        """
        Integrate action density to compute spectral moments.

        Follows WW3 algorithm from w3iogomd.F90:1504-1575
        Converts action density to energy density using group velocity.

        Returns:
            dict: Computed moments and directional parameters
                - m0, m1, m2, m_1: Spectral moments
                - mom_x, mom_y: Directional moments (cosine, sine components)
                - mom_wn: Wavelength moment (for mean wavelength)
        """

        # Initialize moment arrays
        m0 = 0.0   # Total energy
        m1 = 0.0   # First moment (frequency-weighted)
        m2 = 0.0   # Second moment (frequency²-weighted)
        m_1 = 0.0  # Inverse moment (used for energy period)

        mom_x = 0.0  # E-W directional moment (cosine)
        mom_y = 0.0  # N-S directional moment (sine)
        mom_wn = 0.0  # Wavelength moment

        # Loop over frequencies
        for ifreq in range(self.freq_dim):
            # Integrate action over all directions
            action_band = np.sum(self.action[:, ifreq])

            # Energy conversion factor: DDEN / CG
            # Converts action density to energy density
            factor = self.dintegral[ifreq] / (self.group_velocity[ifreq] + SMALL)

            # Energy in frequency band
            energy_band = action_band * factor

            # Accumulate spectral moments
            m0 += energy_band
            m1 += energy_band * self.omega[ifreq]
            m2 += energy_band * self.omega[ifreq] ** 2
            m_1 += energy_band / (self.omega[ifreq] + SMALL)

            # Directional components
            cos_term = np.sum(self.action[:, ifreq] * np.cos(self.directions)) * factor
            sin_term = np.sum(self.action[:, ifreq] * np.sin(self.directions)) * factor
            mom_x += cos_term
            mom_y += sin_term

            # Wavelength moment
            mom_wn += energy_band / (self.wavenumber[ifreq] + SMALL)

        # Add tail extension (Pierson-Moskowitz f⁻⁵ tail)
        # Extends spectrum beyond cutoff frequency
        eband_tail = np.sum(self.action[:, -1]) / (
            self.group_velocity[-1] + SMALL
        )

        m0 += self.fte * eband_tail
        m1 += self.fte * 0.20 * eband_tail  # Approximate tail first moment
        m2 += self.fte * 0.5 * self.omega[-1] ** 4 * self.ddir * eband_tail
        m_1 += self.fttr * eband_tail
        mom_wn += self.ftwl * eband_tail

        # Directional tail components
        cos_tail = self.fte * np.sum(
            self.action[:, -1] * np.cos(self.directions)
        ) / (self.group_velocity[-1] + SMALL)
        sin_tail = self.fte * np.sum(
            self.action[:, -1] * np.sin(self.directions)
        ) / (self.group_velocity[-1] + SMALL)
        mom_x += cos_tail
        mom_y += sin_tail

        return {
            'm0': m0,
            'm1': m1,
            'm2': m2,
            'm_1': m_1,
            'mom_x': mom_x,
            'mom_y': mom_y,
            'mom_wn': mom_wn,
        }

    def compute_parameters(self):
        """
        Compute all wave parameters from spectral moments.

        Uses spectral moment method as implemented in WW3:
        - HS: Significant wave height = 4√m₀
        - T01: Mean period = 2π(m₀/m₁)
        - T02: Zero-crossing period = 2π√(m₀/m₂)
        - T0M1: Energy period = 2π(m₋₁/m₀)
        - TP: Peak period = 2π/ω_peak
        - THM: Mean direction = atan2(mom_y, mom_x)
        - THS: Directional spread = √(2(1 - √((mom_x² + mom_y²)/m₀²)))
        - WLM: Mean wavelength = 2π(mom_wn/m₀)

        Returns:
            dict: Wave parameters
                - hs: Significant wave height [m]
                - t01, t02, t0m1: Wave periods [s]
                - tp: Peak period [s]
                - thm: Mean direction [rad]
                - ths: Directional spread [rad]
                - wlm: Mean wavelength [m]
                - width: Spectral width
                - moments: Dictionary of raw spectral moments
                - depth: Water depth [m]
        """

        # Compute moments
        moments = self.integrate_moments()

        m0 = moments['m0']
        m1 = moments['m1']
        m2 = moments['m2']
        m_1 = moments['m_1']
        mom_x = moments['mom_x']
        mom_y = moments['mom_y']
        mom_wn = moments['mom_wn']

        results = {}

        # ===== WAVE HEIGHTS =====
        # Significant Wave Height: HS = 4√m₀
        if m0 > SMALL:
            results['hs'] = 4.0 * np.sqrt(m0)
        else:
            results['hs'] = 0.0

        # ===== WAVE PERIODS =====
        if m0 > SMALL and m1 > SMALL and m2 > SMALL:
            # Mean Period: T01 = 2π(m₀/m₁)
            results['t01'] = TPI * m0 / m1

            # Zero-crossing Period: T02 = 2π√(m₀/m₂)
            results['t02'] = TPI * np.sqrt(m0 / m2)

            # Energy Period: T0M1 = 2π(m₋₁/m₀)
            if m_1 > SMALL:
                results['t0m1'] = TPI * m_1 / m0
            else:
                # Fallback: use last frequency
                results['t0m1'] = TPI / self.omega[-1]

            # Peak Period: Tp = 2π/ω_peak
            peak_freq_idx = np.argmax(self.action.sum(axis=0))
            results['tp'] = TPI / self.omega[peak_freq_idx]

            # Mean Wavelength: λ = 2π(mom_wn/m₀)
            results['wlm'] = TPI * mom_wn / m0

        else:
            results['t01'] = 0.0
            results['t02'] = 0.0
            results['t0m1'] = 0.0
            results['tp'] = 0.0
            results['wlm'] = 0.0

        # ===== DIRECTIONAL PROPERTIES =====
        # Mean Direction: θ = atan2(mom_y, mom_x)
        if (np.abs(mom_x) + np.abs(mom_y)) > SMALL:
            results['thm'] = np.arctan2(mom_y, mom_x)
        else:
            results['thm'] = 0.0

        # Directional Spread: σθ = √(2(1 - √((mom_x² + mom_y²)/m₀²)))
        if m0 > SMALL:
            dir_spread_arg = (mom_x**2 + mom_y**2) / (m0**2 + SMALL)
            dir_spread_arg = np.clip(dir_spread_arg, 0.0, 1.0)
            results['ths'] = np.sqrt(2.0 * (1.0 - np.sqrt(dir_spread_arg)))
        else:
            results['ths'] = 0.0

        # ===== SPECTRAL WIDTH =====
        # Spectral Width: σf = √(m₂m₀/m₁² - 1)
        if m0 > SMALL and m1 > SMALL:
            width_arg = m2 * m0 / (m1**2 + SMALL) - 1.0
            width_arg = max(0.0, width_arg)
            results['width'] = np.sqrt(width_arg)
        else:
            results['width'] = 0.0

        # ===== STORE INTERNAL RESULTS =====
        results['moments'] = moments
        results['depth'] = self.depth

        return results

    def get_spectrum_1d(self):
        """
        Extract 1D frequency spectrum by integrating over directions.

        Follows w3iogomd.F90 frequency spectrum computation:
        E(f) = ∑_θ [A(θ,f) × (DDEN/CG)] × DTH

        Converts action density to energy density and integrates over directions.

        Returns:
            tuple:
                - freq (ndarray): Frequency array [Hz]
                - e_freq (ndarray): 1D energy spectrum [m²/Hz]
        """
        # Direction-integrated action per frequency bin
        action_band = np.sum(self.action, axis=0)

        # Energy conversion factor: DDEN / CG
        factor = self.dintegral / (self.group_velocity + SMALL)

        # Energy per frequency band
        ebd = action_band * factor

        # Integrate over directions with directional bin width
        # e_freq = EBD × DTH
        e_freq = ebd * self.ddir

        return self.frequencies, e_freq

    def get_spectrum_directional_1d(self):
        """
        Extract 1D directional spectrum by integrating over frequencies.

        Follows w3iogomd.F90 directional spectrum computation:
        E(θ) = ∑_f [A(θ,f) × (DDEN/CG)]

        Converts action density to energy density and integrates over frequencies.

        Returns:
            tuple:
                - dirs (ndarray): Direction array [radians]
                - e_dir (ndarray): 1D directional spectrum [m²/rad]
        """
        # Energy conversion factor: DDEN / CG
        factor = self.dintegral / (self.group_velocity + SMALL)

        # Integrate over frequencies with frequency bandwidth for each
        # e_dir = ∑_f [A(θ,f) × factor × DSII]
        # Where DSII is computed from dfreq (frequency bandwidth)
        e_dir = np.sum(self.action * factor[np.newaxis, :] * self.dfreq[np.newaxis, :], axis=1)

        return self.directions, e_dir
