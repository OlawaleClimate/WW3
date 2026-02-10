"""
Wave spectrum processing and moment calculations.

Core module for computing spectral moments from action density spectrum.
Reference: w3iogomd.F90 (W3OUTG subroutine)
"""

import numpy as np
from .constants import (
    GRAV, TPI, TPIINV, SMALL, UNDEF,
    DEFAULT_NK, DEFAULT_NTH, DEFAULT_FR1, DEFAULT_XFR, DEFAULT_DMIN
)
from .dispersion import solve_dispersion


class WaveSpectrum:
    """
    Wave spectrum object for computing wave parameters from action density.

    Can accept either:
    1. Pre-computed WW3 grid parameters (SIG, DSII/DDEN, FTE, etc.)
    2. Raw frequency/direction grids (will auto-compute factors)

    Attributes:
        action (ndarray): Action density spectrum A(f,θ) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s]
        frequencies (ndarray): Frequency grid [Hz]
        directions (ndarray): Direction grid [radians]
        depth (float): Water depth [m]
        freq_dim (int): Number of frequency bins
        dir_dim (int): Number of directional bins
        dintegral (ndarray): Integration factors DDEN = DTH × DSII × SIG
        fte, fttr, ftwl (float): Tail factors
    """

    def __init__(self, action, depth=1000.0,
                 # WW3 pre-computed parameters (preferred)
                 omega=None, dintegral=None, fte=None, fttr=None, ftwl=None,
                 wavenumber=None, group_velocity=None,
                 # Alternative: raw grid parameters
                 frequencies=None, directions=None,
                 freq_ratio=None, fr1=None):
        """
        Initialize wave spectrum object.

        Arguments:
            action (ndarray): Action density (ndir, nfreq) or (nfreq, ndir)
            depth (float or ndarray): Water depth [m]

            ===== WW3 PRE-COMPUTED PARAMETERS (Preferred) =====
            omega (ndarray): Angular frequencies [rad/s] from SIG in w3gridmd
            dintegral (ndarray): Integration factors DDEN = DTH × DSII × SIG
            fte (float): Energy tail factor (from w3gridmd)
            fttr (float): Period tail factor (from w3gridmd)
            ftwl (float): Wavelength tail factor (from w3gridmd)
            wavenumber (ndarray): Pre-computed WN(IK) from WAVNU1 [1/m]
            group_velocity (ndarray): Pre-computed CG(IK) from WAVNU1 [m/s]

            ===== ALTERNATIVE: Raw grid (auto-computed) =====
            frequencies (ndarray): Frequency grid [Hz]
            directions (ndarray): Direction grid [rad]
            freq_ratio (float): Frequency ratio for generation (default XFR=1.1)
            fr1 (float): First frequency [Hz] for generation

        Examples:
            # Method 1: Using WW3 pre-computed parameters (RECOMMENDED)
            >>> omega_ww3 = np.array([...])  # from w3gridmd SIG array
            >>> dintegral_ww3 = np.array([...])  # from w3gridmd DDEN
            >>> fte_ww3 = 0.25 * omega_ww3[-1] * dth * omega_ww3[-1]
            >>> spectrum = WaveSpectrum(action_from_ww3, depth=100.0,
            ...                         omega=omega_ww3,
            ...                         dintegral=dintegral_ww3,
            ...                         fte=fte_ww3, fttr=fttr_ww3, ftwl=ftwl_ww3,
            ...                         wavenumber=wn_from_ww3,
            ...                         group_velocity=cg_from_ww3)

            # Method 2: Using raw frequency grid (auto-computed)
            >>> spectrum = WaveSpectrum(action, depth=100.0,
            ...                         frequencies=freq_array,
            ...                         directions=dir_array)
        """

        self.action = np.atleast_2d(np.asarray(action, dtype=float))

        # Ensure action is (ndir, nfreq)
        if self.action.shape[0] < self.action.shape[1]:
            self.action = self.action.T

        self.dir_dim, self.freq_dim = self.action.shape
        self.depth = np.atleast_1d(np.asarray(depth, dtype=float))[0]

        # ===== USE WW3 PRE-COMPUTED PARAMETERS IF PROVIDED =====
        if omega is not None:
            # Using WW3 pre-computed grid parameters
            self.omega = np.atleast_1d(np.asarray(omega, dtype=float))
            self.frequencies = self.omega / (2.0 * np.pi)

            # Use provided integration factors
            if dintegral is not None:
                self.dintegral = np.atleast_1d(np.asarray(dintegral, dtype=float))
            else:
                # Compute from omega if dintegral not provided
                self._compute_integration_factors_from_omega()

            # Use provided tail factors
            self.fte = fte if fte is not None else self._default_fte()
            self.fttr = fttr if fttr is not None else self._default_fttr()
            self.ftwl = ftwl if ftwl is not None else self._default_ftwl()

            # Generate or use provided directional grid
            if directions is None:
                self.directions = self._generate_dir_grid(self.dir_dim)
            else:
                self.directions = np.atleast_1d(np.asarray(directions, dtype=float))

            self.ddir = 2.0 * np.pi / self.dir_dim

            # Use provided or compute wavenumber/group velocity
            if wavenumber is not None and group_velocity is not None:
                self.wavenumber = np.atleast_1d(np.asarray(wavenumber, dtype=float))
                self.group_velocity = np.atleast_1d(np.asarray(group_velocity, dtype=float))
            else:
                # Compute from dispersion relation
                self._compute_dispersion()

        # ===== USE RAW GRID PARAMETERS (Auto-compute factors) =====
        else:
            # Generate frequency grid if not provided
            if frequencies is None:
                fr1 = fr1 or DEFAULT_FR1
                freq_ratio = freq_ratio or DEFAULT_XFR
                self.frequencies = self._generate_freq_grid(
                    self.freq_dim, fr1, freq_ratio
                )
            else:
                self.frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))

            # Compute angular frequencies
            self.omega = 2.0 * np.pi * self.frequencies

            # Generate direction grid if not provided
            if directions is None:
                self.directions = self._generate_dir_grid(self.dir_dim)
            else:
                self.directions = np.atleast_1d(np.asarray(directions, dtype=float))

            self.ddir = 2.0 * np.pi / self.dir_dim

            # Compute integration factors
            self._compute_integration_factors()

            # Compute wavenumber and group velocity
            self._compute_dispersion()

    @staticmethod
    def _generate_freq_grid(nfreq, fr1, freq_ratio):
        """Generate logarithmic frequency grid."""
        return fr1 * freq_ratio ** np.arange(nfreq)

    @staticmethod
    def _generate_dir_grid(ndir):
        """Generate regular directional grid."""
        return np.linspace(0, 2*np.pi, ndir, endpoint=False)

    def _default_fte(self):
        """Default FTE from peak omega (Pierson-Moskowitz)."""
        omega_peak = self.omega[-1]
        return 0.25 * omega_peak * self.ddir * omega_peak

    def _default_fttr(self):
        """Default FTTR (tail period factor)."""
        omega_peak = self.omega[-1]
        return 0.20 * self.ddir * omega_peak

    def _default_ftwl(self):
        """Default FTWL (tail wavelength factor)."""
        omega_peak = self.omega[-1]
        return (GRAV / 6.0) / omega_peak * self.ddir * omega_peak

    def _compute_integration_factors_from_omega(self):
        """Compute DDEN from omega when only omega is provided."""
        # Compute frequency bandwidth from omega
        dfreq = np.zeros(len(self.omega))
        dfreq[0] = (self.omega[1] - self.omega[0]) / 2.0
        for i in range(1, len(self.omega) - 1):
            dfreq[i] = (self.omega[i+1] - self.omega[i-1]) / 2.0
        dfreq[-1] = (self.omega[-1] - self.omega[-2]) / 2.0

        # Full integration factor
        self.dintegral = self.ddir * dfreq * self.omega

    def _compute_integration_factors(self):
        """Compute frequency bandwidth and integration factors."""

        # Frequency bandwidth (DSII in WW3)
        self.dfreq = np.zeros(self.freq_dim)
        self.dfreq[0] = (self.frequencies[1] - self.frequencies[0]) / 2.0
        for i in range(1, self.freq_dim - 1):
            self.dfreq[i] = (self.frequencies[i+1] - self.frequencies[i-1]) / 2.0
        self.dfreq[-1] = (self.frequencies[-1] - self.frequencies[-2]) / 2.0

        # Directional step
        self.ddir = 2.0 * np.pi / self.dir_dim

        # Full integration factor: DDEN = ddir × dfreq × omega
        self.dintegral = self.ddir * self.dfreq * self.omega

        # Tail factors (Pierson-Moskowitz f⁻⁵ tail)
        omega_peak = self.omega[-1]
        self.fte = 0.25 * omega_peak * self.ddir * omega_peak
        self.ftwl = (GRAV / 6.0) / omega_peak * self.ddir * omega_peak
        self.fttr = 0.20 * self.ddir * omega_peak

    def _compute_dispersion(self):
        """Compute wavenumber and group velocity."""

        self.wavenumber = np.zeros(self.freq_dim)
        self.group_velocity = np.zeros(self.freq_dim)

        for ifreq in range(self.freq_dim):
            wn, cg = solve_dispersion(self.omega[ifreq], self.depth)
            self.wavenumber[ifreq] = wn
            self.group_velocity[ifreq] = cg

    def integrate_moments(self):
        """
        Integrate action density to compute spectral moments.

        Returns a dictionary with computed moments and parameters.
        """

        # Initialize moment arrays
        m0 = 0.0   # Total energy
        m1 = 0.0   # First moment
        m2 = 0.0   # Second moment
        m_1 = 0.0  # Inverse moment

        mom_x = 0.0  # E-W directional moment
        mom_y = 0.0  # N-S directional moment
        mom_wn = 0.0  # Wavelength moment

        # Loop over frequencies
        for ifreq in range(self.freq_dim):
            # Integrate action over directions
            action_band = np.sum(self.action[:, ifreq])

            # Energy conversion factor
            factor = self.dintegral[ifreq] / (self.group_velocity[ifreq] + SMALL)

            # Energy in frequency band
            energy_band = action_band * factor

            # Accumulate moments
            m0 += energy_band
            m1 += energy_band * self.omega[ifreq]
            m2 += energy_band * self.omega[ifreq]**2
            m_1 += energy_band / (self.omega[ifreq] + SMALL)

            # Directional components
            mom_x += np.sum(self.action[:, ifreq] * np.cos(self.directions)) * factor
            mom_y += np.sum(self.action[:, ifreq] * np.sin(self.directions)) * factor

            # Wavelength moment
            mom_wn += energy_band / (self.wavenumber[ifreq] + SMALL)

        # Add tail extension (Pierson-Moskowitz f⁻⁵ tail)
        eband_tail = np.sum(self.action[:, -1]) / (self.group_velocity[-1] + SMALL)
        m0 += self.fte * eband_tail
        m1 += self.fte * 0.20 * eband_tail  # Approximate tail first moment
        m2 += self.fte * 0.5 * self.omega[-1]**4 * self.ddir * eband_tail
        m_1 += self.fttr * eband_tail
        mom_wn += self.ftwl * eband_tail
        mom_x += self.fte * np.sum(self.action[:, -1] * np.cos(self.directions)) / (
            self.group_velocity[-1] + SMALL
        )
        mom_y += self.fte * np.sum(self.action[:, -1] * np.sin(self.directions)) / (
            self.group_velocity[-1] + SMALL
        )

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

        Returns a dictionary with wave parameters:
            hs: Significant wave height [m]
            t01: Mean period [s]
            t02: Zero-crossing period [s]
            t0m1: Energy period [s]
            thm: Mean direction [rad]
            ths: Directional spread [rad]
            wlm: Mean wavelength [m]
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

        # Significant Wave Height: HS = 4√m₀
        if m0 > SMALL:
            results['hs'] = 4.0 * np.sqrt(m0)
        else:
            results['hs'] = 0.0

        # Periods (only if energy is significant)
        if m0 > SMALL and m1 > SMALL and m2 > SMALL:
            # Mean Period: T01 = 2π(m0/m1)
            results['t01'] = TPI * m0 / m1

            # Zero-crossing Period: T02 = 2π√(m0/m2)
            results['t02'] = TPI * np.sqrt(m0 / m2)

            # Energy Period: T0M1 = 2π(m_1/m0)
            if m_1 > SMALL:
                results['t0m1'] = TPI * m_1 / m0
            else:
                results['t0m1'] = TPI / self.omega[-1]

            # Mean Wavelength: λ = 2π(mom_wn/m0)
            results['wlm'] = TPI * mom_wn / m0

            # Peak Period: Tp = 2π/ω_peak
            results['tp'] = TPI / self.omega[np.argmax(self.action.sum(axis=0))]

        else:
            results['t01'] = 0.0
            results['t02'] = 0.0
            results['t0m1'] = 0.0
            results['wlm'] = 0.0
            results['tp'] = 0.0

        # Mean Direction: θ = atan2(mom_y, mom_x)
        if (np.abs(mom_x) + np.abs(mom_y)) > SMALL:
            results['thm'] = np.arctan2(mom_y, mom_x)
        else:
            results['thm'] = 0.0

        # Directional Spread: σθ = √(2(1 - √((mom_x² + mom_y²)/m0²)))
        if m0 > SMALL:
            dir_spread_arg = (mom_x**2 + mom_y**2) / (m0**2 + SMALL)
            dir_spread_arg = np.clip(dir_spread_arg, 0.0, 1.0)
            results['ths'] = np.sqrt(2.0 * (1.0 - np.sqrt(dir_spread_arg)))
        else:
            results['ths'] = 0.0

        # Spectral Width: σf = √(m2*m0/m1² - 1)
        if m0 > SMALL and m1 > SMALL:
            width_arg = m2 * m0 / (m1**2 + SMALL) - 1.0
            width_arg = max(0.0, width_arg)
            results['width'] = np.sqrt(width_arg)
        else:
            results['width'] = 0.0

        # Store moments for reference
        results['moments'] = moments
        results['depth'] = self.depth

        return results

    def get_spectrum_1d(self):
        """
        Get 1D frequency spectrum by integrating over directions.

        Returns:
            freq (ndarray): Frequency [Hz]
            e_freq (ndarray): Energy density [m²/Hz]
        """
        # Integrate over directions
        e_freq = np.sum(self.action, axis=0) * self.ddir

        return self.frequencies, e_freq

    def get_spectrum_directional_1d(self):
        """
        Get 1D directional spectrum by integrating over frequencies.

        Returns:
            dirs (ndarray): Direction [radians]
            e_dir (ndarray): Energy density [m²/rad]
        """
        # Integrate over frequencies (multiply by dfreq for each frequency)
        e_dir = np.sum(self.action * self.dfreq[np.newaxis, :], axis=1)

        return self.directions, e_dir
