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

    Attributes:
        action (ndarray): Action density spectrum A(f,θ) [m²·s·rad⁻¹]
        frequencies (ndarray): Frequency grid [Hz]
        directions (ndarray): Direction grid [radians]
        depth (float or ndarray): Water depth [m]
        freq_dim (int): Number of frequency bins
        dir_dim (int): Number of directional bins
    """

    def __init__(self, action, frequencies=None, directions=None, depth=1000.0,
                 freq_ratio=None, fr1=None):
        """
        Initialize wave spectrum object.

        Arguments:
            action (ndarray): Action density (ndir, nfreq) or (nfreq, ndir)
            frequencies (ndarray): Frequency grid [Hz], generated if None
            directions (ndarray): Direction grid [rad], generated if None
            depth (float or ndarray): Water depth [m]
            freq_ratio (float): Frequency ratio for grid generation
            fr1 (float): First frequency [Hz] for grid generation
        """

        self.action = np.atleast_2d(np.asarray(action, dtype=float))

        # Ensure action is (ndir, nfreq)
        if self.action.shape[0] < self.action.shape[1]:
            self.action = self.action.T

        self.dir_dim, self.freq_dim = self.action.shape

        # Generate frequency grid if not provided
        if frequencies is None:
            fr1 = fr1 or DEFAULT_FR1
            freq_ratio = freq_ratio or DEFAULT_XFR
            self.frequencies = self._generate_freq_grid(
                self.freq_dim, fr1, freq_ratio
            )
        else:
            self.frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))

        # Generate direction grid if not provided
        if directions is None:
            self.directions = self._generate_dir_grid(self.dir_dim)
        else:
            self.directions = np.atleast_1d(np.asarray(directions, dtype=float))

        self.depth = np.atleast_1d(np.asarray(depth, dtype=float))[0]

        # Compute angular frequencies
        self.omega = 2.0 * np.pi * self.frequencies

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
        # Integrate over frequencies
        e_dir = np.sum(self.action, axis=1) * self.dfreq

        return self.directions, e_dir
