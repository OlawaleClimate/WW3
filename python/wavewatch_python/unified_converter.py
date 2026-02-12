"""
Unified WaveWatch III Spectrum Converter

Single converter that intelligently handles both WW3 and non-WW3 data:
- Auto-computes omega and group velocity from grid parameters if needed
- Accepts pre-computed values for efficiency
- Works with any action density spectrum

Use Case 1 - Generic Data (auto-compute everything):
  converter = SpectrumConverter(action, depth=100.0, fr1=0.04, xfr=1.1, nk=30, nth=36)
  energy_2d, freq, dirs = converter.to_energy_2d()

Use Case 2 - WW3 Data (use pre-computed parameters):
  converter = SpectrumConverter(action, omega=omega_ww3, group_velocity=cg_ww3)
  energy_2d, freq, dirs = converter.to_energy_2d()

Use Case 3 - WW3 Optimized (use all pre-computed parameters):
  converter = SpectrumConverter(action, omega=omega_ww3, dintegral=dden_ww3,
                               group_velocity=cg_ww3, wavenumber=wn_ww3)
  params = converter.compute_wave_parameters()
"""

import numpy as np
from .constants import GRAV, TPI, SMALL
from .dispersion import solve_dispersion
from . import spectrum_converter as sc


class SpectrumConverter:
    """
    Unified spectrum converter for WW3 and non-WW3 action density data.

    Intelligently computes required parameters based on what's provided,
    while allowing pre-computed values for efficiency.
    """

    def __init__(self, action, depth=None,
                 # Grid parameters for auto-generation
                 fr1=None, xfr=None, nk=None, nth=None,
                 # Pre-computed WW3 parameters
                 omega=None, group_velocity=None,
                 dintegral=None, wavenumber=None,
                 directions=None,
                 # Convenience: pass all WW3 params in a dict
                 ww3_params=None):
        """
        Initialize unified spectrum converter.

        Parameters can be provided in order of computational efficiency:

        Option 1 (Most efficient - all WW3 pre-computed):
          action, depth, omega, group_velocity, dintegral, wavenumber

        Option 2 (Moderate - frequency parameters):
          action, depth, omega, group_velocity

        Option 3 (Generic - grid parameters):
          action, depth, fr1, xfr, nk, nth

        Arguments:
            action (ndarray): Action density spectrum (ndir, nfreq) [m²·s·rad⁻¹]
            depth (float): Water depth [m] - required for cg computation if not provided
            fr1 (float): First frequency [Hz] - default 0.04 Hz if not provided
            xfr (float): Frequency increment - default 1.1 if not provided
            nk (int): Number of frequency bins - inferred from action shape if not provided
            nth (int): Number of directional bins - inferred from action shape if not provided
            omega (ndarray): Pre-computed angular frequencies [rad/s]
            group_velocity (ndarray): Pre-computed group velocity [m/s]
            dintegral (ndarray): Pre-computed DDEN factors
            wavenumber (ndarray): Pre-computed wavenumber [1/m]
            directions (ndarray): Direction grid [radians]

        How fr1 is determined (in order of precedence):
            1. Provided explicitly as parameter
            2. Inferred from omega array: fr1 = omega[0]/(2π) [Hz]
            3. Default value: 0.04 Hz (typical WW3 value)

        ww3_params Dictionary:
            If ww3_params is provided, it can contain:
            - 'dden' or 'dintegral': Pre-computed DDEN factors (DDEN = DTH × DSII × SIG)
            - 'wn' or 'wavenumber': Pre-computed wavenumber [1/m]
            - 'fte': Energy tail factor
            - 'fttr': Period tail factor (defaults to fte if not provided)
            - 'ftwl': Wavelength tail factor

            Example:
            converter = SpectrumConverter(
                action,
                omega=omega,
                group_velocity=cg,
                depth=depth,
                ww3_params={
                    'dden': dden,
                    'wn': wn,
                    'fte': fte,
                    'fttr': fttr,
                    'ftwl': ftwl
                }
            )
        """

        # Validate and store action
        self.action = np.atleast_2d(np.asarray(action, dtype=float))
        if self.action.shape[0] < self.action.shape[1]:
            self.action = self.action.T

        self.ndir, self.nfreq = self.action.shape
        self.depth = depth

        # Extract WW3 parameters from dict if provided
        # ww3_params keys can be: 'dden'/'dintegral', 'wn'/'wavenumber', 'fte', 'fttr', 'ftwl'
        if ww3_params is not None:
            # Extract DDEN (support both key names)
            if dintegral is None:
                dintegral = ww3_params.get('dden') if ww3_params.get('dden') is not None else ww3_params.get('dintegral')

            # Extract wavenumber (support both key names)
            if wavenumber is None:
                wavenumber = ww3_params.get('wn') if ww3_params.get('wn') is not None else ww3_params.get('wavenumber')

            # Store tail factors for later use
            self.tail_fte = ww3_params.get('fte')
            self.tail_fttr = ww3_params.get('fttr')
            self.tail_ftwl = ww3_params.get('ftwl')
        else:
            # Initialize tail factors as None (will be computed if needed)
            self.tail_fte = None
            self.tail_fttr = None
            self.tail_ftwl = None

        # Store or compute frequency grid
        if omega is not None:
            self.omega = np.atleast_1d(np.asarray(omega, dtype=float))
            if len(self.omega) != self.nfreq:
                raise ValueError(f"omega length {len(self.omega)} != nfreq {self.nfreq}")
            # Infer fr1 and xfr from omega
            if len(self.omega) > 1:
                # Estimate xfr from frequency ratio
                self.xfr = self.omega[1] / self.omega[0]
                # Estimate fr1 (first frequency in Hz)
                # omega[0] = 2π × fr1
                self.fr1 = self.omega[0] / TPI
            else:
                # Single frequency - use defaults
                self.fr1 = 0.04
                self.xfr = 1.1
        elif fr1 is not None and xfr is not None and nk is not None:
            if nk != self.nfreq:
                raise ValueError(f"nk {nk} != nfreq {self.nfreq}")
            # Build frequency grid
            sigma_start = fr1 * TPI / (xfr**2)
            self.omega = np.zeros(nk)
            for ik in range(nk):
                sigma_start *= xfr
                self.omega[ik] = sigma_start
            self.fr1 = fr1
            self.xfr = xfr
        else:
            # Default: Use typical WW3 values if minimal info provided
            # Use WW3 defaults: fr1=0.04 Hz, xfr=1.1
            DEFAULT_FR1 = 0.04
            DEFAULT_XFR = 1.1

            # nk comes from action shape (already determined)
            nk = self.nfreq

            self.fr1 = DEFAULT_FR1
            self.xfr = DEFAULT_XFR

            # Build frequency grid with defaults
            sigma_start = self.fr1 * TPI / (self.xfr**2)
            self.omega = np.zeros(nk)
            for ik in range(nk):
                sigma_start *= self.xfr
                self.omega[ik] = sigma_start

        # Store or compute group velocity
        if group_velocity is not None:
            self.cg = np.atleast_1d(np.asarray(group_velocity, dtype=float))
            if len(self.cg) != self.nfreq:
                raise ValueError(f"group_velocity length {len(self.cg)} != nfreq {self.nfreq}")
        else:
            if depth is None:
                raise ValueError("Must provide depth to compute group velocity")
            self.cg = np.zeros(self.nfreq)
            for i in range(self.nfreq):
                _, self.cg[i] = solve_dispersion(self.omega[i], depth)

        # Store or compute wavenumber
        if wavenumber is not None:
            self.wn = np.atleast_1d(np.asarray(wavenumber, dtype=float))
        else:
            if depth is None:
                raise ValueError("Must provide depth to compute wavenumber")
            self.wn = np.zeros(self.nfreq)
            for i in range(self.nfreq):
                self.wn[i], _ = solve_dispersion(self.omega[i], depth)

        # Store or compute direction grid
        if directions is not None:
            self.directions = np.atleast_1d(np.asarray(directions, dtype=float))
        else:
            if nth is not None and nth != self.ndir:
                raise ValueError(f"nth {nth} != ndir {self.ndir}")
            self.directions = np.linspace(0, 2*np.pi, self.ndir, endpoint=False)

        # Compute grid spacing
        self.dth = 2.0 * np.pi / self.ndir

        # Compute or use provided DSII
        if hasattr(self, 'xfr'):
            sxfr = 0.5 * (self.xfr - 1.0 / self.xfr)
            self.dsii = self.omega * sxfr
            self.dsii[0] = 0.5 * self.omega[0] * (self.xfr - 1.0)
            self.dsii[-1] = 0.5 * self.omega[-1] * (self.xfr - 1.0) / self.xfr
        else:
            # Estimate DSII from frequency spacing
            self.dsii = np.zeros(self.nfreq)
            self.dsii[0] = (self.omega[1] - self.omega[0]) / 2.0
            for i in range(1, self.nfreq - 1):
                self.dsii[i] = (self.omega[i+1] - self.omega[i-1]) / 2.0
            self.dsii[-1] = (self.omega[-1] - self.omega[-2]) / 2.0

        # Compute or use provided DDEN
        if dintegral is not None:
            self.dden = np.atleast_1d(np.asarray(dintegral, dtype=float))
        else:
            self.dden = self.dth * self.dsii * self.omega

        # Compute frequencies [Hz]
        self.frequencies = self.omega / TPI

    def to_energy_2d(self):
        """
        Convert 2D action density to 2D energy density.

        Formula: E(f,θ) = A(k,θ) × SIG × (2π / CG)

        Returns:
            tuple: (energy_2d, frequencies, directions)
        """
        # Delegate to existing spectrum_converter function
        energy_2d, freq, dirs, _ = sc.action_to_energy_2d(
            self.action, self.omega, self.cg,
            directions=self.directions, ddir=self.dth, dsii=self.dsii
        )
        return energy_2d, freq, dirs

    def to_frequency_spectrum_1d(self):
        """
        Get 1D frequency spectrum.

        Formula: E(f) = ∑_θ [A(θ,f) × (DDEN/CG)] × DTH

        Returns:
            tuple: (frequencies, e_freq)
        """
        # Delegate to existing spectrum_converter function
        freq, e_freq = sc.action_to_frequency_spectrum_1d(
            self.action, self.omega, self.cg,
            dintegral=self.dden, directions=self.directions, ddir=self.dth, dsii=self.dsii
        )
        return freq, e_freq

    def to_directional_spectrum_1d(self):
        """
        Get 1D directional spectrum.

        Formula: E(θ) = ∑_f [A(θ,f) × (DDEN/CG) × DSII]

        Returns:
            tuple: (directions, e_dir)
        """
        # Delegate to existing spectrum_converter function
        dirs, e_dir = sc.action_to_directional_spectrum_1d(
            self.action, self.omega, self.cg,
            dintegral=self.dden, directions=self.directions, ddir=self.dth, dsii=self.dsii
        )
        return dirs, e_dir

    def get_peak_frequency(self):
        """Get peak frequency from spectrum."""
        # Delegate to existing spectrum_converter function
        return sc.get_peak_frequency(
            self.action, self.omega, self.cg,
            dintegral=self.dden, directions=self.directions, ddir=self.dth, dsii=self.dsii
        )

    def get_peak_direction(self):
        """Get peak direction from spectrum."""
        # Delegate to existing spectrum_converter function
        return sc.get_peak_direction(
            self.action, self.omega, self.cg,
            dintegral=self.dden, directions=self.directions, ddir=self.dth, dsii=self.dsii
        )

    def normalize(self, target_energy=None):
        """Normalize 2D energy spectrum to target energy."""
        energy_2d, _, _ = self.to_energy_2d()
        total = np.sum(energy_2d)

        if total < SMALL:
            return energy_2d

        normalized = energy_2d / total
        if target_energy is not None:
            normalized = normalized * target_energy

        return normalized

    def compute_wave_parameters(self, fte=None, fttr=None, ftwl=None):
        """
        Compute wave parameters (Hs, Tp, mean direction, etc.).

        Requires depth to be set. Optional tail factors for PM extension.

        Arguments:
            fte (float): Energy tail factor (auto-computed if not provided)
            fttr (float): Period tail factor (defaults to fte)
            ftwl (float): Wavelength tail factor

        Returns:
            dict: Wave parameters

        Precedence for tail factors:
            1. Explicitly provided as arguments to this method
            2. Provided via ww3_params dict during initialization
            3. Auto-computed from grid parameters
        """
        if self.depth is None:
            raise ValueError("Depth required to compute wave parameters")

        # Determine tail factors with 3-tier precedence
        # Tier 1: Explicitly provided to this method (highest priority)
        # Tier 2: From ww3_params dict (if provided during init)
        # Tier 3: Auto-computed from grid (lowest priority/default)

        if fte is None:
            fte = self.tail_fte if self.tail_fte is not None else 0.25 * self.omega[-1] * self.dth * self.omega[-1]

        if fttr is None:
            fttr = self.tail_fttr if self.tail_fttr is not None else fte * 0.20  # Approximate

        if ftwl is None:
            ftwl = self.tail_ftwl if self.tail_ftwl is not None else (GRAV / 6.0) / self.omega[-1] * self.dth * self.omega[-1]

        # Initialize moments
        m0 = m1 = m2 = m_1 = 0.0
        mom_x = mom_y = mom_wn = 0.0

        # Integrate over spectrum
        factor = self.dden / (self.cg + SMALL)
        cos_dirs = np.cos(self.directions)
        sin_dirs = np.sin(self.directions)

        for ik in range(self.nfreq):
            action_band = np.sum(self.action[:, ik])
            ebd = action_band * factor[ik]

            m0 += ebd
            m1 += ebd * self.omega[ik]
            m2 += ebd * self.omega[ik]**2
            m_1 += ebd / (self.omega[ik] + SMALL)

            mom_x += np.sum(self.action[:, ik] * cos_dirs) * factor[ik]
            mom_y += np.sum(self.action[:, ik] * sin_dirs) * factor[ik]
            mom_wn += ebd / (self.wn[ik] + SMALL)

        # Add tail extension
        eband_tail = np.sum(self.action[:, -1]) / (self.cg[-1] + SMALL)
        m0 += fte * eband_tail
        m1 += fte * 0.20 * eband_tail
        m2 += fte * 0.5 * self.omega[-1]**4 * self.dth * eband_tail
        m_1 += fttr * eband_tail
        mom_wn += ftwl * eband_tail

        mom_x += fte * np.sum(self.action[:, -1] * cos_dirs) / (self.cg[-1] + SMALL)
        mom_y += fte * np.sum(self.action[:, -1] * sin_dirs) / (self.cg[-1] + SMALL)

        # Compute wave parameters
        params = {}

        # Wave height
        params['hs'] = 4.0 * np.sqrt(max(0.0, m0))

        # Periods
        if m0 > SMALL and m1 > SMALL and m2 > SMALL:
            params['t01'] = TPI * m0 / m1
            params['t02'] = TPI * np.sqrt(m0 / m2)
            params['t0m1'] = TPI * m_1 / m0 if m_1 > SMALL else 0.0
        else:
            params['t01'] = params['t02'] = params['t0m1'] = 0.0

        # Peak
        freq, e_freq = self.to_frequency_spectrum_1d()
        peak_idx = np.argmax(e_freq)
        params['tp'] = 1.0 / freq[peak_idx] if freq[peak_idx] > 0 else 0.0
        params['fp'] = freq[peak_idx]

        # Direction
        if (np.abs(mom_x) + np.abs(mom_y)) > SMALL:
            params['thm'] = np.arctan2(mom_y, mom_x)
        else:
            params['thm'] = 0.0

        # Directional spread
        if m0 > SMALL:
            dir_spread_arg = (mom_x**2 + mom_y**2) / (m0**2 + SMALL)
            dir_spread_arg = np.clip(dir_spread_arg, 0.0, 1.0)
            params['ths'] = np.sqrt(2.0 * (1.0 - np.sqrt(dir_spread_arg)))
        else:
            params['ths'] = 0.0

        # Mean wavelength
        params['wlm'] = TPI * mom_wn / m0 if m0 > SMALL else 0.0

        # Spectral width
        if m0 > SMALL and m1 > SMALL:
            width_arg = m2 * m0 / (m1**2 + SMALL) - 1.0
            params['width'] = np.sqrt(max(0.0, width_arg))
        else:
            params['width'] = 0.0

        # Store moments
        params['m0'] = m0
        params['m1'] = m1
        params['m2'] = m2
        params['m_1'] = m_1

        return params
