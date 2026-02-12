"""
Spectrum conversion utilities for WaveWatch III.

Converts between different spectral representations:
- Action density to energy density (frequency-direction)
- Action density to 1D frequency spectrum
- Action density to 1D directional spectrum

This implementation mirrors the Fortran w3iogomd.F90 and w3gridmd.F90 algorithms,
properly handling the coordinate transformation from action density A(k,θ) in
wavenumber-direction space to energy density E(f,θ) in frequency-direction space.
"""

import numpy as np
from .constants import GRAV, TPI, SMALL


def build_ww3_grid(fr1, xfr, nk, nth, depth=None):
    """
    Build WW3 spectral grid arrays (mirrors w3gridmd.F90).

    Constructs the frequency and directional grids with proper frequency
    bandwidth (DSII) computation including special cases for first and last bins.

    Arguments:
        fr1 (float): First frequency [Hz]
        xfr (float): Frequency increment factor (geometric ratio)
        nk (int): Number of frequency bins
        nth (int): Number of directional bins
        depth (float, optional): Water depth [m]. If provided, automatically
                                computes wavenumber and group velocity.

    Returns:
        dict with keys:
            'freq' : Frequency array [Hz]
            'sigma' : Angular frequency array SIG [rad/s]
            'theta' : Direction array [radians]
            'dth' : Directional bin width [radians]
            'dsii' : Frequency bandwidth array [rad/s]
            'dden' : DDEN = DTH × DSII × SIG (conversion factor)
            'fte' : Tail energy factor
            'wavenumber' : Wavenumber array [1/m] (only if depth provided)
            'group_velocity' : Group velocity array [m/s] (only if depth provided)
    """
    dth = TPI / nth

    # Directions: evenly spaced, centred
    theta = np.array([(ith - 0.5) * dth for ith in range(1, nth + 1)])

    # Frequencies: geometric series (mimics Fortran w3gridmd.F90)
    sigma_start = fr1 * TPI / (xfr**2)
    sigma = np.zeros(nk)
    for ik in range(nk):
        sigma_start *= xfr
        sigma[ik] = sigma_start

    freq = sigma / TPI

    # Frequency bandwidths DSII
    sxfr = 0.5 * (xfr - 1.0 / xfr)
    dsii = sigma * sxfr

    # Special cases for first and last bins (w3gridmd.F90 lines 1305-1306)
    dsii[0] = 0.5 * sigma[0] * (xfr - 1.0)
    dsii[-1] = 0.5 * sigma[-1] * (xfr - 1.0) / xfr

    # DDEN: combined factor for discrete integration
    dden = dth * dsii * sigma

    # Tail energy factor
    fte = 0.25 * sigma[-1] * dth * sigma[-1]

    result = {
        'freq': freq,
        'sigma': sigma,
        'theta': theta,
        'dth': dth,
        'dsii': dsii,
        'dden': dden,
        'fte': fte,
    }

    # Optionally compute wavenumber and group velocity if depth is provided
    if depth is not None:
        from .dispersion import solve_dispersion
        wn = np.zeros(nk)
        cg = np.zeros(nk)
        for ik in range(nk):
            wn[ik], cg[ik] = solve_dispersion(sigma[ik], depth)
        result['wavenumber'] = wn
        result['group_velocity'] = cg

    return result


def action_to_energy_2d(action, omega, group_velocity, dintegral=None,
                        directions=None, ddir=None, dsii=None):
    """
    Convert 2D action density spectrum to 2D energy density spectrum.

    Converts from action density A(k,θ) to energy density E(f,θ) using:
      E(f,θ) = A(k,θ) × σ × (∂k/∂f)
              = A(k,θ) × SIG × (2π / CG)

    Where:
      - A(k,θ) = F(k,θ) / σ (action = energy / intrinsic frequency)
      - σ = SIG = intrinsic (angular) frequency [rad/s]
      - ∂k/∂f = 2π / CG = Jacobian of (k,θ) → (f,θ) transformation
      - SIG is the angular frequency array from WW3 w3gridmd.F90

    This is a pure spectral conversion with NO bin-width factors included.
    Bin widths (DTH, DSII) are applied only when integrating to 1D spectra.

    Arguments:
        action (ndarray): Action density spectrum (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s], length nfreq
        group_velocity (ndarray): Group velocity [m/s], length nfreq
        dintegral (ndarray, optional): DDEN factors (not used in core conversion)
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]
        dsii (ndarray, optional): Frequency bandwidths [rad/s]

    Returns:
        tuple:
            - energy_2d (ndarray): Energy density (ndir, nfreq) [m²/Hz/rad]
            - frequencies (ndarray): Frequency array [Hz]
            - directions (ndarray): Direction array [radians]
            - grid_params (dict): Grid parameters including dth, dsii, dden
    """

    # Validate and reshape action
    action = np.atleast_2d(np.asarray(action, dtype=float))
    if action.shape[0] < action.shape[1]:
        action = action.T

    ndir, nfreq = action.shape
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    group_velocity = np.atleast_1d(np.asarray(group_velocity, dtype=float))

    # Validate dimensions
    if len(omega) != nfreq:
        raise ValueError(f"omega length ({len(omega)}) != nfreq ({nfreq})")
    if len(group_velocity) != nfreq:
        raise ValueError(f"group_velocity length ({len(group_velocity)}) != nfreq ({nfreq})")

    # Generate or use provided direction grid
    if directions is None:
        directions = np.linspace(0, 2*np.pi, ndir, endpoint=False)
    else:
        directions = np.atleast_1d(np.asarray(directions, dtype=float))

    # Compute or use provided directional spacing
    if ddir is None:
        ddir = 2.0 * np.pi / ndir

    # Compute or use provided DSII (frequency bandwidths)
    if dsii is None:
        sxfr = 0.5 * (1.1 - 1.0/1.1)  # Default XFR = 1.1
        dsii_array = omega * sxfr
    else:
        dsii_array = np.atleast_1d(np.asarray(dsii, dtype=float))

    # DDEN factor for later use
    dden = ddir * dsii_array * omega

    # CORRECT CONVERSION FORMULA (no bin widths):
    # E(f,θ) = A(k,θ) × σ × (2π / CG)
    # Where: σ = SIG = intrinsic (angular) frequency = omega
    # This combines the action-to-energy relationship (×σ) with the Jacobian (×2π/CG)
    conversion_factor = omega * (2.0 * np.pi) / (group_velocity + SMALL)

    # Apply conversion factor to all directions at each frequency
    energy_2d = action * conversion_factor[np.newaxis, :]

    # Convert omega to frequency
    frequencies = omega / (2.0 * np.pi)

    # Return grid parameters for use in 1D integration
    grid_params = {
        'dth': ddir,
        'dsii': dsii_array,
        'dden': dden,
    }

    return energy_2d, frequencies, directions, grid_params


def action_to_frequency_spectrum_1d(action, omega, group_velocity, dintegral=None,
                                     directions=None, ddir=None, dsii=None):
    """
    Convert 2D action density to 1D frequency spectrum.

    Integrates action density over all directions to produce a frequency spectrum E(f).
    Follows w3iogomd.F90 frequency spectrum computation.

    The conversion is:
      E(f) = ∑_θ E(f,θ) × Δθ = ∑_θ [A(θ,f) × SIG × (2π/CG)] × DTH

    Arguments:
        action (ndarray): Action density (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        dintegral (ndarray, optional): DDEN factors (w3gridmd DDEN array)
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]
        dsii (ndarray, optional): Frequency bandwidths [rad/s]

    Returns:
        tuple:
            - freq (ndarray): Frequency array [Hz]
            - e_freq (ndarray): 1D frequency spectrum [m²/Hz]
    """

    # Get 2D energy spectrum and grid parameters
    energy_2d, freq, dirs, grid_params = action_to_energy_2d(
        action, omega, group_velocity, dintegral, directions, ddir, dsii
    )

    dth = grid_params['dth']

    # Integrate over directions with directional bin width
    # e_freq = ∑_θ E(f,θ) × DTH
    e_freq = np.sum(energy_2d, axis=0) * dth

    return freq, e_freq


def action_to_directional_spectrum_1d(action, omega, group_velocity, dintegral=None,
                                       directions=None, ddir=None, dsii=None):
    """
    Convert 2D action density to 1D directional spectrum.

    Integrates action density over all frequencies to produce a directional spectrum E(θ).
    Follows w3iogomd.F90 directional spectrum computation.

    The conversion is:
      E(θ) = ∑_f [A(θ,f) × SIG × (2π/CG) × DSII] = ∑_f [A(θ,f) × DSII × SIG / CG]

    Arguments:
        action (ndarray): Action density (ndir, nfreq) [m²·s·rad⁻¹]
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        dintegral (ndarray, optional): DDEN factors (w3gridmd DDEN array)
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]
        dsii (ndarray, optional): Frequency bandwidths [rad/s]

    Returns:
        tuple:
            - dirs (ndarray): Direction array [radians]
            - e_dir (ndarray): 1D directional spectrum [m²/rad]
    """

    # Validate and reshape action
    action = np.atleast_2d(np.asarray(action, dtype=float))
    if action.shape[0] < action.shape[1]:
        action = action.T

    ndir, nfreq = action.shape
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    group_velocity = np.atleast_1d(np.asarray(group_velocity, dtype=float))

    # Validate dimensions
    if len(omega) != nfreq:
        raise ValueError(f"omega length ({len(omega)}) != nfreq ({nfreq})")
    if len(group_velocity) != nfreq:
        raise ValueError(f"group_velocity length ({len(group_velocity)}) != nfreq ({nfreq})")

    # Generate or use provided direction grid
    if directions is None:
        directions = np.linspace(0, 2*np.pi, ndir, endpoint=False)
    else:
        directions = np.atleast_1d(np.asarray(directions, dtype=float))

    # Compute or use provided DSII (frequency bandwidths)
    if dsii is None:
        sxfr = 0.5 * (1.1 - 1.0/1.1)  # Default XFR = 1.1
        dsii_array = omega * sxfr
    else:
        dsii_array = np.atleast_1d(np.asarray(dsii, dtype=float))

    # CORRECT CONVERSION FORMULA with frequency bin width:
    # E(θ) = ∑_f [A(θ,f) × SIG × (2π/CG) × DSII]
    conversion_factor = omega * (2.0 * np.pi) / (group_velocity + SMALL)

    # Integrate over frequencies with frequency bin width
    # e_dir = ∑_f [action * conversion_factor * DSII]
    e_dir = np.sum(action * conversion_factor[np.newaxis, :] * dsii_array[np.newaxis, :], axis=1)

    return directions, e_dir


def normalize_spectrum(energy_2d, target_energy=None):
    """
    Normalize 2D energy spectrum to specific total energy.

    Arguments:
        energy_2d (ndarray): 2D energy spectrum (ndir, nfreq)
        target_energy (float, optional): Target total energy [m²]
                                        If None, returns normalized spectrum

    Returns:
        ndarray: Normalized spectrum
    """

    total_energy = np.sum(energy_2d)

    if total_energy < SMALL:
        return energy_2d

    normalized = energy_2d / total_energy

    if target_energy is not None:
        normalized = normalized * target_energy

    return normalized


def get_peak_frequency(action_2d, omega, group_velocity, directions=None, ddir=None):
    """
    Get peak frequency from action density spectrum.

    Finds the frequency with maximum energy in the frequency spectrum.

    Arguments:
        action_2d (ndarray): 2D action spectrum (ndir, nfreq)
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]

    Returns:
        float: Peak frequency [Hz]
    """

    # Get 1D frequency spectrum
    freq, e_freq = action_to_frequency_spectrum_1d(
        action_2d, omega, group_velocity, directions=directions, ddir=ddir
    )

    # Find peak
    peak_idx = np.argmax(e_freq)

    return freq[peak_idx]


def get_peak_direction(action_2d, omega, group_velocity, directions=None, ddir=None):
    """
    Get peak direction from action density spectrum.

    Finds the direction with maximum energy in the directional spectrum.

    Arguments:
        action_2d (ndarray): 2D action spectrum (ndir, nfreq)
        omega (ndarray): Angular frequency [rad/s]
        group_velocity (ndarray): Group velocity [m/s]
        directions (ndarray, optional): Direction grid [radians]
        ddir (float, optional): Directional spacing [radians]

    Returns:
        float: Peak direction [radians]
    """

    # Get 1D directional spectrum
    dirs, e_dir = action_to_directional_spectrum_1d(
        action_2d, omega, group_velocity, directions=directions, ddir=ddir
    )

    # Find peak
    peak_idx = np.argmax(e_dir)

    return dirs[peak_idx]


class SpectrumConverter:
    """
    Unified spectrum converter for WW3 and non-WW3 action density data.

    Intelligently computes required parameters based on what's provided,
    while allowing pre-computed values for efficiency.

    Supports three usage patterns:
    1. Generic: auto-compute everything from grid parameters
    2. WW3: provide pre-computed omega and group velocity
    3. Optimized: provide all parameters for maximum efficiency
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
        # Import here to avoid circular dependency
        from .dispersion import solve_dispersion

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
        energy_2d, freq, dirs, _ = action_to_energy_2d(
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
        freq, e_freq = action_to_frequency_spectrum_1d(
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
        dirs, e_dir = action_to_directional_spectrum_1d(
            self.action, self.omega, self.cg,
            dintegral=self.dden, directions=self.directions, ddir=self.dth, dsii=self.dsii
        )
        return dirs, e_dir

    def get_peak_frequency(self):
        """Get peak frequency from spectrum."""
        # Delegate to existing spectrum_converter function
        return get_peak_frequency(
            self.action, self.omega, self.cg,
            dintegral=self.dden, directions=self.directions, ddir=self.dth, dsii=self.dsii
        )

    def get_peak_direction(self):
        """Get peak direction from spectrum."""
        # Delegate to existing spectrum_converter function
        return get_peak_direction(
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

