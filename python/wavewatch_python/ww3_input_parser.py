"""
Parser for WW3 grid input file (ww3_grid.inp).

Extracts spectral parameters from line 15 of the input file for use with
SpectrumConverter and build_ww3_grid().

Reference: ww3_grid.inp documentation
"""

import re
from typing import Dict, Tuple, Optional


def parse_ww3_grid_line(filename: str) -> Dict[str, float]:
    """
    Parse line 15 of ww3_grid.inp to extract spectral parameters.

    Line 15 format:
        XFR  FR1       NK  NTH  DIR_OFFSET
        1.1  0.04118   25  24   0.

    Arguments:
        filename (str): Path to ww3_grid.inp file

    Returns:
        dict: Dictionary with keys 'xfr', 'fr1', 'nk', 'nth', 'dir_offset'

    Raises:
        ValueError: If line 15 cannot be parsed or has incorrect format
        FileNotFoundError: If file does not exist

    Example:
        >>> params = parse_ww3_grid_line('ww3_grid.inp')
        >>> print(params)
        {'xfr': 1.1, 'fr1': 0.04118, 'nk': 25, 'nth': 24, 'dir_offset': 0.0}

        >>> from wavewatch_python import build_ww3_grid
        >>> grid = build_ww3_grid(
        ...     fr1=params['fr1'],
        ...     xfr=params['xfr'],
        ...     nk=int(params['nk']),
        ...     nth=int(params['nth']),
        ...     depth=100.0
        ... )
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the spectral parameters line
    # It should be around line 15 (after frequency/direction comments)
    spectral_line = None
    line_number = None

    for i, line in enumerate(lines):
        # Skip comment lines (start with $)
        if line.strip().startswith('$'):
            continue

        # Skip empty lines
        if not line.strip():
            continue

        # Look for line containing frequency parameters
        # Should have 5 values: XFR FR1 NK NTH DIR_OFFSET
        values = line.split()

        # The spectral line should have exactly 5 numeric values
        if len(values) == 5:
            try:
                float_values = [float(v) for v in values]
                spectral_line = values
                line_number = i + 1
                break
            except ValueError:
                # Not numeric values, skip
                continue

    if spectral_line is None:
        raise ValueError(
            f"Could not find spectral parameters line in {filename}. "
            "Expected line with 5 numeric values: XFR FR1 NK NTH DIR_OFFSET"
        )

    try:
        xfr, fr1, nk, nth, dir_offset = [float(v) for v in spectral_line]
    except ValueError as e:
        raise ValueError(
            f"Failed to parse spectral parameters on line {line_number}: {spectral_line}\n"
            f"Expected: XFR FR1 NK NTH DIR_OFFSET (all numeric)"
        ) from e

    return {
        'xfr': xfr,
        'fr1': fr1,
        'nk': int(nk),
        'nth': int(nth),
        'dir_offset': dir_offset,
    }


def extract_grid_params(filename: str, depth: float) -> Dict:
    """
    Extract WW3 grid parameters and prepare for SpectrumConverter.

    Convenience function that parses ww3_grid.inp and returns parameters
    ready to use with build_ww3_grid().

    Arguments:
        filename (str): Path to ww3_grid.inp file
        depth (float): Water depth [m] (not in input file, must be provided separately)

    Returns:
        dict: Dictionary with 'fr1', 'xfr', 'nk', 'nth', 'depth' ready for build_ww3_grid()

    Example:
        >>> params = extract_grid_params('ww3_grid.inp', depth=100.0)
        >>> from wavewatch_python import build_ww3_grid
        >>> grid = build_ww3_grid(**params)
    """
    parsed = parse_ww3_grid_line(filename)
    return {
        'fr1': parsed['fr1'],
        'xfr': parsed['xfr'],
        'nk': parsed['nk'],
        'nth': parsed['nth'],
        'depth': depth,
    }


def print_ww3_params(filename: str, depth: float = None) -> None:
    """
    Print extracted WW3 grid parameters in readable format.

    Useful for verification and debugging.

    Arguments:
        filename (str): Path to ww3_grid.inp file
        depth (float, optional): Water depth [m] (for reference)

    Example:
        >>> print_ww3_params('ww3_grid.inp', depth=100.0)
        WW3 Grid Parameters
        ===================
        Frequency ratio (XFR):    1.1
        First frequency (FR1):    0.04118 Hz
        Number of frequencies (NK):  25
        Number of directions (NTH):  36
        Direction offset:         0.0
        Water depth:              100.0 m
    """
    params = parse_ww3_grid_line(filename)

    print("WW3 Grid Parameters")
    print("=" * 40)
    print(f"Frequency ratio (XFR):      {params['xfr']:.4f}")
    print(f"First frequency (FR1):      {params['fr1']:.5f} Hz")
    print(f"Number of frequencies (NK): {params['nk']}")
    print(f"Number of directions (NTH): {params['nth']}")
    print(f"Direction offset:           {params['dir_offset']:.1f}")
    if depth is not None:
        print(f"Water depth:                {depth:.1f} m")
    print()


# Example usage
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python ww3_input_parser.py <ww3_grid.inp> [depth]")
        print()
        print("Example:")
        print("  python ww3_input_parser.py ww3_grid.inp 100.0")
        sys.exit(1)

    input_file = sys.argv[1]
    depth = float(sys.argv[2]) if len(sys.argv) > 2 else None

    try:
        print_ww3_params(input_file, depth)
        params = parse_ww3_grid_line(input_file)
        print("Parsed parameters (copy-paste ready):")
        print(f"  fr1={params['fr1']}, xfr={params['xfr']}, nk={params['nk']}, nth={params['nth']}")
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)
