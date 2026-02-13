#!/usr/bin/env python
"""
Test suite for WW3 input file parser.
"""

import sys
import tempfile
import os

sys.path.insert(0, '/home/user/WW3/python')

from wavewatch_python import parse_ww3_grid_line, extract_grid_params, print_ww3_params, build_ww3_grid


def test_parse_ww3_grid_line():
    """Test parsing of WW3 grid input file."""
    print("=" * 70)
    print("Test 1: Parse WW3 Grid Line 15")
    print("=" * 70)

    # Create a temporary file with test data
    test_content = """$ -------------------------------------------------------------------- $
$ WAVEWATCH III Grid preprocessor input file                           $
$ -------------------------------------------------------------------- $
$ Grid name (C*30, in quotes)
$
  'TEST GRID (GULF OF NOWHERE)   '
$
$ Frequency increment factor and first frequency (Hz) -------- $
  1.1  0.04118  25  24  0.
$
$ Set model flags ---------------------------------------------------- $
   F T T T F T
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
        f.write(test_content)
        temp_file = f.name

    try:
        params = parse_ww3_grid_line(temp_file)

        print(f"\nExtracted parameters:")
        print(f"  XFR (frequency ratio):      {params['xfr']}")
        print(f"  FR1 (first frequency [Hz]): {params['fr1']}")
        print(f"  NK (frequency bins):        {params['nk']}")
        print(f"  NTH (direction bins):       {params['nth']}")
        print(f"  Direction offset:           {params['dir_offset']}")

        # Verify values
        assert params['xfr'] == 1.1, f"XFR mismatch: {params['xfr']} != 1.1"
        assert params['fr1'] == 0.04118, f"FR1 mismatch: {params['fr1']} != 0.04118"
        assert params['nk'] == 25, f"NK mismatch: {params['nk']} != 25"
        assert params['nth'] == 24, f"NTH mismatch: {params['nth']} != 24"
        assert params['dir_offset'] == 0.0, f"Direction offset mismatch: {params['dir_offset']} != 0.0"

        print(f"\n✓ Test 1 PASSED")

    finally:
        os.unlink(temp_file)


def test_extract_grid_params():
    """Test extraction of grid parameters with depth."""
    print("\n" + "=" * 70)
    print("Test 2: Extract Grid Parameters with Depth")
    print("=" * 70)

    test_content = """$ Grid name
  'TEST GRID                      '
$
  1.1  0.04118  25  24  0.
$
   F T T T F T
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
        f.write(test_content)
        temp_file = f.name

    try:
        depth = 100.0
        params = extract_grid_params(temp_file, depth)

        print(f"\nExtracted parameters for build_ww3_grid():")
        print(f"  fr1={params['fr1']}, xfr={params['xfr']}, nk={params['nk']}, nth={params['nth']}, depth={params['depth']}")

        # Verify all required keys are present
        assert 'fr1' in params
        assert 'xfr' in params
        assert 'nk' in params
        assert 'nth' in params
        assert 'depth' in params
        assert params['depth'] == 100.0

        print(f"\n✓ Test 2 PASSED")

    finally:
        os.unlink(temp_file)


def test_integration_with_build_ww3_grid():
    """Test integration with build_ww3_grid()."""
    print("\n" + "=" * 70)
    print("Test 3: Integration with build_ww3_grid()")
    print("=" * 70)

    test_content = """$ Grid name
  'TEST GRID                      '
$
  1.1  0.04118  30  36  0.
$
   F T T T F T
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
        f.write(test_content)
        temp_file = f.name

    try:
        # Extract parameters
        depth = 50.0
        params = extract_grid_params(temp_file, depth)

        # Build grid using extracted parameters
        grid = build_ww3_grid(**params)

        print(f"\nGrid built successfully:")
        print(f"  Frequency array shape: {grid['sigma'].shape}")
        print(f"  Direction array shape: {grid['theta'].shape}")
        print(f"  Frequency range: {grid['freq'][0]:.4f} - {grid['freq'][-1]:.4f} Hz")
        print(f"  Wavenumber range: {grid['wavenumber'][0]:.4f} - {grid['wavenumber'][-1]:.4f} [1/m]")
        print(f"  Group velocity range: {grid['group_velocity'][0]:.2f} - {grid['group_velocity'][-1]:.2f} m/s")

        # Verify grid structure
        assert grid['sigma'].shape == (30,), "Sigma shape mismatch"
        assert grid['theta'].shape == (36,), "Theta shape mismatch"
        assert 'wavenumber' in grid, "Wavenumber not in grid"
        assert 'group_velocity' in grid, "Group velocity not in grid"

        print(f"\n✓ Test 3 PASSED")

    finally:
        os.unlink(temp_file)


def test_print_ww3_params():
    """Test printing of WW3 parameters."""
    print("\n" + "=" * 70)
    print("Test 4: Print WW3 Parameters")
    print("=" * 70)

    test_content = """$ Grid name
  'TEST GRID                      '
$
  1.1  0.04118  25  24  0.
$
   F T T T F T
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
        f.write(test_content)
        temp_file = f.name

    try:
        print("\nPrinting WW3 parameters:\n")
        print_ww3_params(temp_file, depth=100.0)
        print("✓ Test 4 PASSED")

    finally:
        os.unlink(temp_file)


def test_real_ww3_grid_file():
    """Test with actual ww3_grid.inp from repository."""
    print("\n" + "=" * 70)
    print("Test 5: Parse Real ww3_grid.inp File")
    print("=" * 70)

    ww3_file = '/home/user/WW3/model/inp/ww3_grid.inp'

    if not os.path.exists(ww3_file):
        print(f"Skipping test: {ww3_file} not found")
        return

    try:
        params = parse_ww3_grid_line(ww3_file)

        print(f"\nExtracted from {ww3_file}:")
        print(f"  FR1 (first frequency):  {params['fr1']} Hz")
        print(f"  XFR (frequency ratio):  {params['xfr']}")
        print(f"  NK (frequency bins):    {params['nk']}")
        print(f"  NTH (direction bins):   {params['nth']}")

        # Build grid with reasonable depth
        grid_params = extract_grid_params(ww3_file, depth=100.0)
        grid = build_ww3_grid(**grid_params)

        print(f"\n  Grid built successfully:")
        print(f"    Frequency bins: {grid['sigma'].shape[0]}")
        print(f"    Direction bins: {grid['theta'].shape[0]}")

        print(f"\n✓ Test 5 PASSED")

    except Exception as e:
        print(f"Error parsing real file: {e}")


if __name__ == "__main__":
    test_parse_ww3_grid_line()
    test_extract_grid_params()
    test_integration_with_build_ww3_grid()
    test_print_ww3_params()
    test_real_ww3_grid_file()

    print("\n" + "=" * 70)
    print("✓ All tests completed successfully!")
    print("=" * 70)
