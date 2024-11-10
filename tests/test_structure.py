import os
from pathlib import Path

import pytest

from complexbuilder.common.structure import find_complexes_from_iptm


def test_find_complexes_from_iptm():
    """Test find_complexes_from_iptm function."""
    testdir = Path("/data2/moriwaki/BGCcomplex/BGC0001818")
    complexes = find_complexes_from_iptm(testdir, 0.5)
    positive_control_complex_name = "BAV57458.1_BAV57459.1"
    predicted_complexes_name = [complex.name for complex in complexes]
    assert positive_control_complex_name in predicted_complexes_name
