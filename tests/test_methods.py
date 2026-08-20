import pytest
import numpy as np
from numintegrate.methods import midpoint, trapezoidal, simpson13, simpson38, gauss_quadrature

def test_midpoint():
    f = lambda x: 3 * x**2
    assert np.isclose(midpoint(f, 0, 1, 1000), 1.0, atol=1e-4)

def test_trapezoidal():
    f = lambda x: 3 * x**2
    assert np.isclose(trapezoidal(f, 0, 1, 1000), 1.0, atol=1e-4)

def test_simpson13():
    f = lambda x: x**3
    # Simpson's 1/3 is exact for polynomials up to degree 3
    assert np.isclose(simpson13(f, 0, 1, 2), 0.25, atol=1e-10)

def test_simpson13_odd_n():
    f = lambda x: x
    with pytest.raises(ValueError):
        simpson13(f, 0, 1, 3)

def test_simpson38():
    f = lambda x: x**3
    # Simpson's 3/8 is exact for polynomials up to degree 3
    assert np.isclose(simpson38(f, 0, 1, 3), 0.25, atol=1e-10)

def test_simpson38_invalid_n():
    f = lambda x: x
    with pytest.raises(ValueError):
        simpson38(f, 0, 1, 4)

def test_gauss_quadrature():
    f = lambda x: x**5
    # 3-point Gauss-Legendre is exact for polynomials up to degree 2n-1 = 5
    assert np.isclose(gauss_quadrature(f, 0, 1, deg=3), 1/6, atol=1e-10)
