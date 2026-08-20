import numpy as np
from typing import Callable, Union

def _validate_inputs(f: Callable, a: Union[int, float], b: Union[int, float], n: int = None):
    """Common input validation."""
    if not callable(f):
        raise TypeError("The parameter 'f' must be a callable function.")
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        raise TypeError("Integration limits 'a' and 'b' must be real numbers.")
    if n is not None:
        if not isinstance(n, int) or n <= 0:
            raise ValueError("Number of subintervals 'n' must be a positive integer.")
    if a == b:
        return 0.0
    return None

def midpoint(f: Callable[[Union[float, np.ndarray]], Union[float, np.ndarray]], a: float, b: float, n: int) -> float:
    """
    Evaluate a definite integral using the composite midpoint rule.

    Args:
        f (Callable): The function to integrate.
        a (float): The lower limit of integration.
        b (float): The upper limit of integration.
        n (int): The number of subintervals.

    Returns:
        float: The approximate value of the integral.
    """
    early_return = _validate_inputs(f, a, b, n)
    if early_return is not None:
        return early_return

    x = np.linspace(a, b, n + 1)
    midpoints = (x[:-1] + x[1:]) / 2
    return float(np.sum(f(midpoints)) * (b - a) / n)

def trapezoidal(f: Callable[[Union[float, np.ndarray]], Union[float, np.ndarray]], a: float, b: float, n: int) -> float:
    """
    Evaluate a definite integral using the composite trapezoidal rule.

    Args:
        f (Callable): The function to integrate.
        a (float): The lower limit of integration.
        b (float): The upper limit of integration.
        n (int): The number of subintervals.

    Returns:
        float: The approximate value of the integral.
    """
    early_return = _validate_inputs(f, a, b, n)
    if early_return is not None:
        return early_return

    x = np.linspace(a, b, n + 1)
    y = f(x)
    return float((b - a) / (2 * n) * (y[0] + 2 * np.sum(y[1:-1]) + y[-1]))

def simpson13(f: Callable[[Union[float, np.ndarray]], Union[float, np.ndarray]], a: float, b: float, n: int) -> float:
    """
    Evaluate a definite integral using Simpson's 1/3 rule.

    Args:
        f (Callable): The function to integrate.
        a (float): The lower limit of integration.
        b (float): The upper limit of integration.
        n (int): The number of subintervals (must be even).

    Returns:
        float: The approximate value of the integral.
        
    Raises:
        ValueError: If `n` is not an even number.
    """
    early_return = _validate_inputs(f, a, b, n)
    if early_return is not None:
        return early_return

    if n % 2 != 0:
        raise ValueError("n must be an even integer for Simpson's 1/3 rule.")
    
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return float((b - a) / (3 * n) * (y[0] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2]) + y[-1]))

def simpson38(f: Callable[[Union[float, np.ndarray]], Union[float, np.ndarray]], a: float, b: float, n: int) -> float:
    """
    Evaluate a definite integral using Simpson's 3/8 rule.

    Args:
        f (Callable): The function to integrate.
        a (float): The lower limit of integration.
        b (float): The upper limit of integration.
        n (int): The number of subintervals (must be a multiple of 3).

    Returns:
        float: The approximate value of the integral.
        
    Raises:
        ValueError: If `n` is not a multiple of 3.
    """
    early_return = _validate_inputs(f, a, b, n)
    if early_return is not None:
        return early_return

    if n % 3 != 0:
        raise ValueError("n must be a multiple of 3 for Simpson's 3/8 rule.")
        
    x = np.linspace(a, b, n + 1)
    y = f(x)
    
    result = y[0] + y[-1]
    for i in range(1, n):
        if i % 3 == 0:
            result += 2 * y[i]
        else:
            result += 3 * y[i]
            
    return float((b - a) * 3 / (8 * n) * result)

def gauss_quadrature(f: Callable[[Union[float, np.ndarray]], Union[float, np.ndarray]], a: float, b: float, deg: int = 5) -> float:
    """
    Evaluate a definite integral using Gauss-Legendre quadrature.

    Args:
        f (Callable): The function to integrate.
        a (float): The lower limit of integration.
        b (float): The upper limit of integration.
        deg (int): The number of sample points and weights. Defaults to 5.

    Returns:
        float: The approximate value of the integral.
    """
    early_return = _validate_inputs(f, a, b)
    if early_return is not None:
        return early_return

    if not isinstance(deg, int) or deg <= 0:
        raise ValueError("Degree 'deg' must be a positive integer.")

    x, w = np.polynomial.legendre.leggauss(deg)
    # Transform roots from [-1, 1] to [a, b]
    t = 0.5 * (x + 1) * (b - a) + a
    return float(0.5 * (b - a) * np.sum(w * f(t)))