# numintegrate

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

**numintegrate** is a Python library containing fundamental numerical integration methods. 

## Features

The library includes the following five integration methods:
- Composite Midpoint Rule (`midpoint`)
- Composite Trapezoidal Rule (`trapezoidal`)
- Composite Simpson's 1/3 Rule (`simpson13`)
- Composite Simpson's 3/8 Rule (`simpson38`)
- Gauss Quadrature (`gauss_quadrature`)

## Installation

You can install the package by cloning the repository and using `pip`:

```bash
git clone https://github.com/silanurcayir/numintegrate.git
cd numintegrate
pip install .
```

## Usage

```python
import numpy as np
from numintegrate import trapezoidal, simpson13, gauss_quadrature

# Define a function to integrate
def f(x):
    return np.sin(x)

# Calculate the integral of sin(x) from 0 to pi
a, b = 0, np.pi
n = 100 # Number of subintervals (must be even for Simpson's 1/3)

result_trapz = trapezoidal(f, a, b, n)
result_simp = simpson13(f, a, b, n)
result_gauss = gauss_quadrature(f, a, b, deg=5)

print(f"Trapezoidal: {result_trapz}")
print(f"Simpson's 1/3: {result_simp}")
print(f"Gauss Quadrature: {result_gauss}")
```

## Development and Testing
To run tests, make sure you have `pytest` installed and simply run:
```bash
pytest tests/
```

## Author
**Sılanur Çayır**
- GitHub: [@silanurcayir](https://github.com/silanurcayir)
- Email: [silanurcayir@gmail.com](mailto:silanurcayir@gmail.com)

## Contributing
Contributions are welcome! You can open an *issue* or submit a *pull request* directly.

## Citation
This library was developed as part of the following graduation project:

> Çayır, S. (2026). *Python Tabanlı Nümerik İntegrasyon Kütüphanesinin
> Geliştirilmesi.* Graduation Project. Bursa Technical University,
> Faculty of Engineering and Natural Sciences, Department of Mathematics. Advisor:
> Assoc. Prof. Dr. Burhan Alveroğlu.
