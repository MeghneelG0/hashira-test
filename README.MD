# Polynomial Interpolation Solver

This program calculates the coefficients of a polynomial that passes through a given set of points. The points are provided in a JSON file, where the y-coordinates can be in various numerical bases. The program uses **Lagrange interpolation** to find the unique polynomial for the given points.

---

## Prerequisites

Before you can compile and run this code, ensure the following dependencies are installed:

- **C++ Compiler**: A modern C++ compiler that supports C++17 (e.g., `g++`).
- **Boost C++ Libraries**: The program uses the Boost Multiprecision library (`cpp_int`) to handle very large numbers.
- **json.hpp**: The [nlohmann/json](https://github.com/nlohmann/json) single-header library for parsing JSON.

---

## Installation

### On Debian/Ubuntu (or WSL)
```bash
sudo apt update
sudo apt install -y g++ libboost-all-dev

output: 
testcase1:
Polynomial coefficients (c0..c2):
3 0 1
Secret (f(0)) = 3

testcas2: 
Polynomial coefficients (c0..c6):
-6290016743746469796 13164366091085008239 -9806825133911860179 3544566661514267415 -673833059248572449 65027537304912471 -2290267902684209
Secret (f(0)) = -6290016743746469796
