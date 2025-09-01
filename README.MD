
Polynomial Interpolation Solver
This program calculates the coefficients of a polynomial that passes through a given set of points. The points are provided in a JSON file, where the y-coordinates can be in various numerical bases. The program uses Lagrange interpolation to find the unique polynomial for the given points.

Prerequisites
Before you can compile and run this code, you need the following dependencies installed:

A C++ Compiler: A modern C++ compiler that supports C++17. g++ is recommended.

Boost C++ Libraries: The program uses the Boost Multiprecision library (cpp_int) to handle very large numbers.

json.hpp: The nlohmann/json single-header library for parsing JSON.

Installation
On Debian/Ubuntu (or WSL):

sudo apt update
sudo apt install -y g++ libboost-all-dev

json.hpp library:
Download the json.hpp file from the official repository and place it in the same directory as solve_poly.cpp.

File Structure
Your project folder should be organized as follows for the commands to work correctly:

/hashira/
├── solve_poly.cpp      # The main C++ source code
├── json.hpp            # The JSON library header
└── testcases.json      # Your JSON input file

How to Run
Follow these two steps to compile and execute the program.

1. Compile the Code
Open your terminal, navigate to the /hashira/ directory, and run the following command to compile the C++ source file.

g++ -std=c++17 solve_poly.cpp -o solver_poly

This command will create an executable file named solver_poly.

2. Execute the Program
Run the compiled program by passing the path to your JSON test case file as an argument:

./solver_poly testcases.json

Example
To run with the second test case file:

# Ensure testcases2.json is in the directory
./solver_poly testcases2.json

Output Format
The program will print the results directly to the terminal. The output consists of two parts:

Polynomial Coefficients: A list of numbers representing the coefficients of the polynomial, starting from the constant term (c_0) up to the highest degree term (c_k−1).

Secret (f(0)): The value of the polynomial at x=0, which is always equal to the constant term (c_0).
