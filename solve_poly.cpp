#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <boost/multiprecision/cpp_int.hpp>
#include "json.hpp" // Download from https://github.com/nlohmann/json

using namespace std;
using namespace boost::multiprecision;
using json = nlohmann::json;

using BigInt = cpp_int;

// Convert string in base `base` to BigInt
BigInt convertToDecimal(const string &s, int base) {
    BigInt result = 0;
    for (char c : s) {
        int digit;
        if (isdigit(c)) {
            digit = c - '0';
        } else if (isalpha(c)) {
            digit = 10 + (tolower(c) - 'a');
        } else {
            continue;
        }
        result = result * base + digit;
    }
    return result;
}

// Lagrange interpolation: returns coefficients of polynomial degree (k-1)
// The coefficients are returned in order [c0, c1, c2, ...]
vector<BigInt> interpolate(const vector<pair<BigInt, BigInt>> &points, int k) {
    int m = k - 1;
    vector<BigInt> final_coeffs(m + 1, 0);

    for (int i = 0; i < k; i++) {
        BigInt xi = points[i].first;
        BigInt yi = points[i].second;

        // Start with the numerator polynomial L_i(x) = 1
        vector<BigInt> term_coeffs = {1};
        BigInt denom = 1;

        for (int j = 0; j < k; j++) {
            if (i == j) continue;
            BigInt xj = points[j].first;

            // Multiply term_coeffs by (x - xj)
            // A polynomial [c0, c1] * (x - xj) becomes [-xj*c0, c0 - xj*c1, c1]
            vector<BigInt> new_term(term_coeffs.size() + 1, 0);
            for (size_t a = 0; a < term_coeffs.size(); a++) {
                new_term[a] += -xj * term_coeffs[a];
                new_term[a + 1] += term_coeffs[a];
            }
            term_coeffs = new_term;
            denom *= (xi - xj);
        }

        // Add the scaled term polynomial to the final result
        // final_coeffs += (yi / denom) * term_coeffs
        for (size_t a = 0; a < term_coeffs.size(); a++) {
            final_coeffs[a] += term_coeffs[a] * yi / denom;
        }
    }

    return final_coeffs;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: ./solver_poly <json_file>\n";
        return 1;
    }

    // Read JSON file
    ifstream inFile(argv[1]);
    if (!inFile) {
        cerr << "Error: Cannot open file " << argv[1] << "\n";
        return 1;
    }

    json data;
    inFile >> data;

    int k = data["keys"]["k"];

    vector<pair<BigInt, BigInt>> points;

    // --- CORRECTED CODE BLOCK ---
    // Explicitly fetch points for x=1, 2, ..., k to ensure correct numerical order.
    for (int i = 1; i <= k; ++i) {
        string key = to_string(i);
        BigInt x = i;

        // Check if the key exists before accessing
        if (data.find(key) == data.end()) {
            cerr << "Error: JSON key '" << key << "' not found.\n";
            return 1;
        }

        // Get the data for the specific key
        int base = stoi(data[key]["base"].get<string>());
        string valueStr = data[key]["value"];
        BigInt y = convertToDecimal(valueStr, base);

        points.push_back({x, y});
    }
    // --- END OF CORRECTION ---


    // Now `points` has exactly k elements in the correct order.
    // We can pass it directly to the interpolation function.
    vector<BigInt> coeff = interpolate(points, k);

    cout << "Polynomial coefficients (c0..c" << k - 1 << "):\n";
    for (size_t i = 0; i < coeff.size(); i++) {
        cout << coeff[i] << (i + 1 == coeff.size() ? "" : " ");
    }
    cout << "\n";

    // Optional: secret = f(0) = c0
    cout << "Secret (f(0)) = " << coeff[0] << "\n";

    return 0;
}