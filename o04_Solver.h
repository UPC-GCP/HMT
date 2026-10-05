#ifndef SOLVER_H_
#define SOLVER_H_

// Imports
#include <math.h>
#include <cstddef>
#include <complex.h>
#include <json/json.h>

// Class
template <size_t Dim> class Solver {
private:
    
public:
    // Variables
    double tolNum{}, lastRes{};
    size_t maxIter{}, lastIter{};
};

// Functions

inline double calcErr(std::vector<double> matA, std::vector<double> matB) {
    // Control
    double rsNew{}, errVal{};

    // Error
    for (size_t k = 0; k < matA.size(); k++) {
        errVal = abs(matA[k] - matB[k]); rsNew += errVal * errVal;
    }

    return std::sqrt(rsNew);
}


inline double calcErr(std::vector<double> matA, std::vector<double> matB, std::vector<bool> bObs) {
    // Control
    double rsNew{}, errVal{};

    // Error
    for (size_t k = 0; k < matA.size(); k++) {
        if (bObs[k]) { continue; }
        errVal = abs(matA[k] - matB[k]); rsNew += errVal * errVal;
    }

    return std::sqrt(rsNew);
}

inline double calcErr(std::vector<std::complex<double>> matA, std::vector<std::complex<double>> matB) {
    // Control
    double rsNew{}; std::complex<double> errVal{};

    // Error
    for (size_t k = 0; k < matA.size(); k++) {
        errVal = abs(matA[k] - matB[k]); rsNew += std::real(errVal * errVal);
    }

    return std::sqrt(rsNew);
}

#endif
