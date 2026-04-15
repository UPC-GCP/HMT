
#include <vector>
#include <cmath>

#include "Discretizer.h"


inline std::vector<double> newProdMatVec(std::vector<Matrix> Mat, std::vector<std::vector<double>> Vec){

    // Control
    int n = Vec.size(), m = Vec[0].size(), k; std::vector<double> aVec(n*m);

    // Calculate
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){

            // Control
            k = i * m + j;
            if ((i == 0 || i == n-1) && (j == 0 || j == m-1)){continue;}

            // Calculate
            if (i > 0 && i < n-1 && j > 0 && j < m-1){
                // Interior Nodes
                aVec[k] = Mat[k].aw * Vec[i-1][j] + Mat[k].ae * Vec[i+1][j] + Mat[k].as * Vec[i][j-1] + Mat[k].an * Vec[i][j+1] + Mat[k].ap * Vec[i][j];
            } else if (i == 0){
                // West Boundary - ae, ap
                aVec[k] = Mat[k].ae * Vec[i+1][j] + Mat[k].ap * Vec[i][j];
            } else if (i == n-1){
                // East Boundary - aw, ap
                aVec[k] = Mat[k].aw * Vec[i-1][j] + Mat[k].ap * Vec[i][j];
            } else if (j == 0){
                // South Boundary - an, ap
                aVec[k] = Mat[k].an * Vec[i][j+1] + Mat[k].ap * Vec[i][j];
            } else if (j == m-1){
                // North Boundary - as, ap
                aVec[k] = Mat[k].as * Vec[i][j-1] + Mat[k].ap * Vec[i][j];
            }

        }
    }

    return aVec;

}


inline std::vector<double> operProdMatVec(std::vector<Matrix> Mat, std::vector<double> Vec, int n, int m){

    // Control
    int k; std::vector<double> aVec(n*m);

    // Calculate
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){

            // Control
            k = i*m + j;
            if ((i == 0 || i == n-1) && (j == 0 || j == m-1)){continue;}

            // Calculate
            if (i > 0 && i < n-1 && j > 0 && j < m-1){
                // Interior Nodes
                aVec[k] = Mat[k].aw * Vec[k-m] + Mat[k].ae * Vec[k+m] + Mat[k].as * Vec[k-1] + Mat[k].an * Vec[k+1] + Mat[k].ap * Vec[k];
            } else if (i == 0){
                // West Boundary - ae, ap
                aVec[k] = Mat[k].ae * Vec[k+m] + Mat[k].ap * Vec[k];
            } else if (i == n-1){
                // East Boundary - aw, ap
                aVec[k] = Mat[k].aw * Vec[k-m] + Mat[k].ap * Vec[k];
            } else if (j == 0){
                // South Boundary - an, ap
                aVec[k] = Mat[k].an * Vec[k+1] + Mat[k].ap * Vec[k];
            } else if (j == m-1){
                // North Boundary - as, ap
                aVec[k] = Mat[k].as * Vec[k-1] + Mat[k].ap * Vec[k];
            }

        }
    }

    return aVec;

}

inline double operDotProd(std::vector<double> v1, std::vector<double> v2){

    // Calculate
    double sum = 0;
    for (size_t i = 0; i < v1.size(); i++){
        sum += v1[i] * v2[i];
    }
    
    return sum;

}

inline std::vector<double> operElementProd(std::vector<double> v1, std::vector<double> v2){

    // Calculate
    std::vector<double> vec(v1.size());
    for (size_t i = 0; i < v1.size(); i++){
        vec[i] = v1[i] * v2[i];
    }

    return vec;
}

inline std::vector<double> operCombLinVec(std::vector<double> v1, std::vector<double> v2, double a = 1, double b = 1){

    // Calculate
    std::vector<double> vec(v1.size());
    for (size_t i = 0; i < v1.size(); i++){
        vec[i] = a * v1[i] + b * v2[i];
    }

    return vec;

}
