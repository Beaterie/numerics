//
//  LinearEquationSystem.hpp
//  advanced-numerics
//
//  Created by Josef Roth on 18.08.20.
//  Copyright © 2020 Josef Roth. All rights reserved.
//

#ifndef LinearEquationSystem_hpp
#define LinearEquationSystem_hpp

#include <stdio.h>
#include "configuration.h"

using namespace std;
using namespace configuration;

class LinearEquationSystem {

public:
    // Members
    // Ax = b
    double m_A[matrixSize][matrixSize]; // A
    double m_b[matrixSize]; // b
    double m_sol_A[matrixSize]; // x (solution for Gauss)
    
    // x = Bx + c
    double m_B[matrixSize][matrixSize]; // B (for SOR method)
    double m_c[matrixSize]; // c
    double m_sol_B[matrixSize]; // x (solution for SOR)
    double m_B_2[matrixSize][matrixSize]; // B_2 (upper triangular matrix)
    double m_epsilon; // epsilon (for SOR)
    
    
    // Constructors
    LinearEquationSystem();
    LinearEquationSystem(double mat[matrixSize][matrixSize], double vec[matrixSize]);
 
    
    // Functions
    void gauss(bool withPivotElement);
    void solve_gauss();
    void SOR(double accuracy, double omega);
    void solve_SOR(double omega);
    
    bool convergence_condition(double accuracy, double omega);
    void swap_rows(int index_a, int index_b);
    void print_Ax_b();
    void print_Bx_c();
    void print_solution();
    
};

#endif
