//
//  LinearEquationSystem.cpp
//  advanced-numerics
//
//  Created by Josef Roth on 18.08.20.
//  Copyright © 2020 Josef Roth. All rights reserved.
//

#include "LinearEquationSystem.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "configuration.h"

using namespace std;
using namespace configuration;


LinearEquationSystem::LinearEquationSystem() {
    for (int i=0; i<matrixSize; i++) {
        for (int j=0; j<matrixSize; j++) {
            m_A[i][j] = 0;
        }
        m_b[i] = 0;
    }
}

LinearEquationSystem::LinearEquationSystem(double mat[matrixSize][matrixSize], double vec[matrixSize]) {
    for (int i=0; i<matrixSize; i++) {
        for (int j=0; j<matrixSize; j++) {
            m_A[i][j] = mat[i][j];
        }
        m_b[i] = vec[i];
    }
}



void LinearEquationSystem::gauss(bool withPivotElement) {
    
    cout << "\n\n\nvvv~~~~~~~~ Gauss Method ~~~~~~~~vvv\n\n";
    int index = 0;
    double pivot = 0;
    print_Ax_b();
    
    for (int k=0; k<(matrixSize-1); k++) {
        
        index = k;
        pivot = m_A[k][k];
        //cout << "try pivot: " << pivot << " at index " << index << "\n";
        
        // find (row with) pivot in the column
        // (only if pivot element is activated)
        if (withPivotElement) {
            for (int i=k; i<matrixSize; i++) {
                if (abs(pivot) < abs(m_A[i][k])) {
                    pivot = m_A[i][k];
                    index = i;
                }
            }
            if (k != index) {
                // swap pivot to k-th topline
                swap_rows(k,index);
                index = k;
            }
        }
        
        cout << "(solving the " << (k+2) << ". column with ";
        //cout << "index: " << index << " | ";
        cout << "pivot: " << pivot << ")\n";
        
        // modify rows
        double multiplier = 1;
        for (int i=k; i<matrixSize; i++) {
            if (i != index) {
                // compute new elements of matrix "A"
                if (pivot == 0) {
                    multiplier = 0;
                }
                else {
                    multiplier = m_A[i][k] / pivot;
                }
                for (int j=k; j<matrixSize; j++) {
                    m_A[i][j] = m_A[i][j] - m_A[index][j]*multiplier;
                }
                // compute vector "b"
                m_b[i] = m_b[i] - m_b[index]*multiplier;
            }
        }
        print_Ax_b();
        index = 0;
        pivot = 0;
    }
    solve_gauss();
    print_solution();
    return;
}

void LinearEquationSystem::solve_gauss() {
    // iteration over all rows
    for (int i=matrixSize-1; i>=0; i--) {
        
        // set basis value
        m_sol_A[i] = m_b[i];
        
        // iteration over all columns
        for (int j=matrixSize-1; j>=i; j--) {
            
            // compute with previous elements
            if (j > i) {
                m_sol_A[i] = m_sol_A[i] - m_sol_A[j]*m_A[i][j];
            }
            
            // solve current (diagonal) element
            else {
                if (m_sol_A[i] != 0 && m_A[i][j] != 0) {
                    m_sol_A[i] = m_sol_A[i] / m_A[i][j];
                }
            }
        }
    }
}

void LinearEquationSystem::SOR(double accuracy, double omega) {
    cout << "\n\n\nvvv~~~~~~~~ Seidel Method ~~~~~~~~vvv\n\n";
    print_Ax_b();
    
    // iteration over all rows
    for (int i=0; i<matrixSize; i++) {
        // iteration over all columns
        for (int j=0; j<matrixSize; j++) {
            // compute elements of matrix "B"
            // empty diagonal
            if (i == j) {
                m_B[i][j] = 0;
            }
            // compute other elements
            else {
                if (m_A[i][i] != 0) {
                    m_B[i][j] = (-1) * m_A[i][j] / m_A[i][i];
                }
            }
        }
        // compute vector "c"
        if (m_A[i][i] != 0) {
            m_c[i] = m_b[i] / m_A[i][i];
        }
    }
    print_Bx_c();
    
    // build upper triangular matrix "B_2"
    // iteration over all rows
    for (int i=0; i<matrixSize; i++) {
        // iteration over all columns
        for (int j=0; j<matrixSize; j++) {
            // empty diagonal
            if (i == j) {
                m_B_2[i][j] = 0;
            }
            // upper elements
            else if (i < j) {
                m_B_2[i][j] = m_B[i][j];
            }
            // lower elements
            else {
                m_B_2[i][j] = 0;
            }
        }
    }
    
    if (convergence_condition(accuracy, omega)) {
        solve_SOR(omega);
    }
}

void LinearEquationSystem::solve_SOR(double omega) {
    cout << "iteration --> ";
    // table print
    for (int i=0; i<matrixSize; i++) {
        cout << "x_" << i << " | ";
    }
    cout << "--> ||x^(n) - x^(n-1)||";
    cout << "\n";
    
    // seed vector and old vector setup
    cout << " n = 0  --> ";
    double old_sol_B[matrixSize];
    for (int i=0; i<matrixSize; i++) {
        m_sol_B[i] = 0;
        old_sol_B[i] = 0;
        cout << "0.000" << " | ";
    }
    cout << "\n";
    
    double vec_norm = m_epsilon;
    double abs_value = 0;
    
    // iteration for approximated solution
    int n = 1;
    while (vec_norm >= m_epsilon) {
        vec_norm = 0;
        cout << " n = " << n << "  --> ";
        // iteration over all rows
        for (int i=0; i<matrixSize; i++) {
            m_sol_B[i] = 0;
            // iteration over all columns
            for (int j=0; j<matrixSize; j++) {
                // diagonal elements
                if (i == j) {
                    m_sol_B[i] += (1-omega) * old_sol_B[j];// * m_B[i][j];
                }
                else {
                    // upper diagonal elements
                    if (i < j) {
                        m_sol_B[i] += omega * old_sol_B[j] * m_B[i][j];
                    }
                    // lower diagonal elements
                    else {
                        m_sol_B[i] += omega * m_sol_B[j] * m_B[i][j];
                    }
                }
            }
            m_sol_B[i] += omega * m_c[i];
            cout << m_sol_B[i] << " | ";
        }
        for (int i=0; i<matrixSize; i++) {
            abs_value = abs(m_sol_B[i] - old_sol_B[i]);
            if (abs_value > vec_norm) {
                vec_norm = abs_value;
            }
            old_sol_B[i] = m_sol_B[i];
        }
        cout << "--> " << vec_norm;
        cout << "\n";
        n++;
    }
    cout << "\n\n";
}



bool LinearEquationSystem::convergence_condition(double accuracy, double omega) {
    double B_norm = 0;
    double B_2_norm = 0;
    
    double dummy = 0;
    double dummy_2 = 0;
    
    // iteration over all rows
    for (int i=0; i<matrixSize; i++) {
        // iteration over all columns
        for (int j=0; j<matrixSize; j++) {
            dummy += abs(m_B[i][j]);
            dummy_2 += abs(m_B_2[i][j]);
        }
        if (dummy > B_norm) {
            B_norm = dummy;
        }
        if (dummy_2 > B_2_norm) {
            B_2_norm = dummy_2;
        }
        dummy = 0;
        dummy_2 = 0;
    }
    
    m_epsilon = accuracy * (1 - B_norm) / B_2_norm;
    
    // check for convergence condition
    if (B_norm >= 1) {
        cout << "Error: Convergence condition is not given: ||B|| = " << B_norm << "\n";
        return false;
    }
    else {
        cout << "||B|| = " << B_norm << "\n";
        cout << "||B_2|| = " << B_2_norm << "\n";
        cout << "accuracy = " << accuracy << "\n";
        cout << "epsilon = " << m_epsilon << "\n";
        cout << "omega = " << omega << "\n\n";
        return true;
    }

}

void LinearEquationSystem::swap_rows(int index_a, int index_b) {
    cout << "(swapping rows " << index_a << " and " << index_b << ")\n";
    double temp;
    // left-hand side swap
    for (int i=0; i<matrixSize; i++) {
        temp = m_A[index_a][i];
        m_A[index_a][i] = m_A[index_b][i];
        m_A[index_b][i] = temp;
    }
    // right-hand side swap
    temp = m_b[index_a];
    m_b[index_a] = m_b[index_b];
    m_b[index_b] = temp;
}



void LinearEquationSystem::print_Ax_b() {
    cout << "System: Ax = b (--> A | b)\n";
    for (int i=0; i<matrixSize; i++) {
        for (int j=0; j<matrixSize; j++) {
            cout << m_A[i][j] << " ";
        }
        cout << "  |  " << m_b[i] << "\n";
    }
    cout << "\n";
}

void LinearEquationSystem::print_Bx_c() {
    cout << "System: x = Bx + c (--> B | c)\n";
    for (int i=0; i<matrixSize; i++) {
        cout << "x_" << i << " = ";
        for (int j=0; j<matrixSize; j++) {
            cout << m_B[i][j];
            if (j == (matrixSize-1)) {
                cout << "  |  ";
            }
            else {
                cout << " + "; // << "x_" << j << " + ";
            }
        }
        cout << m_c[i] << "\n";
    }
    cout << "\n";
}

void LinearEquationSystem::print_solution() {
    cout << "Solution:\n";
    for (int i=0; i<matrixSize; i++) {
        cout << "u_" << i << " = " << m_sol_A[i] << "\n";
    }
    cout << "\n";
}
