//
//  Poisson.cpp
//  advanced-numerics
//
//  Created by Josef Roth on 17.08.20.
//  Copyright © 2020 Josef Roth. All rights reserved.
//

#include "Poisson.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "configuration.h"
#include "LinearEquationSystem.hpp"

using namespace std;
using namespace configuration;

Poisson::Poisson() {
    int index = 0;
    for (int i=0; i<steps; i++) {
        for (int j=0; j<steps; j++) {
            // random values
            boundaries[i][j] = index;
            index++;
        }
    }
}

void Poisson::compute_grid_values() {
    for (int i=0; i<steps; i++) {
        boundaries[i][0] = function_u(i, 0);
        boundaries[i][steps-1] = function_u(i, steps-1);
    }
    for (int j=0; j<steps; j++) {
        boundaries[0][j] = function_u(0, j);
        boundaries[steps-1][j] = function_u(steps-1, j);
    }
}

double Poisson::function_u(int i, int j) {
    double result = (i*stepSize - lowerBoundary) * (j*stepSize - lowerBoundary) *
                    (i*stepSize - upperBoundary) * (j*stepSize - upperBoundary);
    return result;
}

void Poisson::create_LES() {
    int index_range = steps-2;
    int index = 0;
    // iteration over inner stencils
    for (int i=1; i<(steps-1); i++) {
        for (int j=1; j<(steps-1); j++) {
            
            // set right-hand side
            les.m_b[index] = pow(stepSize,2) *
                                    ( function_u(i-1, j) + function_u(i+1, j) - 4*function_u(i, j)
                                    + function_u(i, j-1) + function_u(i, j+1) );
            // check for border stencils
            if (i == 1) {
                les.m_b[index] = les.m_b[index] - function_u(i-1, j);
            }
            else if (i == (steps-2)) {
                les.m_b[index] = les.m_b[index] - function_u(i+1, j);
            }
            if (j == 1) {
                les.m_b[index] = les.m_b[index] - function_u(i, j-1);
            }
            else if (j == (steps-2)) {
                les.m_b[index] = les.m_b[index] - function_u(i, j+1);
            }
            
            // set left-hand side
            les.m_A[index][index] = -4;
            // check upper stencil
            if ((i-1) >= 1) {
                les.m_A[index][(index-index_range)] = 1;
            }
            // check lower stencil
            if ((i+1) < (steps-1)) {
                les.m_A[index][(index+index_range)] = 1;
            }
            // check left stencil
            if ((j-1) >= 1) {
                les.m_A[index][(index-1)] = 1;
            }
            // check right stencil
            if ((j+1) < (steps-1)) {
                les.m_A[index][(index+1)] = 1;
            }
            
            index++;
        }
    }
    les.print_Ax_b();
}

void Poisson::print_grid() {
    cout << "Grid:\n";
    int index = 0;
    for (int i=0; i<steps; i++) {
        for (int j=0; j<steps; j++) {
            if (boundaries[i][j] != 0) {
                cout << " u_" << index << " ";
                index++;
            }
            else {
                cout << " 0 ";
            }
        }
        cout << "\n";
    }
    cout << "\n";
}

void Poisson::gauss(bool withPivotElement) {
    les.gauss(withPivotElement);
}

void Poisson::SOR(double accuracy, double omega) {
    les.SOR(accuracy, omega);
}
