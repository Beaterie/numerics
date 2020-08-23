//
//  main.cpp
//  advanced-numerics
//
//  Created by Josef Roth on 10.08.20.
//  Copyright © 2020 Josef Roth. All rights reserved.
//

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "configuration.h"
#include "LinearEquationSystem.hpp"
#include "Poisson.hpp"
using namespace std;
using namespace configuration;


int main(int argc, const char * argv[]) {
    
    bool withPivotElement = true;
    double accuracy = 0.01;
    double omega = 1.12;
    
    Poisson u = Poisson();
    u.compute_grid_values();
    u.print_grid();
    u.create_LES();
    u.SOR(accuracy, omega);
    u.gauss(withPivotElement);
    
    
    
    // vvv~~~~~~~~~ example equation systems taken from exercise sheets ~~~~~~~~vvv
    
    // for Gauss (from answer sheet of exercise 5)
    //double matrix[matrixSize][matrixSize] = {{2,-9,5}, {1.2,-5.3999,6}, {1,-1,-7.5}};
    //double vec[matrixSize] = {-4, 0.6001, -8.5};
    
    // for Gauss (exercise sheet 5)
    //double matrix[matrixSize][matrixSize] = {{10,6,2,0}, {5,1,-2,4}, {3,5,1,-1}, {0,6,-2,2}};
    //double vec[matrixSize] = {8,7,2,2};
    
    // for Gauss (from answer sheet of exercise 5)
    //double matrix[matrixSize][matrixSize] = {{2,1,-5,1}, {1,-3,0,-6}, {0,2,-1,2}, {1,4,-7,6}};
    //double vec[matrixSize] = {8,9,-5,0};
    
    // for SOR (from exercise sheet 8)
    //double matrix[matrixSize][matrixSize] = {{6.25,-1,0.5}, {-1,5,2.12}, {0.5,2.12,3.6}};
    //double vec[matrixSize] = {7.5,-8.68,-0.24};

    //LinearEquationSystem les = LinearEquationSystem(matrix, vec);
    //les.SOR(accuracy, omega);
    //les.gauss(true);
    
    return 0;
}


