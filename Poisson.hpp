//
//  Poisson.hpp
//  advanced-numerics
//
//  Created by Josef Roth on 17.08.20.
//  Copyright © 2020 Josef Roth. All rights reserved.
//

#ifndef Poisson_hpp
#define Poisson_hpp

#include <stdio.h>
#include "configuration.h"
#include "LinearEquationSystem.hpp"

using namespace std;
using namespace configuration;

class Poisson {
    
public:
    // Members
    double boundaries[steps][steps];
    LinearEquationSystem les = LinearEquationSystem();
    
    
    // Constructor
    Poisson();
    
    
    // Functions
    void compute_grid_values();
    double function_u(int i, int j);
    void create_LES();
    
    void print_grid();
    void gauss(bool withPivotElement);
    void SOR(double accuracy, double omega);
    
};

#endif /* Poisson_hpp */
