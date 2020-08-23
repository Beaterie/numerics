//
//  configuration.h
//  advanced-numerics
//
//  Created by Josef Roth on 17.08.20.
//  Copyright © 2020 Josef Roth. All rights reserved.
//

#ifndef configuration_h
#define configuration_h

using namespace std;


namespace configuration
{
    const int steps = 5; // N
    const double lowerBoundary = 0; // a
    const double upperBoundary = 1; // b
    const double stepSize = (upperBoundary - lowerBoundary) / (steps - 1); // h
    const int matrixSize = (steps-2)*(steps-2); // number of unknowns OR (when testing
                                                // a simple les without Poisson)
                                                // matrix/vector dimension
}

#endif
