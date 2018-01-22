//--------------------------------------------------------------------------
// Algorithm:  Virtual base class for regression algorithms.  Defines
// the interface that must be implemented for use with an Optimizer
//--------------------------------------------------------------------------

#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "array3d.h"

class Algorithm
{

    private:

        virtual void WriteSelf() = 0;
    public:

        //--------------------------------------------------------------------------
        // Interface that inherited classes must implement
        //--------------------------------------------------------------------------

        virtual double ComputeFunctional(const Array3D<double>& X, double &dataMatching, double &regularity) = 0;
        virtual Array3D<double> ComputeGradient(const Array3D<double>& X) = 0;
        virtual ~Algorithm() {};

};

#endif // ALGORITHM_H
