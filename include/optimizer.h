//--------------------------------------------------------------------------
// Optimizer: Virutal base class for optimizers
//--------------------------------------------------------------------------

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "algorithm.h"
#include "array3d.h"

class Optimizer
{

    protected:

        int _maxIters;
        double _breakRatio;
        double _stepsize;

        // Needed for access to functions to compute functional and gradient
        Algorithm* _algorithm;

    public:

        virtual void SetMaxIterations(int maxIters) = 0;
        virtual void SetBreakRatio(double breakRatio) = 0;
        virtual void SetStepsize(double stepsize) = 0;

        // Interface for user of class to start the optimization
        virtual void Optimize(Array3D<double>& X) = 0;
        virtual ~Optimizer() {};
};

#endif // OPTIMIZER_H
