//--------------------------------------------------------------------------
// AdaptiveGradientDescent: An adaptive step size gradient descent
// optimizer
//--------------------------------------------------------------------------

#ifndef ADAPTIVEGRADIENTDESCENT_H
#define ADAPTIVEGRADIENTDESCENT_H

#include "optimizer.h"

class AdaptiveGradientDescent : public Optimizer
{

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructor
        //--------------------------------------------------------------------------

        AdaptiveGradientDescent();
        AdaptiveGradientDescent(Algorithm* algorithm);
        ~AdaptiveGradientDescent();

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        void SetAlgorithm(Algorithm* algorithm);

        virtual void SetMaxIterations(int maxIters);
        virtual void SetBreakRatio(double breakRatio);
        virtual void SetStepsize(double stepsize);

        //--------------------------------------------------------------------------
        // Implementation of Optimzer interface
        //--------------------------------------------------------------------------

        virtual void Optimize(Array3D<double>& X);

};

#endif // ADAPTIVEGRADIENTDESCENT_H
