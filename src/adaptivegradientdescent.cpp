#include "adaptivegradientdescent.h"
#include "shape4dstate.h"

#include <stdio.h>			// For printf

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXAdaptiveGradientDescent
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
AdaptiveGradientDescent::AdaptiveGradientDescent()
{
    this->_algorithm = 0;
    this->_stepsize = 0.5;
    this->_maxIters = 500;
    this->_breakRatio = 1e-6;
}

//----------------------------------------------------------------
// EXAdaptiveGradientDescent
//----------------------------------------------------------------
// Inputs:
//   algorithm - pointer to an algorithm object which is needed
//               so the optimizer can call appropriate methods
//               to compute the functional and gradient
//
// Outputs:
//----------------------------------------------------------------
// Init constructor - Recommended constructor
//----------------------------------------------------------------
AdaptiveGradientDescent::AdaptiveGradientDescent(Algorithm* algorithm)
{
    this->_algorithm = algorithm;
    this->_stepsize = 0.5;
    this->_maxIters = 500;
    this->_breakRatio = 1e-6;
}

//----------------------------------------------------------------
// ~EXAdaptiveGradientDescent
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
AdaptiveGradientDescent::~AdaptiveGradientDescent()
{
    delete this->_algorithm;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetMaxIterations
//----------------------------------------------------------------
// Inputs:
//   maxIters - the maximum number of iterations
//
// Outputs:
//----------------------------------------------------------------
// Setter for maxIters
//----------------------------------------------------------------
void AdaptiveGradientDescent::SetMaxIterations(int maxIters)
{
    this->_maxIters = maxIters;
}

//----------------------------------------------------------------
// SetBreakRatio
//----------------------------------------------------------------
// Inputs:
//   maxIters - the break ratio for convergence
//
// Outputs:
//----------------------------------------------------------------
// Setter for breakRatio
//----------------------------------------------------------------
void AdaptiveGradientDescent::SetBreakRatio(double breakRatio)
{
    this->_breakRatio = breakRatio;
}

//----------------------------------------------------------------
// SetStepsize
//----------------------------------------------------------------
// Inputs:
//   stepsize - the starting stepsize for gradient descent
//
// Outputs:
//----------------------------------------------------------------
// Setter for stepsize
//----------------------------------------------------------------
void AdaptiveGradientDescent::SetStepsize(double stepsize)
{
    this->_stepsize = stepsize;
}

//----------------------------------------------------------------
// SetAlgorithm
//----------------------------------------------------------------
// Inputs:
//   algorithm - pointer to an algorithm object which is needed
//               so the optimizer can call appropriate methods
//               to compute the functional and gradient
//
// Outputs:
//----------------------------------------------------------------
// Setter for the algorithm pointer
//----------------------------------------------------------------
void AdaptiveGradientDescent::SetAlgorithm(Algorithm* algorithm)
{
    this->_algorithm = algorithm;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF EXOptimzer INTERFACE
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// Optimize
//----------------------------------------------------------------
// Inputs:
//   X - 3D array of data which is changed during optimzation
//       and returned by reference
//
// Outputs:
//----------------------------------------------------------------
// Optimize using adaptive step size gradient descent
//----------------------------------------------------------------
void AdaptiveGradientDescent::Optimize(Array3D<double>& X)
{
    Shape4DState::SetHasConverged(false);

    int dim = X.GetLength();
    int nx = X.GetWidth();
    int T = X.GetHeight();

    //printf("OPTIMIZER: Done getting dim, nx, and T...\n");

    // Now we start adaptive gradient descent
    int maxIter = this->_maxIters;
    double* dataMatchingValues = new double[maxIter];
    double* regularityValues = new double[maxIter];
    double* criterionValues = new double[maxIter];

    //printf("OPTIMIZER: Computing functional at 0...\n");

    // Compute the functional by following the algorithm's pointer
    criterionValues[0] = this->_algorithm->ComputeFunctional(X, dataMatchingValues[0], regularityValues[0]);

    //printf("OPTIMIZER: Done computing functional at 0...\n");

    // Constants to increase/decrease step size
    double stepIncrease = 2.0f;
    double stepDecrease = 0.25f;
    // Optimization termination criteria
    double breakratio = this->_breakRatio;
    // The max number of tries to find a 'good' time step per iteration
    int loopbreak = 40;

    _stepsize = 0.0;

    // Initialize the step size if it is zero
    if (_stepsize == 0.0)
    {
        // Compute the gradient by following the algorithm's pointer
        Array3D<double> G = this->_algorithm->ComputeGradient(X);

        //printf("OPTIMIZER: Done computing gradient...\n");

        // The initial step size is the functional value over the sum of squared gradient
        double sumSquared = 0.0f;
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                for (int k=0; k<T; k++)
                {
                    sumSquared += G(i,j,k)*G(i,j,k);
                }
            }
        }
        _stepsize = (criterionValues[0]/sumSquared)*0.01;
    }
    // Else use the current stepsize * stepIncrease
    else
    {
        _stepsize = _stepsize * stepIncrease;

        if (_stepsize < 1e-8)
            _stepsize = 0.0001;
    }

    printf("\n");

    // Temporary values of X
    Array3D<double> Xnew(dim,nx,T);

    for (int i=1; i<maxIter; i++)
    {
        Array3D<double> G = this->_algorithm->ComputeGradient(X);

        // Boolean which determines if we keep searching for a 'good' time step
        bool minimtest = 0;
        // Inner loop counter
        int loop = 0;

        // We keep decreasing the step size until one is found that lowers
        // the functional between iterations.  Or we stop after a maximium
        // number of iterations is reached.
        while (!minimtest && (loop<loopbreak))
        {
            // Lower the step size
            if (loop != 0)
            {
                _stepsize = _stepsize * stepDecrease;
            }

            for (int j=0; j<dim; j++)
            {
                for (int k=0; k<nx; k++)
                {
                    for (int l=0; l<T; l++)
                    {
                        Xnew(j,k,l) = X(j,k,l) - _stepsize*G(j,k,l);
                    }
                }
            }

            // Compute the functional values at the new positions
            criterionValues[i] = this->_algorithm->ComputeFunctional(Xnew, dataMatchingValues[i], regularityValues[i]);

            // Did this stepsize lead to a reduction in the functional?
            minimtest = (criterionValues[i] < criterionValues[i-1]);
            // Increase the loop counter
            loop = loop + 1;

        }

        // Update the positions of source points
        X = Xnew;

        // Update the iteration values for saving
        Shape4DState::UpdateIteration(Shape4DState::GetIteration()+1, dataMatchingValues[i-1], regularityValues[i-1], _stepsize);
        printf("Iteration %3d   funct = %0.4f   data = %0.4f   reg = %0.4f   step = %0.10lf\n",
               Shape4DState::GetIteration(), criterionValues[i-1], dataMatchingValues[i-1], regularityValues[i-1], _stepsize);

        // If it is time to save intermediate shapes
        if ((i%Shape4DState::GetSaveProgressEveryNIterations() == 0) && (Shape4DState::GetShouldSave()))
        {
            Shape4DState::SetShouldWriteShapesAndVectors(true);
        }

        // Update the step size (increase the step size from previous iteration)
        _stepsize = _stepsize * stepIncrease;

        // Termination criteria
        if ((criterionValues[i-1]-criterionValues[i]) <
                (breakratio*(criterionValues[0]-criterionValues[i])) || (loop==loopbreak))
        {
            break;
        }
    }

    delete[] dataMatchingValues;
    delete[] regularityValues;
    delete[] criterionValues;

    Shape4DState::SetHasConverged(true);
    return;

}
