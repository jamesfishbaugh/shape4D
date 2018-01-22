//--------------------------------------------------------------------------
// EXRegressionAcceleration:  Shape regression parameterized by
// acceleration.
//--------------------------------------------------------------------------

#include "regressionacceleration.h"

#include <stdio.h>						// For printf
#include <math.h>						// For pow
#include <stdlib.h>						// For free

#include "optimizer.h"                  // The optimizer
#include "adaptivegradientdescent.h"	// Specific implementation of optimizer
#include "shape4dstate.h"               // To save intermediate results during optimization

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXRegressionAcceleration
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
RegressionAcceleration::RegressionAcceleration()
{
    this->_iteration = 0;
    this->SetNumTargets(0);
    this->_targets = 0;
    this->_continueRegression = false;
}

//----------------------------------------------------------------
// EXRegressionAcceleration
//----------------------------------------------------------------
// Inputs: 
//   source - reference to a source object
//
// Outputs:
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
RegressionAcceleration::RegressionAcceleration(const RegressionParams& source)
{
    this->_iteration = 0;
    this->SetSource(source);
    this->SetNumTargets(0);
    this->_targets = 0;
    this->_continueRegression = false;
}

//----------------------------------------------------------------
// EXRegressionAcceleration
//----------------------------------------------------------------
// Inputs: 
//   source - reference to a source object
//   numTargets - number of target objects
//   targets - array of pointers to target objects
//
// Outputs:
//----------------------------------------------------------------
// Init constructor - recommended constructor
//----------------------------------------------------------------
RegressionAcceleration::RegressionAcceleration(const RegressionParams& source, int numTargets, TargetData** targets)
{
    this->_iteration = 0;
    this->SetSource(source);
    this->SetNumTargets(numTargets);
    this->SetTargets(targets);
    this->_continueRegression = false;
}

//----------------------------------------------------------------
// ~EXRegressionAcceleration
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
RegressionAcceleration::~RegressionAcceleration()
{
    delete [] this->_targets;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

void RegressionAcceleration::SetImpulse(Array3D<double> &impulse)
{
    this->_impulse = impulse;
    this->_continueRegression = true;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

Array3D<double> RegressionAcceleration::GetX()
{
    return this->_X;
}

//----------------------------------------------------------------
// GetVelocity
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 3D array (dim, num_pts, time) of velocity
//----------------------------------------------------------------
// Returns velocity vectors over time
//----------------------------------------------------------------
Array3D<double> RegressionAcceleration::GetVelocity()
{
    return this->_dX;
}

//----------------------------------------------------------------
// GetAcceleration
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 3D array (dim, num_pts, time) of acceleration
//----------------------------------------------------------------
// Returns acceleration vectors over time
//----------------------------------------------------------------
Array3D<double> RegressionAcceleration::GetAcceleration()
{
    return this->_accel;
}

//----------------------------------------------------------------
// GetImpulse
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 3D array (dim, num_pts, time) of impulse
//----------------------------------------------------------------
// Returns impulse vectors over time
//----------------------------------------------------------------
Array3D<double> RegressionAcceleration::GetImpulse()
{
    return this->_impulse;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// MAIN ALGORITHM METHOD
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// Run
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - A 3D array (dim, num_pts, time) of point 
//            trajectories
//----------------------------------------------------------------
// Starts the regression algorithm
//----------------------------------------------------------------
Array3D<double> RegressionAcceleration::Run()
{
    this->_tau = (double)1.0f/(double)(this->_source.GetT()-1.0f);

    // Initialize member arrays
    this->InitX();
    this->InitInitV0();
    this->InitVelocity();

    if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
    {
        // Create the grids for optimization
        this->SetSourceGrids();
        this->SetTargetGrids();
    }

    //printf("Done setting grids\n");

    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();
    Array2D<double> X0 = this->_source.GetX();

    Array2D<double> sobolevGrad(dim,nx);
    sobolevGrad.FillArray(0.0);
    this->_sobolevGrad = sobolevGrad;


    // Allocate memory for graddesc (stores the positions and velocities adding 2 dimensions to T)
    Array3D<double> graddesc(dim, nx, T+2);

    // The graddesc variable has x0 and v0 concatenated at the end of the impulses
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            for (int k=0; k<T+2; k++)
            {
                if (k==T)
                {
                    graddesc(i,j,k) = this->_initV0(i,j);
                }
                else if (k==(T+1))
                {
                    graddesc(i,j,k) = X0(i,j);
                }
                else
                {
                    if (this->_continueRegression)
                    {
                        graddesc(i,j,k) =  this->_impulse(i,j,k);
                    }
                    else
                    {
                        graddesc(i,j,k) =  0.0f;
                    }
                }
            }
        }
    }

    // Run the optimization
    Array3D<double> impulse(dim, nx, T);
    if (this->_continueRegression)
    {
        impulse = this->_impulse;
    }
    else
    {
        impulse.FillArray(0.0);
    }

    Array2D<double> V0(dim, nx);
    V0 = this->_initV0;

    if (this->_source.ShouldUseFista())
    {
        this->FISTA(impulse, X0, V0);
    }
    else
    {
        this->GradientDescent(impulse, X0, V0);
//        Optimizer* optimizer = (Optimizer*) new AdaptiveGradientDescent(this);
//        optimizer->SetMaxIterations(this->_source.GetMaxIters());
//        optimizer->SetBreakRatio(this->_source.GetBreakRatio());
//        optimizer->SetStepsize(0.05f);
//        optimizer->Optimize(graddesc);

//        this->SplitImpulseAndV0(graddesc, impulse, V0);

//        // Clean up memory
//        delete optimizer;
    }

    // Compute the final trajectories now that we are done
    this->ComputeTrajectories(impulse, this->_X.Get2DSliceAtHeight(0), V0);

    // Compute the final acceleration
    Array3D<double> accel(dim, nx, T);
    accel.FillArray(0.0f);
    for (int t=0; t<T-1; t++)
    {
        // Compute acceleration at time t
        Array2D<double> accelt;

        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            accelt = this->GridKernelSum(t, impulse);
        }
        else
        {
            accelt = this->RegKernelSum(t, impulse);
        }

        accel.Set2DSliceAtHeight(accelt,t);
    }

    // Save the final values
    this->_accel = accel;
    this->_impulse = impulse;

    // Save the shapes and vectors
    vector<char *> names;
    names.push_back((char*)"velocity");
    names.push_back((char*)"acceleration");
    names.push_back((char*)"impulse");

    vector< Array3D<double> > vectors;
    vectors.push_back(this->_dX);
    vectors.push_back(this->_accel);
    vectors.push_back(this->_impulse);

    Shape4DState::SetX(_X);
    Shape4DState::SetVectors(names, vectors);
    Shape4DState::SaveShapesAndVectors();
    Shape4DState::SaveStateAccel();

    Shape4DState::SetShouldWriteShapesAndVectors(false);

    // Return the final trajectores
    return this->_X;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// INITIALIZATION METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// InitInitV0
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Initializes the 2D array (dim, num_pts) this->_initV0 which 
// stores initial velocity
//----------------------------------------------------------------	
void RegressionAcceleration::InitInitV0()
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();

    this->_initV0 = this->_source.GetInitV0();

    // If there is no initV0, return zeros
    if ((this->_initV0.GetLength() == 0) || (this->_initV0.GetWidth() == 0))
    {
        Array2D<double> tempInitV0(dim,nx);
        tempInitV0.FillArray(0.0f);

        this->_initV0 = tempInitV0;
    }
}

//----------------------------------------------------------------
// InitVelocity
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Initializes the 3D array (dim, num_pts, time) this->_dX which 
// stores point velocity
//----------------------------------------------------------------	
void RegressionAcceleration::InitVelocity()
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    Array3D<double> tempVelocity(dim,nx,T);
    tempVelocity.FillArray(0.0f);
    this->_dX = tempVelocity;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// COMPUTATION OF REGRESSION FUNCTIONAL
//----------------------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------
// ComputeRegularity
//----------------------------------------------------------------
// Inputs:
//   impulse - 3D array (dim, num_pts, time) of impulse vectors
//
// Outputs:
//   return - regularity of deformation
//----------------------------------------------------------------
// Computes the regularity of the deformation, interpreted as 
// the total amount of acceleration
//----------------------------------------------------------------	
double RegressionAcceleration::ComputeRegularity(const Array3D<double>& impulse)
{
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();
    double E = 0.0f;

    for (int t=0; t<T; t++)
    {
        // Compute acceleration at time t
        Array2D<double> accelt;

        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            accelt = this->GridKernelSum(t, impulse);
        }
        else
        {
            accelt = this->RegKernelSum(t, impulse);
        }

        // Regularity is the dot product between impulse and acceleration
        //#pragma omp parallel for reduction(+:E)
        for (int i=0; i<nx; i++)
        {
            double dot = (impulse(0,i,t)*accelt(0,i)) +
                    (impulse(1,i,t)*accelt(1,i)) +
                    (impulse(2,i,t)*accelt(2,i));

            E += dot;
        }
    }

    return E;
}


//----------------------------------------------------------------
// ComputeFunctional
//----------------------------------------------------------------
// Inputs:
//   impulse -
//   X0 -
//   V0 -
//   dataMatching - double for data matching value return
//   regularity - double for regularity value return
//
// Outputs:
//   dataMatching - data matching value
//   regularity - regularity value
//   return - total value of the regression functional
//----------------------------------------------------------------
// Computes the value of the regression criteria by computing
// the total data matching and regularity terms.  Returns
// data matching and regularity independently through reference
// parameters as well as returning the functional value.
//----------------------------------------------------------------
double RegressionAcceleration::ComputeFunctional(const Array3D<double>& impulse, const Array2D<double> X0, const Array2D<double> V0, double &dataMatching, double &regularity)
{
    // Apply the current deformation
    this->ComputeTrajectories(impulse, X0, V0);

    //printf("    Computing data matching...\n");
    dataMatching = ComputeDataMatching(this->_X);
    //printf("    Computing regularity...\n");
    regularity = ComputeRegularity(impulse);

    //printf("Data = %0.4f  Reg = %0.4f\n", dataMatching, regularity);

    //printf("Done computing functional\n");

    return dataMatching + this->_source.GetGamma()*regularity;
}

//----------------------------------------------------------------
// ComputeFunctional
//----------------------------------------------------------------
// Inputs:
//   impulseAndV0 - 3D array (dim, num_pts, time+1) of impulse
//                  vectors (:,:,0:time) and initial velocity 
//                  (:,:,time+1)
//   dataMatching - double for data matching value return
//   regularity - double for regularity value return
//
// Outputs:
//   dataMatching - data matching value
//   regularity - regularity value
//   return - total value of the regression functional
//----------------------------------------------------------------
// Computes the value of the regression criteria by computing
// the total data matching and regularity terms.  Returns
// data matching and regularity independently through reference
// parameters as well as returning the functional value.
//----------------------------------------------------------------	
double RegressionAcceleration::ComputeFunctional(const Array3D<double>& impulseAndV0, double &dataMatching, double &regularity)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    // Split out the impulse and initial velocity
    Array3D<double> impulse(dim,nx,T);
    Array2D<double> v0(dim,nx);
    this->SplitImpulseAndV0(impulseAndV0, impulse, v0);

    //printf("  Computing functional value\n");
    //printf("    Computing trajectories...\n");

    // Apply the current deformation
    this->ComputeTrajectories(impulse, v0);

    //printf("    Computing data matching...\n");
    dataMatching = ComputeDataMatching(this->_X);
    //printf("    Computing regularity...\n");
    regularity = ComputeRegularity(impulse);

    //printf("Data = %0.4f  Reg = %0.4f\n", dataMatching, regularity);

    //printf("Done computing functional\n");

    return dataMatching + this->_source.GetGamma()*regularity;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// COMPUTATION OF GRADIENT OF REGRESSION FUNCTIONAL
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeGradient
//----------------------------------------------------------------
// Inputs:
//   impulseAndV0 - 3D array (dim, num_pts, time+1) of impulse
//                  vectors (:,:,0:time) and initial velocity 
//                  (:,:,time+1)
//
// Outputs:
//   return - 3D array (dim, num_pts, time) of the gradient of
//            the regression functional
//----------------------------------------------------------------
// Computes the gradient of regression functional.  Used to 
// update impulse and initial velocity.
//----------------------------------------------------------------
Array3D<double> RegressionAcceleration::ComputeGradient(const Array3D<double>& impulseAndV0)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    // Split out impulse and initial velocity
    Array3D<double> impulse(dim,nx,T);
    Array2D<double> V0 (dim,nx);
    this->SplitImpulseAndV0(impulseAndV0, impulse, V0);
    Array2D<double> X0 (dim,nx);
    X0 = this->_X.Get2DSliceAtHeight(0);

    Array3D<double> gradImpulse(dim,nx,T);
    gradImpulse.FillArray(0.0);
    Array2D<double> gradX0 (dim,nx);
    gradX0.FillArray(0.0);
    Array2D<double> gradV0 (dim,nx);
    gradV0.FillArray(0.0);

    this->ComputeGradient(impulse, X0, V0, gradImpulse, gradX0, gradV0);

    // Package of the final gradient calculation
    Array3D<double> G(dim,nx,T+2);
    G.FillArray(0.0f);

    // Compute the final value of the gradient
    //#pragma omp parallel for collapse(3)
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            for (int k=0; k<T+2; k++)
            {
                if (k==T)
                    // Initial velocity
                {
                    G(i,j,k) = gradV0(i,j);

                }
                else if (k==(T+1))
                    // Initial position
                {
                    G(i,j,k) = gradX0(i,j);
                }
                // Impulse
                else
                {
                    G(i,j,k) = gradImpulse(i,j,k);
                }
            }
        }
    }

    return G;
}

//----------------------------------------------------------------
// ComputeGradient
//----------------------------------------------------------------
// Inputs:
//   impulse -
//   X0 -
//   V0 -
//   gradImpulse -
//   gradX0 -
//   gradV0 -
//
// Outputs:
//   return - 3D array (dim, num_pts, time) of the gradient of
//            the regression functional
//----------------------------------------------------------------
// Computes the gradient of regression functional.  Used to
// update impulse and initial velocity.
//----------------------------------------------------------------
void RegressionAcceleration::ComputeGradient(const Array3D<double>& impulse, const Array2D<double>&, const Array2D<double>& V0,
                                             Array3D<double>& gradImpulse, Array2D<double>& gradX0, Array2D<double>& gradV0)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();
    double sigmaV2 = pow(this->_source.GetSigmaV(),2);

    //printf("  Computing gradient of functional\n");
    //printf("    Computing trajectories...\n");

    // Apply the deformation
    this->ComputeTrajectories(impulse, V0);

    // Save shape and vector data when needed
    if (Shape4DState::GetShouldWriteShapesAndVectors())
    {
        // Compute the final acceleration
        Array3D<double> accel(dim, nx, T);
        accel.FillArray(0.0f);
        for (int t=0; t<T-1; t++)
        {
            // Compute acceleration at time t
            Array2D<double> accelt;

            if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
            {
                accelt = this->GridKernelSum(t, impulse);
            }
            else
            {
                accelt = this->RegKernelSum(t, impulse);
            }

            accel.Set2DSliceAtHeight(accelt,t);
        }

        // Save the final values
        this->_accel = accel;
        this->_impulse = impulse;

        // Save the shapes and vectors
        vector<char *> names;
        names.push_back((char*)"velocity");
        names.push_back((char*)"acceleration");
        names.push_back((char*)"impulse");

        vector< Array3D<double> > vectors;
        vectors.push_back(this->_dX);
        vectors.push_back(this->_accel);
        vectors.push_back(this->_impulse);

        Shape4DState::SetX(_X);
        Shape4DState::SetVectors(names, vectors);
        Shape4DState::SaveShapesAndVectors();
        Shape4DState::SaveStateAccel();

        Shape4DState::SetShouldWriteShapesAndVectors(false);
    }

    //printf("    Computing gradient...\n");

    // The sum of etaAux represents the contribution of the sum of grad A
    Array3D<double> etaAux = this->ComputeDataMatchingGradient(this->_X);
    double sumetaAux = etaAux.Sum();

    // Create etax and etaxdot
    Array3D<double> etax(dim,nx,T);
    etax.FillArray(0.0f);
    Array3D<double> etaxdot(dim,nx,T);
    etaxdot.FillArray(0.0f);

    // Create sumdeta
    Array2D<double> sumdeta(dim, nx);
    sumdeta.FillArray(0.0f);

    // Create detat and detatm1
    Array2D<double> detat(dim,nx);
    Array2D<double> auxt(dim, nx);
    Array2D<double> detatm1(dim,nx);
    Array2D<double> auxtm1(dim, nx);

    // Computes etax[t] and etaxdot[t] with backward integration centered scheme
    for (int t=T-1; t>=1; t--)
    {
        // If the sum of etaAux is non-zero
        if (sumetaAux != 0)
        {
            // Update the sum of deta
            //#pragma omp parallel for collapse(2)
            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<nx; j++)
                {
                    sumdeta(i,j) += etaAux(i,j,t);
                    // Current value of etax is equal to this sum
                    etax(i,j,t) = sumdeta(i,j);
                }
            }
        }

        // Compute derivaitve of eta at time t
        detat.FillArray(0.0f);
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            detat = this->GridCompDeta(t, etaxdot, impulse);
        }
        else
        {
            detat = this->RegCompDeta(t, etaxdot, impulse);
        }
        auxt.FillArray(0.0f);

        // A small step in the direction of deta
        //#pragma omp parallel for collapse(2)
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                auxt(i,j) = (this->_tau*2/sigmaV2) * detat(i,j);
                etax(i,j,t-1) = etax(i,j,t) + auxt(i,j);
                sumdeta(i,j) = sumdeta(i,j) + auxt(i,j);

                etaxdot(i,j,t-1) = etaxdot(i,j,t) + etax(i,j,t-1)*this->_tau;

            }
        }

        // Compute derivaitve of eta at time t-1
        detatm1.FillArray(0.0f);
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            detatm1 = this->GridCompDeta(t-1, etaxdot, impulse);
        }
        else
        {
            detatm1 = this->RegCompDeta(t-1, etaxdot, impulse);
        }
        auxtm1.FillArray(0.0f);

        // Sum the two detas
        detatm1 = detatm1 + detat;

        //  A small step in the direction of deta
        //#pragma omp parallel for
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                auxtm1(i,j) = (this->_tau/sigmaV2) * detatm1(i,j);
                etax(i,j,t-1) = etax(i,j,t) + auxtm1(i,j);
                sumdeta(i,j) = sumdeta(i,j) + auxtm1(i,j);

                etaxdot(i,j,t-1) = etaxdot(i,j,t) + etax(i,j,t-1)*this->_tau;
            }
        }
    }

    // Storage of the final gradient calculation
    Array3D<double> G(dim,nx,T+2);
    G.FillArray(0.0f);

    double gammaTimesTwo = 2.0f*this->_source.GetGamma();

    // Compute the final value of the gradient
    //#pragma omp parallel for collapse(3)
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            for (int k=0; k<T; k++)
            {
                if (k==0)
                {
                    if (this->_source.ShouldEstimateBaseline())
                    {
                        gradX0(i,j) = etax(i,j,0);
                    }
                    else
                    {
                        gradX0(i,j) = 0.0;
                    }

                    gradV0(i,j) = this->_source.GetV0Weight() * etaxdot(i,j,0);
                }

                gradImpulse(i,j,k) = gammaTimesTwo*impulse(i,j,k) + etaxdot(i,j,k);
            }
        }
    }

    if (this->_source.ShouldEstimateBaseline())
    {
        // We need to ensure the baseline shape doesn't self intersect by convolving the gradient with a kernel
        this->_sobolevGrad.FillArray(0.0);
        double stdV = this->_source.GetStdV();
        double smoothingFactor = this->_source.GetBaselineSmoothing();

        //#pragma omp parallel for collapse(2)
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nx; j++)
            {
                double argin = - ( pow(this->_X(0,i,0)-this->_X(0,j,0),2) +
                                   pow(this->_X(1,i,0)-this->_X(1,j,0),2) +
                                   pow(this->_X(2,i,0)-this->_X(2,j,0),2) ) / (smoothingFactor*sigmaV2);

                double argout = stdV*exp(argin);
                //printf("%0.4f\n", argout);

                this->_sobolevGrad(0,i) = this->_sobolevGrad(0,i) + argout*gradX0(0,j);
                this->_sobolevGrad(1,i) = this->_sobolevGrad(1,i) + argout*gradX0(1,j);
                this->_sobolevGrad(2,i) = this->_sobolevGrad(2,i) + argout*gradX0(2,j);
            }
        }

        if (this->_source.ShouldUseFista() == false)
        {
            gradX0 = this->_sobolevGrad;
        }
    }

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SplitImpulseAndV0
//----------------------------------------------------------------
// Inputs:
//   impulseAndV0 - 3D array (dim, num_pts, time+1) of impulse
//                  vectors (:,:,0:time) and initial velocity 
//                  (:,:,time+1)
//   impulse - 3D array (dim, num_pts, time) of impulse as
//             storage for return by reference
//   v0 - 2D array (dim, num_pts) of initial velocity as 
//        storage for return by reference
//
// Outputs:
//   impulse - 3D array (dim, num_pts, time) of impulse
//   v0 - 2D array (dim, num_pts) of initial velocity
//----------------------------------------------------------------
// Splits the impulse and initial velocity from the 3D array
// which is the concatenation of the two
//----------------------------------------------------------------
void RegressionAcceleration::SplitImpulseAndV0(const Array3D<double>& impulseAndV0, Array3D<double>& impulse, Array2D<double>& v0)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    // Split impulse and v0
    //#pragma omp parallel for collapse(3)
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            for (int k=0; k<T+2; k++)
            {
                // Initial velocity
                if (k==T)
                {
                    v0(i,j) = impulseAndV0(i,j,k);
                }
                else if (k==(T+1))
                {
                    this->_X(i,j,0) = impulseAndV0(i,j,k);
                }
                // Impulse
                else
                {
                    impulse(i,j,k) = impulseAndV0(i,j,k);
                }
            }
        }
    }
}

//----------------------------------------------------------------
// SplitGradient
//----------------------------------------------------------------
// Inputs:
//   G - 3D array (dim, num_pts, time+2) of impulse vectors
//       (:,:,0:time) and V0(:,:,time+1) and X0(:,:,time+2)
//   gradImpulse - 3D array (dim, num_pts, time) of grad
//                  impulse as storage for return by reference
//   gradX0 - 2D array (dim, num_pts) of grad X0 as storage for
//            or return by reference
//   gradV0 - 2D array (dim, num_pts) of grad V0 as storage for
//            or return by reference
//
// Outputs:
//   gradImpulse - 3D array (dim, num_pts, time) of grad impulse
//   gradX0 - 2D array (dim, num_pts) of grad X0
//   gradV0 - 2D array (dim, num_pts) of grad V0
//----------------------------------------------------------------
// Splits the gradient into its 3 parts: impulse, X0, and V0
//----------------------------------------------------------------
void RegressionAcceleration::SplitGradient(const Array3D<double>& G, Array3D<double>& gradImpulse, Array2D<double>& gradX0, Array2D<double>& gradV0)
{
    int dim = G.GetLength();
    int nx = G.GetWidth();
    int T = this->_source.GetT();

    int Tplus2 = G.GetHeight();

    if ((T+2) != Tplus2)
    {
        // There is an error in usage here
        printf("Warning: Gradient is wrong dimension to hold impulse, x0, and v0. Application may crash or produce strange results.\n");
    }

    // Split impulse and v0
    //#pragma omp parallel for collapse(3)
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            for (int k=0; k<T+2; k++)
            {
                // Initial velocity
                if (k==T)
                {
                    gradV0(i,j) = G(i,j,k);
                }
                else if (k==(T+1))
                {
                    gradX0(i,j) = G(i,j,k);
                }
                // Impulse
                else
                {
                    gradImpulse(i,j,k) = G(i,j,k);
                }
            }
        }
    }

}

//----------------------------------------------------------------
// ComputeTrajectories
//----------------------------------------------------------------
// Inputs:
//   impulse - 3D array (dim, num_pts, time) of impulse vectors
//   v0 - 2D array (dim, num_pts) of initial velocity
//
// Outputs:
//----------------------------------------------------------------
// Apply the deformation to the source points this->_X using 
// a verlet integration scheme
//----------------------------------------------------------------
void RegressionAcceleration::ComputeTrajectories(const Array3D<double>& impulse, const Array2D<double>& v0)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();
    double tau2 = pow(this->_tau,2);

    // Set the initial velocity
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            // Should we use the v0 we computed or a v0 the user provided
            if (this->_source.ShouldUseInitV0())
            {
                // User provided initial velocity
                this->_dX(i,j,0) = this->_source.GetInitV0()(i,j);
            }
            else
            {
                // Use the value we computed
                this->_dX(i,j,0) = v0(i,j);
            }
        }
    }

    // Compute trajectories by Verlet integration
    for (int t=0; t<T-1; t++)
    {
        // Compute acceleration at time t
        Array2D<double> accelt;

        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            accelt = this->GridKernelSum(t, impulse);
        }
        else
        {
            accelt = this->RegKernelSum(t, impulse);
        }

        // Use the acceleration to update the positions
        //#pragma omp parallel for collapse(2)
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                this->_X(i,j,t+1) = this->_X(i,j,t) + this->_dX(i,j,t)*this->_tau + 0.5f*accelt(i,j)*tau2;
            }
        }

        // Compute acceleration at time t+1
        Array2D<double> acceltp1;

        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            // Update grids since point positions have changed
            this->UpdateSourceGrids();
            this->UpdateTargetGrids();

            acceltp1 = this->GridKernelSum(t+1, impulse);
        }
        else
        {
            acceltp1 = this->RegKernelSum(t+1, impulse);
        }

        // Now update the velocities
        //#pragma omp parallel for collapse(2)
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                this->_dX(i,j,t+1) = this->_dX(i,j,t) + 0.5f*(accelt(i,j)+acceltp1(i,j))*this->_tau;
            }
        }
    }
}

//----------------------------------------------------------------
// ComputeTrajectories
//----------------------------------------------------------------
// Inputs:
//   impulse - 3D array (dim, num_pts, time) of impulse vectors
//   v0 - 2D array (dim, num_pts) of initial velocity
//
// Outputs:
//----------------------------------------------------------------
// Apply the deformation to the source points this->_X using
// a verlet integration scheme
//----------------------------------------------------------------
void RegressionAcceleration::ComputeTrajectories(const Array3D<double>& impulse, const Array2D<double>& x0, const Array2D<double>& v0)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();

    if (this->_source.ShouldEstimateBaseline())
    {
        // Set the initial velocity
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                // Use the value we computed
                this->_X(i,j,0) = x0(i,j);
            }
        }

        // Update grids since point positions have changed
        this->UpdateSourceGrids();
        this->UpdateTargetGrids();
    }

    this->ComputeTrajectories(impulse, v0);

}

//----------------------------------------------------------------
// WriteSelf
//----------------------------------------------------------------
// Inputs:
//
//
// Outputs:
//----------------------------------------------------------------
// Writes shapes and vectors
//----------------------------------------------------------------
void RegressionAcceleration::WriteSelf()
{

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// OWN OPTIMIZATION METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// FISTA
//----------------------------------------------------------------
// Inputs:
//
//
// Outputs:
//----------------------------------------------------------------
// Fast iterative shrinkage thresholding algorithm
//----------------------------------------------------------------
void RegressionAcceleration::FISTA(Array3D<double>& impulse, Array2D<double>& X0, Array2D<double>& V0)
{
    int maxIters = this->_source.GetMaxIters();
    double breakRatio = this->_source.GetBreakRatio();
    int maxLineIters = 40;
    double stepImpulse = 0.05f;
    double stepX0andV0 = 0.05f;
    double stepIncrease = 1.2f;
    double stepDecrease = 0.25f;

    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    double* dataMatchingValues = new double[maxIters];
    double* regularityValues = new double[maxIters];
    double* criterionValues = new double[maxIters];

    // Compute the functional by following the algorithm's pointer
    criterionValues[0] = this->ComputeFunctional(impulse, X0, V0, dataMatchingValues[0], regularityValues[0]);
    double leastSquaresRef = criterionValues[0];

    Array3D<double> impulsePrev = impulse;
    Array2D<double> V0Prev = V0;
    Array2D<double> X0Prev = X0;

    Array3D<double> impulseTest(dim, nx, T);
    Array2D<double> V0Test(dim, nx);
    Array2D<double> X0Test(dim, nx);

    double tau = 1.0;
    int freezeDirectionCounter = 0;
    int numberItersFreeze = 3;

    Shape4DState::SetHasConverged(false);
    printf("\n");

    for (int iter=0; iter<maxIters; iter++)
    {
        // If it is time to save intermediate shapes
        if (((iter)%Shape4DState::GetSaveProgressEveryNIterations() == 0) && (Shape4DState::GetShouldSave()) && (iter != 0))
        {
            Shape4DState::SetShouldWriteShapesAndVectors(true);
        }

        // Compute the gradient
        Array3D<double> gradImpulse(dim, nx, T);
        Array2D<double> gradX0(dim, nx);
        Array2D<double> gradV0(dim, nx);
        this->ComputeGradient(impulse, X0, V0, gradImpulse, gradX0, gradV0);

        // Compute a good initial stepsize
        if (iter == 0)
        {
            // The initial step size is the functional value over the sum of squared gradient
            double sumSquared = 0.0f;
            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<nx; j++)
                {
                    for (int k=0; k<T; k++)
                    {
                        sumSquared += gradImpulse(i,j,k)*gradImpulse(i,j,k);
                    }
                }
            }
            stepImpulse = (criterionValues[0]/sumSquared)*0.01;
            stepX0andV0 = stepImpulse;
        }

        // Update the iteration values for saving
        Shape4DState::UpdateIteration(iter, dataMatchingValues[iter], regularityValues[iter], stepImpulse);
        printf("Iteration %3d   funct = %0.4f   data = %0.4f   reg = %0.4f   stepImpulse = %0.10lf   stepX0AndV0 = %0.10lf\n",
               iter+1, criterionValues[iter], dataMatchingValues[iter], regularityValues[iter], stepImpulse, stepX0andV0);

        bool minimTest = false;
        unsigned int lineIter = 0;
        for (; lineIter < maxLineIters; lineIter++)
        {
            // Apply the current gradient with current stepsize
            this->GradientDescentStep(impulseTest, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);

            criterionValues[iter+1] = this->ComputeFunctional(impulseTest, X0Test, V0Test, dataMatchingValues[iter+1], regularityValues[iter+1]);

            double QDiffTerm = this->QDiffTerm(impulseTest, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);
            double Q = leastSquaresRef - criterionValues[iter+1] + QDiffTerm;

            //printf("QDiffTerm = %0.4f\n", QDiffTerm);

            if (criterionValues[iter+1] > 1.10*criterionValues[iter])
            {
                Q = -1;
            }

            if (Q >= 0)
            {
                minimTest = true;
            }
            else
            {
                //printf("Negative Q, beginning search...\n");

                double QArray[4];

                // Case 1
                stepImpulse *= stepDecrease;
                Array3D<double> impulsePrime1(dim, nx, T);
                this->GradientDescentStep(impulsePrime1, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);

                double dataMatching1, regularity1;
                double functCase1 = this->ComputeFunctional(impulsePrime1, X0Test, V0Test, dataMatching1, regularity1);
                double QDiffTermCase1 = this->QDiffTerm(impulsePrime1, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);
                QArray[0] = leastSquaresRef - functCase1 + QDiffTermCase1;

                // Case 2
                stepImpulse *= stepDecrease;
                Array3D<double> impulsePrime2(dim, nx, T);
                this->GradientDescentStep(impulsePrime2, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);

                double dataMatching2, regularity2;
                double functCase2 = this->ComputeFunctional(impulsePrime2, X0Test, V0Test, dataMatching2, regularity2);
                double QDiffTermCase2 = this->QDiffTerm(impulsePrime2, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);
                QArray[1] = leastSquaresRef - functCase2 + QDiffTermCase2;

                // Case 3
                stepImpulse /= (stepDecrease*stepDecrease);
                stepX0andV0 *= stepDecrease;
                Array2D<double> X0Prime3(dim, nx);
                Array2D<double> V0Prime3(dim, nx);
                double dataMatching3, regularity3;
                double functCase3;
                double QDiffTermCase3;

                if ((this->_source.GetV0Weight() != 0) || (this->_source.ShouldEstimateBaseline() == true))
                {
                    //printf("Computing case 3 qdiff\n");
                    this->GradientDescentStep(impulseTest, X0Prime3, V0Prime3, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);
                    functCase3 = this->ComputeFunctional(impulseTest, X0Prime3, V0Prime3, dataMatching3, regularity3);
                    QDiffTermCase3 = this->QDiffTerm(impulseTest, X0Prime3, V0Prime3, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);
                    QArray[2] = leastSquaresRef - functCase3 + QDiffTermCase3;

                }
                else
                {
                    //printf("Skipping computing case 3 qdiff\n");
                    QArray[2] = -500000;
                }

                // Case 4
                stepX0andV0 *= stepDecrease;
                Array2D<double> X0Prime4(dim, nx);
                Array2D<double> V0Prime4(dim, nx);
                double dataMatching4, regularity4;
                double functCase4;
                double QDiffTermCase4;

                if ((this->_source.GetV0Weight() != 0) || (this->_source.ShouldEstimateBaseline() == true))
                {
                    //printf("Computing case 4 qdiff\n");
                    this->GradientDescentStep(impulseTest, X0Prime4, V0Prime4, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);
                    functCase4 = this->ComputeFunctional(impulseTest, X0Prime4, V0Prime4, dataMatching4, regularity4);
                    QDiffTermCase4 = this->QDiffTerm(impulseTest, X0Prime4, V0Prime4, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);
                    QArray[3] = leastSquaresRef - functCase4 + QDiffTermCase4;

                }
                else
                {
                    //printf("Skipping computing case 4 qdiff\n");
                    QArray[3] = -500000;
                }

                int indexMax = 0;
                double QMax = QArray[0];
                for (int i=1; i<4; i++)
                {
                    if (QArray[i] > QMax)
                    {
                        QMax = QArray[i];
                        indexMax = i;
                    }
                }

                //printf("%0.4f    %0.4f    %0.4f    %0.4f\n", QArray[0], QArray[1], QArray[2], QArray[3]);

                if (QMax >= 0)
                {
                    minimTest = true;

                    if (indexMax == 0)
                    {
                        //printf("Case 1\n");
                        impulseTest = impulsePrime1;

                        if (functCase1 > 1.10*criterionValues[iter])
                        {
                            stepImpulse *= stepDecrease;
                            stepX0andV0 *= stepDecrease;
                            continue;
                        }

                        criterionValues[iter+1] = functCase1;
                        dataMatchingValues[iter+1] = dataMatching1;
                        regularityValues[iter+1] = regularity1;

                        stepImpulse *= stepDecrease;
                        stepX0andV0 /= (stepDecrease*stepDecrease);

                    }
                    else if (indexMax == 1)
                    {
                        //printf("Case 2\n");
                        impulseTest = impulsePrime2;

                        if (functCase2 > 1.10*criterionValues[iter])
                        {
                            stepImpulse *= stepDecrease;
                            stepX0andV0 *= stepDecrease;
                            continue;
                        }

                        criterionValues[iter+1] = functCase2;
                        dataMatchingValues[iter+1] = dataMatching2;
                        regularityValues[iter+1] = regularity2;

                        stepImpulse *= (stepDecrease*stepDecrease);
                        stepX0andV0 /= (stepDecrease*stepDecrease);
                    }
                    else if (indexMax == 2)
                    {
                        //printf("Case 3\n");
                        X0Test = X0Prime3;
                        V0Test = V0Prime3;

                        if (functCase3 > 1.10*criterionValues[iter])
                        {
                            stepImpulse *= stepDecrease;
                            stepX0andV0 *= stepDecrease;
                            continue;
                        }

                        criterionValues[iter+1] = functCase3;
                        dataMatchingValues[iter+1] = dataMatching3;
                        regularityValues[iter+1] = regularity3;

                        stepX0andV0 /= stepDecrease;

                    }
                    else
                    {
                        //printf("Case 4\n");
                        X0Test = X0Prime4;
                        V0Test = V0Prime4;

                        if (functCase4 > 1.10*criterionValues[iter])
                        {
                            stepImpulse *= stepDecrease;
                            stepX0andV0 *= stepDecrease;
                            continue;
                        }

                        criterionValues[iter+1] = functCase4;
                        dataMatchingValues[iter+1] = dataMatching4;
                        regularityValues[iter+1] = regularity4;
                    }
                }
                // Update failed, continue line search
                else
                {
                    stepImpulse *= stepDecrease;
                    stepX0andV0 *= stepDecrease;
                }
            }

            if (minimTest == true)
            {
                break;
            }

        }    // End for (line search)

        if (minimTest == true)
        {
            double tauNext = 1.0;

            if ((lineIter == 0) && ((freezeDirectionCounter == 0) || (freezeDirectionCounter > numberItersFreeze)))
            {
                tauNext = (1 + sqrt(1.0 + 4.0 * tau * tau / stepIncrease)) / 2.0;
                if (freezeDirectionCounter > numberItersFreeze)
                {
                    freezeDirectionCounter = 0;
                }
            }
            else
            {
                tauNext = (1 + sqrt(1.0 + 4.0 * tau * tau)) / 2.0;
                freezeDirectionCounter++;
            }

            double tauScale = (tau - 1.0) / tauNext;

            impulse = impulseTest + (impulseTest - impulsePrev) * tauScale;
            X0 = X0Test + (X0Test - X0Prev) * tauScale;
            V0 = V0Test + (V0Test - V0Prev) * tauScale;

            impulsePrev = impulseTest;
            X0Prev = X0Test;
            V0Prev = V0Test;

            tau = tauNext;

            if ((freezeDirectionCounter == 0) || (freezeDirectionCounter > numberItersFreeze))
            {
                stepImpulse *= stepIncrease;
                if ((this->_source.GetV0Weight() != 0) || (this->_source.ShouldEstimateBaseline() == true))
                {
                    stepX0andV0 *= stepIncrease;
                }
            }

            leastSquaresRef = criterionValues[iter+1];

        }
        // We have exhausted line searching so we're done
        else
        {
            break;
        }

        // Check for regular old convergence (termination criteria)
        //
        double deltaFCurr = fabs(criterionValues[iter]-criterionValues[iter+1]);
        double deltaFRef = fabs(criterionValues[0]-criterionValues[iter+1]);
        if (deltaFCurr < (breakRatio*deltaFRef))
        {
            break;
        }

    }   // End for (iterations of FISTA)

    delete[] dataMatchingValues;
    delete[] regularityValues;
    delete[] criterionValues;

    Shape4DState::SetHasConverged(true);

}

//----------------------------------------------------------------
// GradientDescent
//----------------------------------------------------------------
// Inputs:
//
//
// Outputs:
//----------------------------------------------------------------
// Gradient descent with adaptive step sizes
//----------------------------------------------------------------
void RegressionAcceleration::GradientDescent(Array3D<double>& impulse, Array2D<double>& X0, Array2D<double>& V0)
{
    int maxIters = this->_source.GetMaxIters();
    double breakRatio = this->_source.GetBreakRatio();
    int maxLineIters = 40;
    double stepImpulse = 0.05f;
    double stepX0andV0 = 0.05f;
    double stepIncrease = 1.2f;
    double stepDecrease = 0.25f;

    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    double* dataMatchingValues = new double[maxIters];
    double* regularityValues = new double[maxIters];
    double* criterionValues = new double[maxIters];

    // Compute the functional by following the algorithm's pointer
    criterionValues[0] = this->ComputeFunctional(impulse, X0, V0, dataMatchingValues[0], regularityValues[0]);
    double leastSquaresRef = criterionValues[0];

    Array3D<double> impulsePrev = impulse;
    Array2D<double> V0Prev = V0;
    Array2D<double> X0Prev = X0;

    Array3D<double> impulseTest(dim, nx, T);
    Array2D<double> V0Test(dim, nx);
    Array2D<double> X0Test(dim, nx);

    Shape4DState::SetHasConverged(false);
    printf("\n");

    for (int iter=0; iter<maxIters; iter++)
    {
        // If it is time to save intermediate shapes
        if ((iter%Shape4DState::GetSaveProgressEveryNIterations() == 0) && (Shape4DState::GetShouldSave()))
        {
            Shape4DState::SetShouldWriteShapesAndVectors(true);
        }

        // Compute the gradient
        Array3D<double> gradImpulse(dim, nx, T);
        Array2D<double> gradX0(dim, nx);
        Array2D<double> gradV0(dim, nx);
        this->ComputeGradient(impulse, X0, V0, gradImpulse, gradX0, gradV0);

        // Compute a good initial stepsize
        if (iter == 0)
        {
            // The initial step size is the functional value over the sum of squared gradient
            double sumSquared = 0.0f;
            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<nx; j++)
                {
                    for (int k=0; k<T; k++)
                    {
                        sumSquared += gradImpulse(i,j,k)*gradImpulse(i,j,k);
                    }
                }
            }
            stepImpulse = (criterionValues[0]/sumSquared)*0.01;
            stepX0andV0 = stepImpulse;
        }

        // Update the iteration values for saving
        Shape4DState::UpdateIteration(iter, dataMatchingValues[iter], regularityValues[iter], stepImpulse);
        printf("Iteration %3d   funct = %0.4f   data = %0.4f   reg = %0.4f   stepImpulse = %0.10lf   stepX0AndV0 = %0.10lf\n",
               iter+1, criterionValues[iter], dataMatchingValues[iter], regularityValues[iter], stepImpulse, stepX0andV0);

        bool minimTest = false;
        unsigned int lineIter = 0;
        for (; lineIter < maxLineIters; lineIter++)
        {
            // Apply the current gradient with current stepsize
            this->GradientDescentStep(impulseTest, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);

            criterionValues[iter+1] = this->ComputeFunctional(impulseTest, X0Test, V0Test, dataMatchingValues[iter+1], regularityValues[iter+1]);

            double Q = leastSquaresRef - criterionValues[iter+1];

            //printf("%0.4f\n", Q);
            //printf("%0.10lf     %0.10lf\n", stepImpulse, stepX0andV0);

            if (Q >= 0)
            {
                minimTest = true;
                break;
            }
            else
            {
                // Case 1
                stepImpulse *= stepDecrease;
                Array3D<double> impulsePrime1(dim, nx, T);
                this->GradientDescentStep(impulsePrime1, X0Test, V0Test, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);

                double dataMatching1, regularity1;
                double functCase1 = this->ComputeFunctional(impulsePrime1, X0Test, V0Test, dataMatching1, regularity1);
                double Q1 = leastSquaresRef - functCase1;

                // Case 2
                stepImpulse /= stepDecrease;
                stepX0andV0 *= stepDecrease;
                Array2D<double> X0Prime2(dim, nx);
                Array2D<double> V0Prime2(dim, nx);
                this->GradientDescentStep(impulseTest, X0Prime2, V0Prime2, impulse, X0, V0, gradImpulse, gradX0, gradV0, stepImpulse, stepX0andV0);

                double dataMatching2, regularity2;
                double functCase2 = this->ComputeFunctional(impulseTest, X0Prime2, V0Prime2, dataMatching2, regularity2);
                double Q2 = leastSquaresRef - functCase2;

                if ( (Q1 >= 0) || (Q2 >= 0) )
                {
                    if (Q1 >= Q2)
                    {
                        impulseTest = impulsePrime1;

                        criterionValues[iter+1] = functCase1;
                        dataMatchingValues[iter+1] = dataMatching1;
                        regularityValues[iter+1] = regularity1;

                        stepImpulse *= stepDecrease;
                        stepX0andV0 /= (stepDecrease);
                    }
                    else
                    {
                        X0Test = X0Prime2;
                        V0Test = V0Prime2;

                        criterionValues[iter+1] = functCase2;
                        dataMatchingValues[iter+1] = dataMatching2;
                        regularityValues[iter+1] = regularity2;
                    }

                    minimTest = true;

                }
                // Update failed, continue line search
                else
                {
                    stepImpulse *= stepDecrease;
                    stepX0andV0 *= stepDecrease;
                }
            }


            if (minimTest == true)
            {
                break;
            }

        }   // End line search

        if (minimTest == true)
        {
            impulse = impulseTest;
            X0 = X0Test;
            V0 = V0Test;

            stepImpulse *= stepIncrease;
            if ((this->_source.GetV0Weight() != 0) || (this->_source.ShouldEstimateBaseline() == true))
            {
                stepX0andV0 *= stepIncrease;
            }

            leastSquaresRef = criterionValues[iter+1];

        }
        // We have exhausted line searching so we're done
        else
        {
            break;
        }

        // Check for regular old convergence (termination criteria)
        //
        double deltaFCurr = fabs(criterionValues[iter]-criterionValues[iter+1]);
        double deltaFRef = fabs(criterionValues[0]-criterionValues[iter+1]);
        if (deltaFCurr < (breakRatio*deltaFRef))
        {
            break;
        }

    }   // End main iters loop


    delete[] dataMatchingValues;
    delete[] regularityValues;
    delete[] criterionValues;

    Shape4DState::SetHasConverged(true);
}

//----------------------------------------------------------------
// GradientDescentStep
//----------------------------------------------------------------
// Inputs:
//   impulseTest -
//   X0Test -
//   V0Test -
//   impulse -
//   X0 -
//   V0 -
//   G - gradient of the criterion (wrt impulse, v0, and x0)
//   stepImpulse -
//   stepX0andV0 -
//
//
// Outputs:
//----------------------------------------------------------------
// Takes a step in the gradient direction
//----------------------------------------------------------------
void RegressionAcceleration::GradientDescentStep(Array3D<double>& impulseTest, Array2D<double>& X0Test, Array2D<double>& V0Test,
                                                 const Array3D<double>& impulse, const Array2D<double>& X0, const Array2D<double>& V0,
                                                 const Array3D<double>& gradImpulse, const Array2D<double>&, const Array2D<double>& gradV0,
                                                 double stepImpulse, double stepX0andV0)
{
    // Update the impulse
    impulseTest = impulse - gradImpulse*stepImpulse;

    // Update V0
    if ((this->_source.ShouldUseInitV0()) || (this->_source.GetV0Weight() == 0.0))
    {
        V0Test = V0;
    }
    else
    {
        V0Test = V0 - gradV0*stepX0andV0;
    }

    // Update X0
    if (this->_source.ShouldEstimateBaseline())
    {
        X0Test = X0 - this->_sobolevGrad*stepX0andV0;
        //X0Test = X0 - gradX0*stepX0andV0;
    }
    else
    {
        X0Test = X0;
    }
}

double RegressionAcceleration::QDiffTerm(const Array3D<double> &impulseTest, const Array2D<double> X0Test, const Array2D<double> V0Test,
                                         const Array3D<double> &impulse, const Array2D<double> X0, const Array2D<double> V0,
                                         const Array3D<double> &gradImpulse, const Array2D<double> &gradX0, const Array2D<double> &gradV0,
                                         double stepImpulse, double stepX0andV0)
{
    int dim = impulse.GetLength();
    int nx = impulse.GetWidth();
    int T = impulse.GetHeight();

    Array3D<double> impulseDiff(dim, nx, T);
    Array2D<double> X0Diff(dim, nx);
    Array2D<double> V0Diff(dim, nx);

    impulseDiff = impulseTest - impulse;
    X0Diff = X0Test - X0;
    V0Diff = V0Test - V0;

    // Term 1 is the sum of all the dot products with gradients
    double termVal1 = 0.0;
    // Term 2 is norm of impulse difference
    double termVal2 = 0.0;
    // Term 3 is the sum of norm of x0 difference and v0 difference
    double termVal3 = 0.0;

    for (int i=0; i<nx; i++)
    {
        double diffX0X = X0Diff(0, i);  double diffX0Y = X0Diff(1, i);  double diffX0Z = X0Diff(2, i);
        double gradX0X = gradX0(0, i);  double gradX0Y = gradX0(1, i);  double gradX0Z = gradX0(2, i);

        termVal1 += diffX0X*gradX0X + diffX0Y*gradX0Y + diffX0Z*gradX0Z;
        termVal3 += diffX0X*diffX0X + diffX0Y*diffX0Y + diffX0Z*diffX0Z;

        double diffV0X = V0Diff(0, i);  double diffV0Y = V0Diff(1, i);  double diffV0Z = V0Diff(2, i);
        double gradV0X = gradV0(0, i);  double gradV0Y = gradV0(1, i);  double gradV0Z = gradV0(2, i);

        termVal1 += diffV0X*gradV0X + diffV0Y*gradV0Y + diffV0Z*gradV0Z;
        termVal3 += diffV0X*diffV0X + diffV0Y*diffV0Y + diffV0Z*diffV0Z;

        for (int j=0; j<T; j++)
        {
            double diffImpulseX = impulseDiff(0, i, j);  double diffImpulseY = impulseDiff(1, i, j);  double diffImpulseZ = impulseDiff(2, i, j);
            double gradImpulseX = gradImpulse(0, i, j);  double gradImpulseY = gradImpulse(1, i, j);  double gradImpulseZ = gradImpulse(2, i, j);

            termVal1 += diffImpulseX*gradImpulseX + diffImpulseY*gradImpulseY + diffImpulseZ*gradImpulseZ;
            termVal2 += diffImpulseX*diffImpulseX + diffImpulseY*diffImpulseY + diffImpulseZ*diffImpulseZ;
        }
    }

    //printf("Term 1 = %0.4f\n", termVal1);
    //printf("Term 2 = %0.4f\n", termVal2);
    //printf("Term 3 = %0.4f\n", termVal3);

    double finalVal = termVal1 + termVal2 / (2.0*stepImpulse) + termVal3 / (2.0*stepX0andV0);


    return finalVal;
}
