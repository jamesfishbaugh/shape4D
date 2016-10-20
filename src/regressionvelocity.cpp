//--------------------------------------------------------------------------
// EXRegressionVelocity:  Shape regression parameterized by velocity.
//--------------------------------------------------------------------------

#include "regressionvelocity.h"

#include <stdio.h>		// For printf
#include <math.h>		// For pow
#include <stdlib.h>		// For free

#include "optimizer.h"					// The optimizer
#include "adaptivegradientdescent.h"    // Specific implementation of optimizer
#include "shape4dstate.h"				// To save intermediate results during optimization

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//---------------------------------------------------------------------------F-------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXRegressionVelocity
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
RegressionVelocity::RegressionVelocity()
{
    this->_iteration = 0;
    this->SetNumTargets(0);
    this->_targets = 0;
    this->_continueRegression = false;
}

//----------------------------------------------------------------
// EXRegressionVelocity
//----------------------------------------------------------------
// Inputs: 
//   source - reference to a source object
//
// Outputs:
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
RegressionVelocity::RegressionVelocity(const RegressionParams& source)
{
    this->_iteration = 0;
    this->SetSource(source);
    this->SetNumTargets(0);
    this->_targets = 0;
    this->_continueRegression = false;
}

//----------------------------------------------------------------
// EXRegressionVelocity
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
RegressionVelocity::RegressionVelocity(const RegressionParams& source, int numTargets, TargetData** targets)
{
    this->_iteration = 0;
    this->SetSource(source);
    this->SetNumTargets(numTargets);
    this->SetTargets(targets);
    this->_continueRegression = false;
}

//----------------------------------------------------------------
// ~EXRegressionVelocityAccel
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
RegressionVelocity::~RegressionVelocity()
{
    delete [] this->_targets;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

void RegressionVelocity::SetMomenta(Array3D<double> &momenta)
{
    this->_momenta = momenta;
    this->_continueRegression = true;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetV0
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - a 2D array (dim, num_pts) of initial velocity
//----------------------------------------------------------------
// Returns initial velocity
//----------------------------------------------------------------
Array2D<double> RegressionVelocity::GetV0()
{
    if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
    {
        this->UpdateSourceGrids();
        this->UpdateTargetGrids();

        return this->GridKernelSum(0, this->_momenta);
    }
    else
    {
        return this->RegKernelSum(0, this->_momenta);
    }
}

//----------------------------------------------------------------
// GetVelocity
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 3D array (dim, num_pts) of velocity
//----------------------------------------------------------------
// Returns velocity vectors over time
//----------------------------------------------------------------
Array3D<double> RegressionVelocity::GetVelocity()
{
    return this->_dX;
}

//----------------------------------------------------------------
// GetMomenta
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 3D array (dim, num_pts) of momenta
//----------------------------------------------------------------
// Returns momenta vectors over time
//----------------------------------------------------------------
Array3D<double> RegressionVelocity::GetMomenta()
{
    return this->_momenta;
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
Array3D<double> RegressionVelocity::Run()
{
    //printf("Beginning of run\n");

    this->_tau = (double)1.0f/(double)(this->_source.GetT()-1.0f);

    // Initialize member arrays
    this->InitX();

    //printf("Done initializing X\n");

    if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
    {
        // Create the grids for optimization
        this->SetSourceGrids();
        this->SetTargetGrids();
    }

    //printf("Done creating target grid\n");

    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    // Allocate memory for momenta
    Array3D<double> momenta(dim, nx, T);
    momenta.FillArray(0.0f);

    if (this->_continueRegression)
    {
        momenta = this->_momenta;
    }

    //printf("Starting optimizer...\n");

    // Run the optimization
    Optimizer* optimizer = (Optimizer*) new AdaptiveGradientDescent(this);
    optimizer->SetMaxIterations(this->_source.GetMaxIters());
    optimizer->SetBreakRatio(this->_source.GetBreakRatio());
    optimizer->Optimize(momenta);

    // Compute the final trajectories now that we are done
    this->ComputeTrajectories(momenta);

    // Compute the final velocites
    Array3D<double> dX(dim, nx, T);
    // Compute trajectories by Verlet integration
    for (int t=0; t<T; t++)
    {
        Array2D<double> dXt;

        // Compute velocity at time t
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            dXt = this->GridKernelSum(t, momenta);
        }
        else
        {
            dXt = this->RegKernelSum(t, momenta);
        }

        dX.Set2DSliceAtHeight(dXt,t);
    }
    // Save the final values
    this->_dX = dX;
    this->_momenta = momenta;

    if (Shape4DState::GetShouldSave())
    {
        // Save the shapes and vectors
        vector<char *> names;
        names.push_back((char*)"momenta");
        names.push_back((char*)"velocity");

        vector< Array3D<double> > vectors;
        vectors.push_back(this->_momenta);
        vectors.push_back(this->_dX);

        Shape4DState::SetX(_X);
        Shape4DState::SetVectors(names, vectors);
        Shape4DState::SaveShapesAndVectors();
    }

    delete optimizer;

    // Return the final trajectories
    return this->_X;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// COMPUTATION OF REGRESSION FUNCTIONAL
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeRegularity
//----------------------------------------------------------------
// Inputs:
//   momenta - 3D array (dim, num_pts, time) of momenta vectors
//
// Outputs:
//   return - regularity of deformation
//----------------------------------------------------------------
// Computes the regularity of the deformation, interpreted as 
// the total kinetic energy of the system
//----------------------------------------------------------------	
double RegressionVelocity::ComputeRegularity(const Array3D<double>& momenta)
{
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();
    double E = 0.0f;

    for (int t=0; t<T-1; t++)
    {
        for (int i=0; i<nx; i++)
        {
            double dot = (momenta(0,i,t)+momenta(0,i,t+1))*(this->_X(0,i,t+1)-this->_X(0,i,t)) +
                    (momenta(1,i,t)+momenta(1,i,t+1))*(this->_X(1,i,t+1)-this->_X(1,i,t)) +
                    (momenta(2,i,t)+momenta(2,i,t+1))*(this->_X(2,i,t+1)-this->_X(2,i,t));
            E += dot;
        }
    }

    // Because we use momenta(t)+momenta(t+1) above
    return 0.5f*E;
}

//----------------------------------------------------------------
// ComputeFunctional
//----------------------------------------------------------------
// Inputs:
//   moment - 3D array (dim, num_pts, time) of momenta vectors
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
double RegressionVelocity::ComputeFunctional(const Array3D<double>& momenta, double &dataMatching, double &regularity)
{
    //printf("  Computing functional value\n");
    //printf("    Computing trajectories...\n");

    // Compute trajectories
    this->ComputeTrajectories(momenta);

    //printf("    Computing data matching...\n");
    dataMatching = ComputeDataMatching(this->_X);
    //printf("    Computing regularity...\n");
    regularity = ComputeRegularity(momenta);

    //printf("Data = %0.4f  Reg = %0.4f\n", dataMatching, regularity);

    return dataMatching + this->_source.GetGamma()*regularity;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// COMPUTATION OF GRADIENT OF REGRESSION FUNCTIONAL
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeGradient
//----------------------------------------------------------------
// Inputs:
//   momenta - 3D array (dim, num_pts, time) of momenta vectors
// Outputs:
//   return - 3D array (dim, num_pts, time) of the gradient of
//            the regression functional
//----------------------------------------------------------------
// Computes the gradient of regression functional.  Used to 
// update impulse and initial velocity.
//----------------------------------------------------------------
Array3D<double> RegressionVelocity::ComputeGradient(const Array3D<double>& momenta)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();
    double sigmaV2 = pow(this->_source.GetSigmaV(),2);

    //printf("  Computing gradient of functional\n");
    //printf("    Computing trajectories...\n");

    // Compute trajectories
    this->ComputeTrajectories(momenta);

    // Save shape and vector data when needed
    if (Shape4DState::GetShouldWriteShapesAndVectors())
    {
        this->_momenta = momenta;
        this->WriteSelf();
        Shape4DState::SetShouldWriteShapesAndVectors(false);
    }

    //printf("    Computing gradient...\n");

    // The sum of etaAux represents the contribution of the sum of grad A
    Array3D<double> etaAux = this->ComputeDataMatchingGradient(this->_X);
    double sumetaAux = etaAux.Sum();

    // Create etax and etaxdot
    Array3D<double> eta(dim,nx,T);
    eta.FillArray(0.0f);

    // Create sumdeta
    Array2D<double> sumdeta(dim, nx);
    sumdeta.FillArray(0.0f);

    // Create detat and detatm1
    Array2D<double> detat(dim,nx);
    Array2D<double> auxt(dim, nx);
    Array2D<double> detatm1(dim,nx);
    Array2D<double> auxtm1(dim, nx);

    //printf("    Starting etax and etaxdot...\n");

    // Computes etax[t] and etaxdot[t] with backward integration centered scheme
    for (int t=T-1; t>=1; t--)
    {
        //sumetaAux = etaAux->Sum();
        // If the sum of etaAux is non-zero
        if (sumetaAux != 0)
        {
            // Update the sum of deta
            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<nx; j++)
                {
                    sumdeta(i,j) += etaAux(i,j,t);
                    // Current value of etax is equal to this sum
                    eta(i,j,t) = sumdeta(i,j);
                }
            }
        }

        // Compute the summation inside the integral from equation 4.3.17 from thesis (deta)
        // Compute derivaitve of eta at time t
        detat.FillArray(0.0f);
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            detat = this->GridCompDeta(t, eta, momenta);
        }
        else
        {
            detat = this->RegCompDeta(t, eta, momenta);
        }

        auxt.FillArray(0.0f);

        //printf("    Done with GridCompDeta...\n");

        //  A small step in the direction of deta
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                double auxtValue = (this->_tau*2/sigmaV2) * detat(i,j);
                auxt(i,j) = auxtValue;
                double etaxUpdate = eta(i,j,t) + auxt(i,j);
                eta(i,j,t-1) = etaxUpdate;
                double sumdetaValue = sumdeta(i,j) + auxt(i,j);
                sumdeta(i,j) = sumdetaValue;
            }
        }

        // Compute the summation inside the integral from equation 4.3.17 from thesis (deta)
        // Compute derivaitve of eta at time t-1
        detatm1.FillArray(0.0f);
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            detatm1 = this->GridCompDeta(t-1, eta, momenta);
        }
        else
        {
            detatm1 = this->RegCompDeta(t-1, eta, momenta);
        }

        for (int i=0; i<3; i++)
        {
            for (int j=0; j<nx; j++)
            {
                detatm1(i,j) = detatm1(i,j) + detat(i,j);
            }
        }

        auxtm1.FillArray(0.0f);

        //  A small step in the direction of deta
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                double auxtValue = (this->_tau/sigmaV2) * detatm1(i,j);
                auxtm1(i,j) = auxtValue;
                double etaxUpdate = eta(i,j,t) + auxtm1(i,j);
                eta(i,j,t-1) = etaxUpdate;
                double sumdetaValue =  sumdeta(i,j) + auxtm1(i,j);
                sumdeta(i,j) = sumdetaValue;
            }
        }
    }

    Array3D<double> G(dim,nx,T);
    G.FillArray(0.0f);

    double gammaTimesTwo = 2.0f*this->_source.GetGamma();

    // Compute the final value of the gradient
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            for (int k=0; k<T; k++)
            {
                G(i,j,k) = gammaTimesTwo*momenta(i,j,k) + eta(i,j,k);
            }
        }
    }

    return G;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeTrajectories
//----------------------------------------------------------------
// Inputs:
//   momenta - 3D array (dim, num_pts, time) of momenta vectors
//
// Outputs:
//----------------------------------------------------------------
// Apply the deformation to the source points this->_X using
// an Euler integration scheme with prediction/correction
//----------------------------------------------------------------
void RegressionVelocity::ComputeTrajectories(const Array3D<double>& momenta)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    for (int t=0; t<T-1; t++)
    {
        Array2D<double> dXt;

        // Compute acceleration at time t
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            dXt = this->GridKernelSum(t, momenta);
        }
        else
        {
            dXt = this->RegKernelSum(t, momenta);
        }

        // Use the acceleration to update the positions
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                this->_X(i,j,t+1) = this->_X(i,j,t) + dXt(i,j)*this->_tau;
            }
        }

        Array2D<double> dXtp1;
        // Compute acceleration at time t+1
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            this->UpdateSourceGrids();
            this->UpdateTargetGrids();

            dXtp1 = this->GridKernelSum(t+1, momenta);
        }
        else
        {
            dXtp1 = this->RegKernelSum(t+1, momenta);
        }

        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                dXtp1(i,j) = dXtp1(i,j) + dXt(i,j);
            }
        }

        // Now update the velocities
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                this->_X(i,j,t+1) = this->_X(i,j,t) + 0.5f*dXtp1(i,j)*this->_tau;
            }
        }

        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            this->UpdateSourceGrids();
            this->UpdateTargetGrids();
        }
    }
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
void RegressionVelocity::WriteSelf()
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    // Compute the final velocites
    Array3D<double> dX(dim, nx, T);
    // Compute trajectories by Verlet integration
    for (int t=0; t<T; t++)
    {
        Array2D<double> dXt;

        // Compute acceleration at time t
        if (strcmp(this->_source.GetKernelType(), "p3m") == 0)
        {
            dXt = this->GridKernelSum(t, this->_momenta);
        }
        else
        {
            dXt = this->RegKernelSum(t, this->_momenta);
        }
        dX.Set2DSliceAtHeight(dXt,t);
    }

    // Save the final values
    this->_dX = dX;

    // Save the shapes and vectors
    vector<char *> names;
    names.push_back((char*)"momenta");
    names.push_back((char*)"velocity");

    vector< Array3D<double> > vectors;
    vectors.push_back(this->_momenta);
    vectors.push_back(this->_dX);

    Shape4DState::SetX(_X);
    Shape4DState::SetVectors(names, vectors);
    Shape4DState::SaveShapesAndVectors();
    Shape4DState::SaveStateVelocity();
}
