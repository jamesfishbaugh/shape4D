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
    this->FISTA(graddesc);
    Optimizer* optimizer = (Optimizer*) new AdaptiveGradientDescent(this);
	//optimizer->SetStepsize(EXExoshapeState::GetCurrentStepsize());
    optimizer->SetMaxIterations(this->_source.GetMaxIters());
    optimizer->SetBreakRatio(this->_source.GetBreakRatio());
    optimizer->SetStepsize(0.05f);
	optimizer->Optimize(graddesc);
	
	// Split out the impulse and v0
    Array3D<double> impulse(dim, nx, T);
    Array2D<double> v0(dim, nx);
	this->SplitImpulseAndV0(graddesc, impulse, v0);

	// Compute the final trajectories now that we are done
	this->ComputeTrajectories(impulse, v0);

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
	
	// Clean up memory
	delete optimizer;
	
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
	
	double gammaTimesTwo = 2.0f*this->_source.GetGamma();
	
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
void RegressionAcceleration::ComputeGradient(const Array3D<double>& impulse, const Array2D<double>& X0, const Array2D<double>& V0,
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

    // We need to ensure the baseline shape doesn't self intersect by convolving the gradient with a kernel
    Array2D<double> sobolevGradient(dim,nx);
    sobolevGradient.FillArray(0.0f);

    if (this->_source.ShouldEstimateBaseline())
    {
        //#pragma omp parallel for collapse(2)
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nx; j++)
            {
                double argin = - ( pow(this->_X(0,i,0)-this->_X(0,j,0),2) +
                                   pow(this->_X(1,i,0)-this->_X(1,j,0),2) +
                                   pow(this->_X(2,i,0)-this->_X(2,j,0),2) ) / sigmaV2;

                double argout = exp(argin);

                sobolevGradient(0,i) = sobolevGradient(0,i) + argout*gradX0(0,j);
                sobolevGradient(1,i) = sobolevGradient(1,i) + argout*gradX0(1,j);
                sobolevGradient(2,i) = sobolevGradient(2,i) + argout*gradX0(2,j);
            }
        }

        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<nx; j++)
            {
                gradX0(i,j) = sobolevGradient(i,j);
            }
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
void RegressionAcceleration::FISTA(const Array3D<double>& impulseAndV0)
{
    int maxIters = this->_source.GetMaxIters();
    double breakRatio = this->_source.GetBreakRatio();
    int maxLineIters = 40;
    double stepImpulse = 0.05f;
    double stepX0andV0 = 0.05f;

    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();
    double sigmaV2 = pow(this->_source.GetSigmaV(),2);

    // Split out impulse and initial velocity
    Array3D<double> impulse(dim,nx,T);
    Array2D<double> V0 (dim,nx);
    this->SplitImpulseAndV0(impulseAndV0, impulse, V0);
    Array2D<double> X0 (dim,nx);
    X0 = this->_X.Get2DSliceAtHeight(0);

    for (int iter=1; iter<maxIters; iter++)
    {
        bool minimTest = false;
        unsigned int lineIter = 0;
        for (; lineIter < maxLineIters; lineIter++)
        {

        }
    }

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
                                                   double stepImpulse, double stepX0andV0)
{


}
