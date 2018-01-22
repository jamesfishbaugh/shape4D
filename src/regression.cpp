//--------------------------------------------------------------------------
// EXRegression:  Virutal base class that implements common functionality
// for regression
//--------------------------------------------------------------------------

#include <stdio.h>		// For printf
#include <math.h>		// For pow
#include <stdlib.h>		// For free

#include "regression.h"

#include "gridoptimize.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// DESTRUCTOR
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ~EXRegression
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
Regression::~Regression()
{

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetTheSource
//----------------------------------------------------------------
// Inputs:
//   theSource - reference to a source object
//
// Outputs:
//----------------------------------------------------------------
// Setter for the source
//----------------------------------------------------------------
void Regression::SetTheSource(const MultiObjectComplex theSource)
{
    MultiObjectComplex newSource(theSource);
    this->_theSource = newSource;
}

//----------------------------------------------------------------
// SetTargets
//----------------------------------------------------------------
// Inputs:
//   theTargets - vector of target objects
//
// Outputs:
//----------------------------------------------------------------
// Setter for the target vector
//----------------------------------------------------------------
void Regression::SetTargets(const vector<MultiObjectComplex>& theTargets)
{
    vector<MultiObjectComplex> newTargets(theTargets);
    this->_theTarget = newTargets;
}

//----------------------------------------------------------------
// AddTarget
//----------------------------------------------------------------
// Inputs:
//   target - reference to a target object
//
// Outputs:
//----------------------------------------------------------------
// Adds an additional target to the current target vector
//----------------------------------------------------------------
void Regression::AddTarget(const MultiObjectComplex& target)
{
    MultiObjectComplex newTarget(target);
    this->_theTarget.push_back(newTarget);
}

//----------------------------------------------------------------
// SetSource
//----------------------------------------------------------------
// Inputs:
//   source - reference to a source object
//
// Outputs:
//----------------------------------------------------------------
// Setter for the source
//----------------------------------------------------------------
void Regression::SetSource(const RegressionParams& source)
{
    this->_source = source;
}

//----------------------------------------------------------------
// SetNumTargets
//----------------------------------------------------------------
// Inputs:
//   numTargets - number of target objects
//
// Outputs:
//----------------------------------------------------------------
// Setter for the number of targets
//----------------------------------------------------------------
void Regression::SetNumTargets(int numTargets)
{
    this->_numTargets = numTargets;
}

//----------------------------------------------------------------
// SetTargets
//----------------------------------------------------------------
// Inputs:
//   targets - array of pointers to target objects
//
// Outputs:
//----------------------------------------------------------------
// Setter for target object array.  This method assumes the
// number of targets has previously been set.
//----------------------------------------------------------------
void Regression::SetTargets(TargetData** targets)
{
    this->SetTargets(this->_numTargets, targets);
}

//----------------------------------------------------------------
// SetTargets
//----------------------------------------------------------------
// Inputs:
//   numTargets - number of target objects
//   targets - array of pointers to target objects
//
// Outputs:
//----------------------------------------------------------------
// Setter for target object array.
//----------------------------------------------------------------
void Regression::SetTargets(int numTargets, TargetData** targets)
{
    this->_targets = new TargetData*[numTargets];

    for (int i=0; i<numTargets; i++)
    {
        this->_targets[i] = targets[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// INITIALIZATION METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// InitX
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Initializes the 3D array (dim, num_pts, time) this->_X which 
// stores point trajectories
//----------------------------------------------------------------
void Regression::InitX()
{
    Array2D<double> X = this->_source.GetX();
    int dim = X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    Array3D<double> tempX(dim,nx,T);
    tempX.FillArray(0.0f);

    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            tempX(i,j,0) = X(i,j);
        }
    }
    this->_X = tempX;
}

//--------------------------------------------------------------------------
// DATA MATCHING AND GRADIENT OF DATA MATCHING
//--------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeDataMatching
//----------------------------------------------------------------
// Inputs:
//   X - 3D array (dim, num_pts, time) of source points
//
// Outputs:
//   return - total "difference" between source points and all
//            target shapes
//----------------------------------------------------------------
// Computes and sums the shape difference between X and each 
// target.  Each target determines the metric used for matching.
//----------------------------------------------------------------	
double Regression::ComputeDataMatching(const Array3D<double>& X)
{
//    double sum = 0;

//    int numTargets = this->_theTarget.size();

//    for (int i=0; i<numTargets; i++)
//    {
//        MultiObjectComplex curShapeComplex = this->_theTarget[i];

//        for (int j=0; j<curShapeComplex.GetNumberOfShapes(); j++)
//        {
//            ShapeObject* curShapeObject = curShapeComplex.GetShapeAt(j);
//            int curTimeIndex = curShapeObject->GetTimeIndex();
//            Array2D<double> correspondingXt = X.Get2DSliceAtHeight(curTimeIndex);
//            sum += curShapeObject->GetWeight() * curShapeObject->Matching(correspondingXt);
//        }
//    }

//    return sum;

    double sum = 0;

    for (int i=0; i<this->_numTargets; i++)
    {
        int curTimept = this->_targets[i]->GetTimeIndex();
        sum += this->_targets[i]->GetWeight()*this->_targets[i]->Matching(X, curTimept);
    }

    return sum;
}

//----------------------------------------------------------------
// ComputeDataMatchingGradient
//----------------------------------------------------------------
// Inputs:
//   X - 3D array (dim, num_pts, time) of source points
//
// Outputs:
//   return - 3D array (dim, num_pts, time) of the gradient of
//            the data matching metric between the source
//            points and all targets
//----------------------------------------------------------------
// Computes the gradient of the data matching metric between X
// and all targets.  Each target determines the metric used for 
// matching.
//----------------------------------------------------------------
Array3D<double> Regression::ComputeDataMatchingGradient(const Array3D<double>& X)
{	
//    int dim = X.GetLength();
//    int nx = X.GetWidth();
//    int T = X.GetHeight();

//    Array3D<double> G(dim, nx, T);
//    G.FillArray(0.0);

//    int numTargets = this->_theTarget.size();

//    for (int i=0; i<numTargets; i++)
//    {
//        MultiObjectComplex curShapeComplex = this->_theTarget[i];

//        for (int j=0; j<curShapeComplex.GetNumberOfShapes(); j++)
//        {
//            ShapeObject* curShapeObject = curShapeComplex.GetShapeAt(j);
//            int curTimeIndex = curShapeObject->GetTimeIndex();
//            Array2D<double> correspondingXt = X.Get2DSliceAtHeight(curTimeIndex);

//            Array2D<double> curG = curShapeObject->GradMatching(correspondingXt);

//            for (int d=0; d<dim; d++)
//            {
//                for (int p=0; p<nx; p++)
//                {
//                    G(d, p, curTimeIndex) += curShapeObject->GetWeight() * curG(d,p);
//                }
//            }
//        }
//    }

//    return G;

    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    Array3D<double> G(dim,nx,T);
    G.FillArray(0.0f);

    for (int i=0; i<this->_numTargets; i++)
    {
        int curTimept = this->_targets[i]->GetTimeIndex();

        Array2D<double> g = this->_targets[i]->GradMatching(X, curTimept);

        for (int j=0; j<dim; j++)
        {
            for (int k=0; k<nx; k++)
            {
                G(j,k,curTimept) += this->_targets[i]->GetWeight()*g(j,k);
            }
        }
    }

    return G;
} 

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// KERNEL METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// RegKernelSum
//----------------------------------------------------------------
// Inputs:
//   t - time index of kernel evaluation
//   weights - 3D array (dim, num_pts, time) of weights
//
// Outputs:
//   return - 2D array (dim, num_pts) of the kernel computations 
//----------------------------------------------------------------
// Exact kernel computation W*K(this->_Xt,this->Xt) where
// W = impulse (input of function) 
//----------------------------------------------------------------
Array2D<double> Regression::RegKernelSum(int t, const Array3D<double>& weights)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    double stdV = this->_source.GetStdV();
    double sigmaV2 = pow(this->_source.GetSigmaV(),2);

    Array2D<double> output(dim,nx);
    output.FillArray(0.0f);

    //#pragma omp parallel for collapse(2)
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<nx; j++)
        {
            double argin = - ( pow(this->_X(0,i,t)-this->_X(0,j,t),2) +
                               pow(this->_X(1,i,t)-this->_X(1,j,t),2) +
                               pow(this->_X(2,i,t)-this->_X(2,j,t),2) ) / sigmaV2;

            double argout = stdV*exp(argin);

            output(0,i) = output(0,i) + argout*weights(0,j,t);
            output(1,i) = output(1,i) + argout*weights(1,j,t);
            output(2,i) = output(2,i) + argout*weights(2,j,t);
        }
    }

    return output;
}

//----------------------------------------------------------------
// GridKernelSum
//----------------------------------------------------------------
// Inputs:
//   t - time index of kernel evaluation
//   weights - 3D array (dim, num_pts, time) of weights
//
// Outputs:
//   return - 2D array (dim, num_pts) of the kernel computations 
//----------------------------------------------------------------
// Approx kernel computation W*K(this->_Xt,this->Xt) where
// W = impulse (input of function) using a grid based method
// for kernel computations
//----------------------------------------------------------------
Array2D<double> Regression::GridKernelSum(int t, const Array3D<double>& weights)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();

    Array2D<double> output(dim,nx);
    output.FillArray(0.0f);

    double* Xslice = new double[nx*dim];
    double* weightsSlice = new double[nx*dim];

    // Pull out the current slice of X and weights (in row major order)
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            int index = j*dim + i;
            Xslice[index] = this->_X(i,j,t);
            weightsSlice[index] = weights(i,j,t);
        }
    }

    int* gridlong = this->_source.GetGridAt(t)->GetLong();
    double gridpas = this->_source.GetGridAt(t)->GetPas();
    double* gridorigin = this->_source.GetGridAt(t)->GetOrigin();
    double* gridfft3kd = this->_source.GetGridAt(t)->GetFFT3KD();
    double* colfft3kd = new double[((gridlong[0]/2)+1)*gridlong[1]*gridlong[2]];

    for (int i=0; i<((gridlong[0]/2)+1); i++)
    {
        for (int j=0; j<gridlong[1]; j++)
        {
            for (int k=0; k<gridlong[2]; k++)
            {
                int indexRM = i*gridlong[2]*gridlong[1] + j*gridlong[2] + k;
                int indexCM = k*((gridlong[0]/2)+1)*gridlong[1] + j*((gridlong[0]/2)+1) + i;

                colfft3kd[indexCM] = gridfft3kd[indexRM];
            }
        }
    }

    GridOptimize gridOptim(nx, Xslice, dim, weightsSlice, nx, Xslice, gridlong, gridpas, gridorigin, colfft3kd);

    double* cmGridOut = gridOptim.GridOptim(nx, Xslice, dim, weightsSlice, nx, Xslice, gridlong, gridpas, gridorigin, colfft3kd);

    //delete gridOptim;

    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            int indexCM = j*dim+i;

            output(i,j) = cmGridOut[indexCM];
        }
    }

    //free(cmGridOut);
    //delete [] cmGridOut;
    delete [] Xslice;
    delete [] weightsSlice;
    delete [] colfft3kd;

    return output;
}

//----------------------------------------------------------------
// RegCompDeta
//----------------------------------------------------------------
// Inputs:
//   t - time index of kernel evaluation
//   eta- 3D array (dim, num_pts, time) of eta
//   momenta - 3D array (dim, num_pts, time) of momenta
//
// Outputs:
//   return - 2D array (dim, num_pts) of deta 
//----------------------------------------------------------------
// Compute derivative of auxillary variable eta using exact
// kernel computations
//----------------------------------------------------------------
Array2D<double> Regression::RegCompDeta(int t, const Array3D<double>& eta, const Array3D<double>& momenta)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    double sigmaV2 = pow(this->_source.GetSigmaV(),2);
    double gammaTimesTwo =  2*this->_source.GetGamma();

    Array2D<double> output(dim,nx);
    output.FillArray(0.0f);

    //#pragma omp parallel for collapse(2)
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<nx; j++)
        {
            // The kernel evalution in equation 4.3.17 from thesis k(x_p(t), x_q(t))
            double argin = -( pow(this->_X(0,i,t)-this->_X(0,j,t),2) +
                              pow(this->_X(1,i,t)-this->_X(1,j,t),2) +
                              pow(this->_X(2,i,t)-this->_X(2,j,t),2) ) / sigmaV2;
            double argout = this->_source.GetStdV()*exp(argin);

            // Contribution of other terms from equation 4.3.17
            // -argout is grad_1 of k(x_p(t), x_q(t))
            // First line is the dot product between momenta and eta: a_p(u) eta_q(u)
            // Second line is dot product between momenta and eta: a_q(u) eta_p(u)
            // third line is dot product between momenta: 2gamma a_p(u) a_q(u)
            argout = -argout * ( (eta(0,j,t)*momenta(0,i,t) + eta(1,j,t)*momenta(1,i,t) + eta(2,j,t)*momenta(2,i,t)) +
                                 (eta(0,i,t)*momenta(0,j,t) + eta(1,i,t)*momenta(1,j,t) + eta(2,i,t)*momenta(2,j,t)) +
                                 gammaTimesTwo * (momenta(0,i,t)*momenta(0,j,t) + momenta(1,i,t)*momenta(1,j,t) + momenta(2,i,t)*momenta(2,j,t)) );

            // Update the derivative of eta by keeping a sum and multiplied by dx/dt
            output(0,i) = output(0,i) + argout*(this->_X(0,i,t)-this->_X(0,j,t));
            output(1,i) = output(1,i) + argout*(this->_X(1,i,t)-this->_X(1,j,t));
            output(2,i) = output(2,i) + argout*(this->_X(2,i,t)-this->_X(2,j,t));
        }
    }

    return output;
}

//----------------------------------------------------------------
// RegCompDeta
//----------------------------------------------------------------
// Inputs:
//   t - time index of kernel evaluation
//   eta - 3D array (dim, num_pts, time) of eta
//   momenta - 3D array (dim, num_pts, time) of momenta
//
// Outputs:
//   return - 2D array (dim, num_pts) of deta 
//----------------------------------------------------------------
// Compute derivative of auxillary variable eta using approx
// kernel computations using a grid based method.
//----------------------------------------------------------------
Array2D<double> Regression::GridCompDeta(int t, const Array3D<double>& eta, const Array3D<double>& momenta)
{
    int dim = this->_X.GetLength();
    int nx = this->_source.GetNx();
    double gammaR = this->_source.GetGamma();

    Array2D<double> output(dim,nx);
    output.FillArray(0.0f);

    // Pull out the current slice of X (in row major order)
    double* Xslice = new double[nx*dim];
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            int index = j*dim + i;
            Xslice[index] = this->_X(i,j,t);
        }
    }

    // We concatenate several matrices together into a 24 x nx matrix
    int cDim = 8*dim;
    Array2D<double> taux(cDim,nx);
    taux.FillArray(0.0f);

    for (int j=0; j<nx; j++)
    {
        // eta1t;  eta2t;  eta3t
        taux(0,j) = eta(0,j,t);
        taux(1,j) = eta(1,j,t);
        taux(2,j) = eta(2,j,t);

        // momenta1t;  momenta2t;  momenta3t
        taux(3,j) = momenta(0,j,t);
        taux(4,j) = momenta(1,j,t);
        taux(5,j) = momenta(2,j,t);

        // eta1t.*X1t;  eta1t.*X2t;  eta1t.*X3t
        taux(6,j) = eta(0,j,t)*this->_X(0,j,t);
        taux(7,j) = eta(0,j,t)*this->_X(1,j,t);
        taux(8,j) = eta(0,j,t)*this->_X(2,j,t);

        // eta2t.*X1t;  eta2t.*X2t;  eta2t.*X3t
        taux(9, j) = eta(1,j,t)*this->_X(0,j,t);
        taux(10,j) = eta(1,j,t)*this->_X(1,j,t);
        taux(11,j) = eta(1,j,t)*this->_X(2,j,t);

        // eta3t.*X1t;  eta3t.*X2t;  eta3t.*X3t
        taux(12,j) = eta(2,j,t)*this->_X(0,j,t);
        taux(13,j) = eta(2,j,t)*this->_X(1,j,t);
        taux(14,j) = eta(2,j,t)*this->_X(2,j,t);

        // momenta1t.*X1t;  momenta1t.*X2t;  momenta1t.*X3t
        taux(15,j) = momenta(0,j,t)*this->_X(0,j,t);
        taux(16,j) = momenta(0,j,t)*this->_X(1,j,t);
        taux(17,j) = momenta(0,j,t)*this->_X(2,j,t);

        // momenta2t.*X1t;  momenta2t.*X2t;  momenta2t.*X3t
        taux(18,j) = momenta(1,j,t)*this->_X(0,j,t);
        taux(19,j) = momenta(1,j,t)*this->_X(1,j,t);
        taux(20,j) = momenta(1,j,t)*this->_X(2,j,t);

        // momenta3t.*X1t;  momenta3t.*X2t;  momenta3t.*X3t
        taux(21,j) = momenta(2,j,t)*this->_X(0,j,t);
        taux(22,j) = momenta(2,j,t)*this->_X(1,j,t);
        taux(23,j) = momenta(2,j,t)*this->_X(2,j,t);
    }

    int* gridlong = this->_source.GetGridAt(t)->GetLong();
    double gridpas = this->_source.GetGridAt(t)->GetPas();
    double* gridorigin = this->_source.GetGridAt(t)->GetOrigin();
    double* gridfft3kd = this->_source.GetGridAt(t)->GetFFT3KD();

    double* colfft3kd = new double[((gridlong[0]/2)+1)*gridlong[1]*gridlong[2]];

    for (int i=0; i<((gridlong[0]/2)+1); i++)
    {
        for (int j=0; j<gridlong[1]; j++)
        {
            for (int k=0; k<gridlong[2]; k++)
            {
                int indexRM = i*gridlong[2]*gridlong[1] + j*gridlong[2] + k;
                int indexCM = k*((gridlong[0]/2)+1)*gridlong[1] + j*((gridlong[0]/2)+1) + i;

                colfft3kd[indexCM] = gridfft3kd[indexRM];
            }
        }
    }

    double* cmtaux = new double[cDim*nx];
    for (int i=0; i<cDim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            int indexCM = j*cDim+i;
            cmtaux[indexCM] = taux(i,j);
        }
    }

    //printf("    Starting Grid Optim...\n");

    GridOptimize gridOptim(nx, Xslice, cDim, cmtaux, nx, Xslice, gridlong, gridpas, gridorigin, colfft3kd);

    double* gridOut = gridOptim.GridOptim(nx, Xslice, cDim, cmtaux, nx, Xslice, gridlong, gridpas, gridorigin, colfft3kd);

    //delete gridOptim;

    //printf("    Done with Grid Optim...\n");

    for (int i=0; i<nx; i++)
    {
        //int index = 0 + 3*i;
        int index_1 = i*cDim + 0;
        output(0,i) = output(0,i) - (momenta(0,i,t)*gridOut[index_1]*this->_X(0,i,t));
        output(1,i) = output(1,i) - (momenta(0,i,t)*gridOut[index_1]*this->_X(1,i,t));
        output(2,i) = output(2,i) - (momenta(0,i,t)*gridOut[index_1]*this->_X(2,i,t));

        //int index = 1 + 3*i;
        int index_2 = i*cDim + 1;
        output(0,i) = output(0,i) - (momenta(1,i,t)*gridOut[index_2]*this->_X(0,i,t));
        output(1,i) = output(1,i) - (momenta(1,i,t)*gridOut[index_2]*this->_X(1,i,t));
        output(2,i) = output(2,i) - (momenta(1,i,t)*gridOut[index_2]*this->_X(2,i,t));

        //int index = 2 + 3*i;
        int index_3 = i*cDim + 2;
        output(0,i) = output(0,i) - (momenta(2,i,t)*gridOut[index_3]*this->_X(0,i,t));
        output(1,i) = output(1,i) - (momenta(2,i,t)*gridOut[index_3]*this->_X(1,i,t));
        output(2,i) = output(2,i) - (momenta(2,i,t)*gridOut[index_3]*this->_X(2,i,t));

        //int index = 3 + 3*i;
        int index_4 = i*cDim + 3;
        output(0,i) = output(0,i) - (eta(0,i,t)*gridOut[index_4]*this->_X(0,i,t));
        output(1,i) = output(1,i) - (eta(0,i,t)*gridOut[index_4]*this->_X(1,i,t));
        output(2,i) = output(2,i) - (eta(0,i,t)*gridOut[index_4]*this->_X(2,i,t));

        //int index = 3 + 3*i;
        int index_5 = i*cDim + 3;
        output(0,i) = output(0,i) - (2*gammaR*gridOut[index_5]*momenta(0,i,t)*this->_X(0,i,t));
        output(1,i) = output(1,i) - (2*gammaR*gridOut[index_5]*momenta(0,i,t)*this->_X(1,i,t));
        output(2,i) = output(2,i) - (2*gammaR*gridOut[index_5]*momenta(0,i,t)*this->_X(2,i,t));

        //int index = 4 + 3*i;
        int index_6 = i*cDim + 4;
        output(0,i) = output(0,i) - (eta(1,i,t)*gridOut[index_6]*this->_X(0,i,t));
        output(1,i) = output(1,i) - (eta(1,i,t)*gridOut[index_6]*this->_X(1,i,t));
        output(2,i) = output(2,i) - (eta(1,i,t)*gridOut[index_6]*this->_X(2,i,t));

        //int index = 4 + 3*i;
        int index_7 = i*cDim + 4;
        output(0,i) = output(0,i) - (2*gammaR*gridOut[index_7]*momenta(1,i,t)*this->_X(0,i,t));
        output(1,i) = output(1,i) - (2*gammaR*gridOut[index_7]*momenta(1,i,t)*this->_X(1,i,t));
        output(2,i) = output(2,i) - (2*gammaR*gridOut[index_7]*momenta(1,i,t)*this->_X(2,i,t));

        //int index = 5 + 3*i;
        int index_8 = i*cDim + 5;
        output(0,i) = output(0,i) - (eta(2,i,t)*gridOut[index_8]*this->_X(0,i,t));
        output(1,i) = output(1,i) - (eta(2,i,t)*gridOut[index_8]*this->_X(1,i,t));
        output(2,i) = output(2,i) - (eta(2,i,t)*gridOut[index_8]*this->_X(2,i,t));

        //int index = 5 + 3*i;
        int index_9 = i*cDim + 5;
        output(0,i) = output(0,i) - (2*gammaR*gridOut[index_9]*momenta(2,i,t)*this->_X(0,i,t));
        output(1,i) = output(1,i) - (2*gammaR*gridOut[index_9]*momenta(2,i,t)*this->_X(1,i,t));
        output(2,i) = output(2,i) - (2*gammaR*gridOut[index_9]*momenta(2,i,t)*this->_X(2,i,t));

        //int index = 6 + 3*i;
        int index_10 = i*cDim + 6;
        output(0,i) = output(0,i) + (momenta(0,i,t)*gridOut[index_10]);
        output(1,i) = output(1,i) + (momenta(0,i,t)*gridOut[index_10+1]);
        output(2,i) = output(2,i) + (momenta(0,i,t)*gridOut[index_10+2]);

        //int index = 9 + 3*i;
        int index_11 = i*cDim + 9;
        output(0,i) = output(0,i) + (momenta(1,i,t)*gridOut[index_11]);
        output(1,i) = output(1,i) + (momenta(1,i,t)*gridOut[index_11+1]);
        output(2,i) = output(2,i) + (momenta(1,i,t)*gridOut[index_11+2]);

        //int index = 12 + 3*i;
        int index_12 = i*cDim + 12;
        output(0,i) = output(0,i) + (momenta(2,i,t)*gridOut[index_12]);
        output(1,i) = output(1,i) + (momenta(2,i,t)*gridOut[index_12+1]);
        output(2,i) = output(2,i) + (momenta(2,i,t)*gridOut[index_12+2]);

        //int index = 15 + 3*i;
        int index_13 = i*cDim + 15;
        output(0,i) = output(0,i) + (eta(0,i,t)*gridOut[index_13]);
        output(1,i) = output(1,i) + (eta(0,i,t)*gridOut[index_13+1]);
        output(2,i) = output(2,i) + (eta(0,i,t)*gridOut[index_13+2]);

        //int index = 15 + 3*i;
        int index_14 = i*cDim + 15;
        output(0,i) = output(0,i) + (2*gammaR*momenta(0,i,t)*gridOut[index_14]);
        output(1,i) = output(1,i) + (2*gammaR*momenta(0,i,t)*gridOut[index_14+1]);
        output(2,i) = output(2,i) + (2*gammaR*momenta(0,i,t)*gridOut[index_14+2]);

        //int index = 18 + 3*i;
        int index_15 = i*cDim + 18;
        output(0,i) = output(0,i) + (eta(1,i,t)*gridOut[index_15]);
        output(1,i) = output(1,i) + (eta(1,i,t)*gridOut[index_15+1]);
        output(2,i) = output(2,i) + (eta(1,i,t)*gridOut[index_15+2]);

        //int index = 18 + 3*i;
        int index_16 = i*cDim + 18;
        output(0,i) = output(0,i) + (2*gammaR*momenta(1,i,t)*gridOut[index_16]);
        output(1,i) = output(1,i) + (2*gammaR*momenta(1,i,t)*gridOut[index_16+1]);
        output(2,i) = output(2,i) + (2*gammaR*momenta(1,i,t)*gridOut[index_16+2]);

        //int index = 21 + 3*i;
        int index_17 = i*cDim + 21;
        output(0,i) = output(0,i) + (eta(2,i,t)*gridOut[index_17]);
        output(1,i) = output(1,i) + (eta(2,i,t)*gridOut[index_17+1]);
        output(2,i) = output(2,i) + (eta(2,i,t)*gridOut[index_17+2]);

        //int index = 18 + 3*i;
        int index_18 = i*cDim + 21;
        output(0,i) = output(0,i) + (2*gammaR*momenta(2,i,t)*gridOut[index_18]);
        output(1,i) = output(1,i) + (2*gammaR*momenta(2,i,t)*gridOut[index_18+1]);
        output(2,i) = output(2,i) + (2*gammaR*momenta(2,i,t)*gridOut[index_18+2]);
    }

    delete [] Xslice;
    delete [] cmtaux;
    delete [] colfft3kd;
    //delete [] gridOut;
    //free(gridOut);

    return output;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GRID METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetSourceGrids
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Set the source grid initially by using the min and max
// positions of the source at time 0.
//----------------------------------------------------------------
void Regression::SetSourceGrids()
{
    int nx = this->_source.GetNx();
    int T = this->_source.GetT();

    // Compute the min and max positions of the source
    double mini[3];
    mini[0] = 1e10f;   mini[1] = 1e10f;   mini[2] = 1e10f;
    double maxi[3];
    maxi[0] = -1e10f;  maxi[1] = -1e10f;  maxi[2] = -1e10f;
    for (int i=0; i<nx; i++)
    {
        double curX = this->_X(0,i,0);
        double curY = this->_X(1,i,0);
        double curZ = this->_X(2,i,0);

        if (curX < mini[0])
            mini[0] = curX;
        if (curY < mini[1])
            mini[1] = curY;
        if (curZ < mini[2])
            mini[2] = curZ;

        if (curX > maxi[0])
            maxi[0] = curX;
        if (curY > maxi[1])
            maxi[1] = curY;
        if (curZ > maxi[2])
            maxi[2] = curZ;
    }

    //printf("minix %0.4f  miniy %0.4f  miniz %0.4f\n", mini[0], mini[1], mini[2]);
    //printf("maxix %0.4f  maxiy %0.4f  maxiz %0.4f\n", maxi[0], maxi[1], maxi[2]);

    for (int t=0; t<T; t++)
    {
        this->_source.GetGridAt(t)->UpdateGrid(mini, maxi, this->_source.GetSigmaV());
    }
}

//----------------------------------------------------------------
// UpdateSourceGrids
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Update the source grids by using the min and max positions of 
// the source at each time t.
//----------------------------------------------------------------
void Regression::UpdateSourceGrids()
{
    int T = this->_source.GetT();
    int nx = this->_source.GetNx();

    for (int t=0; t<T; t++)
    {
        // Compute the min and max positions of the source
        double mini[3];
        mini[0] = 1e10f;   mini[1] = 1e10f;   mini[2] = 1e10f;
        double maxi[3];
        maxi[0] = -1e10f;  maxi[1] = -1e10f;  maxi[2] = -1e10f;
        for (int i=0; i<nx; i++)
        {
            double curX = this->_X(0,i,t);
            double curY = this->_X(1,i,t);
            double curZ = this->_X(2,i,t);

            if (curX < mini[0])
                mini[0] = curX;
            if (curY < mini[1])
                mini[1] = curY;
            if (curZ < mini[2])
                mini[2] = curZ;

            if (curX > maxi[0])
                maxi[0] = curX;
            if (curY > maxi[1])
                maxi[1] = curY;
            if (curZ > maxi[2])
                maxi[2] = curZ;
        }

        this->_source.GetGridAt(t)->UpdateGrid(mini, maxi, this->_source.GetSigmaV());
    }
}

//----------------------------------------------------------------
// SetTargetGrids
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Set the target grids by using the min and max positions of 
// source and target
//----------------------------------------------------------------
void Regression::SetTargetGrids()
{
    int numTargets = this->_numTargets;

    //printf("Num targets = %d\n", numTargets);

    for (int t=0; t<numTargets; t++)
    {
        //printf("On target %d\n", t+1);

        // Compute the min and max positions of the source
        double mini[3];
        mini[0] = 1e10f;   mini[1] = 1e10f;   mini[2] = 1e10f;
        double maxi[3];
        maxi[0] = -1e10f;  maxi[1] = -1e10f;  maxi[2] = -1e10f;

        Array2D<int> sourceFaces = this->_targets[t]->GetVx();

        //printf("Num faces = %d\n", sourceFaces.GetWidth());

        for (int i=0; i<this->_targets[t]->GetNvx(); i++)
        {
            double curX = this->_X(0,sourceFaces(0,i),0);
            double curY = this->_X(1,sourceFaces(1,i),0);
            double curZ = this->_X(2,sourceFaces(2,i),0);

            if (curX < mini[0])
                mini[0] = curX;
            if (curY < mini[1])
                mini[1] = curY;
            if (curZ < mini[2])
                mini[2] = curZ;

            if (curX > maxi[0])
                maxi[0] = curX;
            if (curY > maxi[1])
                maxi[1] = curY;
            if (curZ > maxi[2])
                maxi[2] = curZ;
        }

        //printf("Done with source\n");

        Array2D<double> targetPts = this->_targets[t]->GetY();

        //printf("Num points = %d\n", targetPts.GetWidth());

        for (int i=0; i<this->_targets[t]->GetNy(); i++)
        {
            double curX = targetPts(0,i);
            double curY = targetPts(1,i);
            double curZ = targetPts(2,i);

            if (curX < mini[0])
                mini[0] = curX;
            if (curY < mini[1])
                mini[1] = curY;
            if (curZ < mini[2])
                mini[2] = curZ;

            if (curX > maxi[0])
                maxi[0] = curX;
            if (curY > maxi[1])
                maxi[1] = curY;
            if (curZ > maxi[2])
                maxi[2] = curZ;
        }

        //printf("Done with target\n");

        //printf("minix %0.4f  miniy %0.4f  miniz %0.4f\n", mini[0], mini[1], mini[2]);
        //printf("maxix %0.4f  maxiy %0.4f  maxiz %0.4f\n", maxi[0], maxi[1], maxi[2]);

        this->_targets[t]->GetGrid()->UpdateGrid(mini, maxi, this->_targets[t]->GetSigmaW());

        //printf("Done Updating grid...\n");
    }
}

//----------------------------------------------------------------
// UpdateTargetGrids
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Updates the target grids by call SetTargetGrids()
//----------------------------------------------------------------
void Regression::UpdateTargetGrids()
{
    int numTargets = this->_numTargets;

    for (int t=0; t<numTargets; t++)
    {
        // Compute the min and max positions of the source
        double mini[3];
        mini[0] = 1e10f;   mini[1] = 1e10f;   mini[2] = 1e10f;
        double maxi[3];
        maxi[0] = -1e10f;  maxi[1] = -1e10f;  maxi[2] = -1e10f;

        Array2D<int> sourceFaces = this->_targets[t]->GetVx();

        for (int i=0; i<this->_targets[t]->GetNvx(); i++)
        {
            double curX = this->_X(0,sourceFaces(0,i),this->_targets[t]->GetTimeIndex());
            double curY = this->_X(1,sourceFaces(1,i),this->_targets[t]->GetTimeIndex());
            double curZ = this->_X(2,sourceFaces(2,i),this->_targets[t]->GetTimeIndex());

            if (curX < mini[0])
                mini[0] = curX;
            if (curY < mini[1])
                mini[1] = curY;
            if (curZ < mini[2])
                mini[2] = curZ;

            if (curX > maxi[0])
                maxi[0] = curX;
            if (curY > maxi[1])
                maxi[1] = curY;
            if (curZ > maxi[2])
                maxi[2] = curZ;
        }

        Array2D<double> targetPts = this->_targets[t]->GetY();
        for (int i=0; i<this->_targets[t]->GetNy(); i++)
        {
            double curX = targetPts(0,i);
            double curY = targetPts(1,i);
            double curZ = targetPts(2,i);

            if (curX < mini[0])
                mini[0] = curX;
            if (curY < mini[1])
                mini[1] = curY;
            if (curZ < mini[2])
                mini[2] = curZ;

            if (curX > maxi[0])
                maxi[0] = curX;
            if (curY > maxi[1])
                maxi[1] = curY;
            if (curZ > maxi[2])
                maxi[2] = curZ;
        }

        this->_targets[t]->GetGrid()->UpdateGrid(mini, maxi, this->_targets[t]->GetSigmaW());
    }
}
