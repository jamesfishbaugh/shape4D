//--------------------------------------------------------------------------
// EXRegressionParams:  Shape regression params
//--------------------------------------------------------------------------

#include <string.h>			// For memcpy
#include <cstdio>

#include "regressionparams.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXRegressionParams
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
RegressionParams::RegressionParams()
{
    this->_nx = 0;
    this->_sigmaV = 0.0f;
    this->_gamma = 0.0f;
    this->_T = 0;
    this->_stdV = 0.0f;
    this->_v0Weight = 0.0f;
    this->_useInitV0 = false;
    this->_estimateBaseline = true;
    this->_kernelType = (char*) "exact";
    this->_maxIters = 500;
    this->_breakRatio = 1e-6;
    this->_useFista = true;

    this->_grids = new Grid*[this->_T];
    for (int t=0; t<this->_T; t++)
    {
        this->_grids[t] = NULL;
    }
}

//----------------------------------------------------------------
// EXRegressionParams
//----------------------------------------------------------------
// Inputs:
//   x - 2D array (dim, num_pts) of source points
//   sigmaV - scale of deformation
//   gamma - weight of regularity
//   xtime - 1D array (time) of time discretization
//
// Outputs:
//----------------------------------------------------------------
// Init constructor - recommened constructor
//----------------------------------------------------------------
RegressionParams::RegressionParams(const Array2D<double>& x, double sigmaV, double gamma, const Array1D<double>& xtime)
{
    this->_nx = x.GetWidth();
    this->_x = x;
    this->_sigmaV = sigmaV;
    this->_gamma = gamma;
    this->_T = xtime.GetLength();
    this->_xtime = xtime;
    this->_stdV = 1.0f;
    this->_v0Weight = 1.0f;
    this->_useInitV0 = false;
    this->_estimateBaseline = true;
    this->_kernelType = (char*) "exact";
    this->_maxIters = 500;
    this->_breakRatio = 1e-6;
    this->_useFista = true;

    this->_grids = new Grid*[this->_T];
    for (int t=0; t<this->_T; t++)
    {
        this->_grids[t] = new Grid();
    }

}

//----------------------------------------------------------------
// EXRegressionParams
//----------------------------------------------------------------
// Inputs:
//   x - 2D array (dim, num_pts) of source points
//   sigmaV - scale of deformation
//   gamma - weight of regularity
//   xtime - 1D array (time) of time discretization
//   initV0 - 2D array (dim, num_pts) of inital velocity
//
// Outputs:
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
RegressionParams::RegressionParams(const Array2D<double>& x, double sigmaV, double gamma, const Array1D<double>& xtime, const Array2D<double>& initV0)
{
    this->_nx = x.GetWidth();
    this->_x = x;
    this->_sigmaV = sigmaV;
    this->_gamma = gamma;
    this->_T = xtime.GetLength();
    this->_xtime = xtime;
    this->_stdV = 1.0f;
    this->_v0Weight = 0.0f;
    this->_useInitV0 = false;
    this->_initV0 = initV0;
    this->_estimateBaseline = true;
    this->_kernelType = (char*) "exact";
    this->_maxIters = 500;
    this->_breakRatio = 1e-6;
    this->_useFista = true;
    this->_baselineSmoothing = 1.0;

    this->_grids = new Grid*[this->_T];
    for (int t=0; t<this->_T; t++)
    {
        this->_grids[t] = new Grid();
    }

}

//----------------------------------------------------------------
// EXRegressionParams
//----------------------------------------------------------------
// Inputs: 
//   source - EXRegressionParams object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
RegressionParams::RegressionParams(const RegressionParams& source)
{
    printf("In copy constructor\n");
    this->_nx = source._nx;
    this->_x = source._x;
    this->_sigmaV = source._sigmaV;
    this->_gamma = source._gamma;
    this->_T = source._T;
    this->_xtime = source._xtime;
    this->_stdV = source._stdV;
    this->_initV0 = source._initV0;
    this->_useInitV0 = source._useInitV0;
    this->_v0Weight = source._v0Weight;
    this->_estimateBaseline = source._estimateBaseline;
    this->_kernelType = source._kernelType;
    this->_maxIters = source._maxIters;
    this->_breakRatio = source._breakRatio;
    this->_useFista = source._useFista;
    this->_baselineSmoothing = 1.0;

    this->_grids = new Grid*[this->_T];
    for (int t=0; t<this->_T; t++)
    {
        this->_grids[t] = source._grids[t];
    }
}

//----------------------------------------------------------------
// ~EXRegressionParams
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
RegressionParams::~RegressionParams()
{
    //for (int t=0; t<this->_T; t++)
    //{
    //	delete [] this->_grids[t];
    //}

    delete [] this->_grids;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetX
//----------------------------------------------------------------
// Inputs: 
//   x - 2D array (dim, num_pts) of source points
//
// Outputs:
//----------------------------------------------------------------
// Setter for source points
//----------------------------------------------------------------
void RegressionParams::SetX(const Array2D<double>& x)
{
    this->_nx = x.GetWidth();
    this->_x = x;
}

//----------------------------------------------------------------
// SetSigmaV
//----------------------------------------------------------------
// Inputs: 
//   sigmaV - scale of deformation
//
// Outputs:
//----------------------------------------------------------------
// Setter for scale of deformation
//----------------------------------------------------------------
void RegressionParams::SetSigmaV(double sigmaV)
{
    this->_sigmaV = sigmaV;
}

//----------------------------------------------------------------
// SetGamma
//----------------------------------------------------------------
// Inputs: 
//   gamma - weight for regularity
//
// Outputs:
//----------------------------------------------------------------
// Setter for weight for regularity
//----------------------------------------------------------------
void RegressionParams::SetGamma(double gamma)
{
    this->_gamma = gamma;
}

//----------------------------------------------------------------
// SetStdV
//----------------------------------------------------------------
// Inputs: 
//   stdV - stdV
//
// Outputs:
//----------------------------------------------------------------
// Setter for stdV
//----------------------------------------------------------------
void RegressionParams::SetStdV(double stdV)
{
    this->_stdV = stdV;
}

//----------------------------------------------------------------
// SetXTime
//----------------------------------------------------------------
// Inputs: 
//   xtime - 1D array (time) of time discretization
//
// Outputs:
//----------------------------------------------------------------
// Setter for time discretization array
//----------------------------------------------------------------
void RegressionParams::SetXTime(const Array1D<double>& xtime)
{
    this->_T = xtime.GetLength();
    this->_xtime = xtime;
}

//----------------------------------------------------------------
// SetInitV0
//----------------------------------------------------------------
// Inputs: 
//   initV0 - 2D array (dim, num_pts) of initial velocity
//
// Outputs:
//----------------------------------------------------------------
// Setter for initial velocity
//----------------------------------------------------------------
void RegressionParams::SetInitV0(const Array2D<double>& initV0)
{
    this->_useInitV0 = true;
    this->_initV0 = initV0;
    this->_v0Weight = 0.0f;
}

//----------------------------------------------------------------
// SetV0Weight
//----------------------------------------------------------------
// Inputs: 
//   v0Weight - importance of v0 in regression
//
// Outputs:
//----------------------------------------------------------------
// Setter for v0 weight
//----------------------------------------------------------------
void RegressionParams::SetV0Weight(double v0Weight)
{
    this->_v0Weight = v0Weight;
}

//----------------------------------------------------------------
// SetShouldUseInitV0
//----------------------------------------------------------------
// Inputs: 
//   yesNo - use only init V0
//
// Outputs:
//----------------------------------------------------------------
// Setter for ShouldUseInitV0
//----------------------------------------------------------------
void RegressionParams::SetShouldUseInitV0(bool yesNo)
{
    this->_useInitV0 = yesNo;
}

//----------------------------------------------------------------
// SetShouldEstimateBaseline
//----------------------------------------------------------------
// Inputs:
//   yesNo - estimate baseline or not
//
// Outputs:
//----------------------------------------------------------------
// Setter for EstimateBaseline
//----------------------------------------------------------------
void RegressionParams::SetShouldEstimateBaseline(bool yesNo)
{
    this->_estimateBaseline = yesNo;
}

//----------------------------------------------------------------
// SetMaxIters
//----------------------------------------------------------------
// Inputs:
//   maxIters - the maximum number of iterations
//
// Outputs:
//----------------------------------------------------------------
// Setter for maxIters
//----------------------------------------------------------------
void RegressionParams::SetMaxIters(int maxIters)
{
    this->_maxIters = maxIters;
}

//----------------------------------------------------------------
// SetBreakRatio
//----------------------------------------------------------------
// Inputs:
//   breakRatio - the break ratio for convergence
//
// Outputs:
//----------------------------------------------------------------
// Setter for breakRatio
//----------------------------------------------------------------
void RegressionParams::SetBreakRatio(double breakRatio)
{
    this->_breakRatio = breakRatio;
}

//----------------------------------------------------------------
// SetKernelType
//----------------------------------------------------------------
// Inputs:
//   kernelType - "exact" or "p3m" kernel type
//
// Outputs:
//----------------------------------------------------------------
// Setter for kernel type
//----------------------------------------------------------------
void RegressionParams::SetKernelType(char* kernelType)
{
    this->_kernelType=  kernelType;

    // If it not a supported kernel type set it to exact
    if ((strcmp(this->_kernelType,"exact") != 0) && (strcmp(this->_kernelType,"p3m") != 0))
    {
        this->_kernelType = (char*) "exact";
    }
}

//----------------------------------------------------------------
// SetShouldUseFista
//----------------------------------------------------------------
// Inputs:
//   yesNo - use use FISTA algorithm
//
// Outputs:
//----------------------------------------------------------------
// Setter for this->_useFista
//----------------------------------------------------------------
void RegressionParams::SetShouldUseFista(bool yesNo)
{
    this->_useFista = yesNo;
}

//----------------------------------------------------------------
// SetBaselineSmoothing
//----------------------------------------------------------------
// Inputs:
//   baselineSmoothing - factor for smoothing baseline shape
//
// Outputs:
//----------------------------------------------------------------
// Setter for this->_baselineSmoothing
//----------------------------------------------------------------
void RegressionParams::SetBaselineSmoothing(double baselineSmoothing)
{
    this->_baselineSmoothing = baselineSmoothing;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetX
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 2D array (dim, num_pts) of source points
//----------------------------------------------------------------
// Returns source points
//----------------------------------------------------------------
const Array2D<double> RegressionParams::GetX() const
{
    return this->_x;
}

//----------------------------------------------------------------
// GetNx
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - number of source points
//----------------------------------------------------------------
// Returns the number of source points
//----------------------------------------------------------------
int RegressionParams::GetNx() const
{
    return this->_nx;
}

//----------------------------------------------------------------
// GetSigmaV
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - scale of deformation
//----------------------------------------------------------------
// Returns the scale of deformation
//----------------------------------------------------------------
double RegressionParams::GetSigmaV() const
{
    return this->_sigmaV;
}

//----------------------------------------------------------------
// GetGamma
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - weight of regularity
//----------------------------------------------------------------
// Returns the weight of regularity
//----------------------------------------------------------------
double RegressionParams::GetGamma() const
{
    return this->_gamma;
}

//----------------------------------------------------------------
// GetStdV
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - stdV
//----------------------------------------------------------------
// Returns stdV
//----------------------------------------------------------------
double RegressionParams::GetStdV() const
{
    return this->_stdV;
}

//----------------------------------------------------------------
// GetT
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - number of time points in the discretization
//----------------------------------------------------------------
// Returns the number of time points in the discretization
//----------------------------------------------------------------
int RegressionParams::GetT() const
{
    return this->_T;
}

//----------------------------------------------------------------
// GetXTime
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 1D array (time) of time discretization
//----------------------------------------------------------------
// Returns the time discretization array
//----------------------------------------------------------------
const Array1D<double> RegressionParams::GetXTime() const
{
    return this->_xtime;
}

//----------------------------------------------------------------
// ShouldUseInitV0
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - true if using an initial velocity, false otherwise
//----------------------------------------------------------------
// Returns true if using an initial velocity, false otherwise
//----------------------------------------------------------------
bool RegressionParams::ShouldUseInitV0() const
{
    return this->_useInitV0;
}

//----------------------------------------------------------------
// ShouldEstimateBaseline
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - true if using an initial velocity, false otherwise
//----------------------------------------------------------------
// Returns true if estimating baseline, false otherwise
//----------------------------------------------------------------
bool RegressionParams::ShouldEstimateBaseline() const
{
    return this->_estimateBaseline;
}

//----------------------------------------------------------------
// GetKernelType
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - string representing the kernel type
//----------------------------------------------------------------
// Returns "exact" or "p3m" kernel type
//----------------------------------------------------------------
const char* RegressionParams::GetKernelType() const
{
    return this->_kernelType;
}

//----------------------------------------------------------------
// GetMaxIters
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - maximum number of iterations
//----------------------------------------------------------------
// Returns maximum number of iterations for the optimzier
//----------------------------------------------------------------
int RegressionParams::GetMaxIters() const
{
    return this->_maxIters;
}

//----------------------------------------------------------------
// GetBreakRatio
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - break ratio for convergence
//----------------------------------------------------------------
// Returns break ratio for the optimzier
//----------------------------------------------------------------
double RegressionParams::GetBreakRatio() const
{
    return this->_breakRatio;
}

//----------------------------------------------------------------
// GetInitV0
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 2D array (dim, num_pts) of initial velocity
//----------------------------------------------------------------
// Returns initial velocity
//----------------------------------------------------------------
const Array2D<double> RegressionParams::GetInitV0() const
{	
    return this->_initV0;
}

//----------------------------------------------------------------
// GetInitV0
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - 2D array (dim, num_pts) of initial velocity
//----------------------------------------------------------------
// Returns initial velocity
//----------------------------------------------------------------
double RegressionParams::GetV0Weight() const
{
    return this->_v0Weight;
}

//----------------------------------------------------------------
// ShouldUseFista
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - true if using FISTA, false otherwise
//----------------------------------------------------------------
// Returns true if using FISTA, false otherwise
//----------------------------------------------------------------
bool RegressionParams::ShouldUseFista() const
{
    return this->_useFista;
}

//----------------------------------------------------------------
// GetBaselineSmoothing
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - baseline shape smoothing factor
//----------------------------------------------------------------
// Returns the baseline shape smoothing factor
//----------------------------------------------------------------
double RegressionParams::GetBaselineSmoothing() const
{
    return this->_baselineSmoothing;
}

//----------------------------------------------------------------
// GetGridAt
//----------------------------------------------------------------
// Inputs:
//   index - index of desired grid
//
// Outputs:
//   return - pointer to a grid
//----------------------------------------------------------------
// Returns a pointer to the grid at index
//----------------------------------------------------------------
Grid* RegressionParams::GetGridAt(int index)
{
    return this->_grids[index];
}

//----------------------------------------------------------------
// GetGridAt
//----------------------------------------------------------------
// Inputs:
//   index - index of desired grid
//
// Outputs:
//   return - pointer to a grid (const)
//----------------------------------------------------------------
// Returns a pointer (const) to the grid at index
//----------------------------------------------------------------
const Grid* RegressionParams::GetGridAt(int index) const
{
    return this->_grids[index];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// OVERLOADED OPERATORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// operator = 
//----------------------------------------------------------------
// Inputs: 
//   source - EXRegressionParams object to copy
//
// Outputs:
//   return - copied EXRegressionParams object 
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
RegressionParams& RegressionParams::operator = (const RegressionParams& source)
{
    if (this != &source)
    {
        this->_nx = source._nx;
        this->_x = source._x;
        this->_sigmaV = source._sigmaV;
        this->_gamma = source._gamma;
        this->_T = source._T;
        this->_xtime = source._xtime;
        this->_stdV = source._stdV;
        this->_initV0 = source._initV0;
        this->_useInitV0 = source._useInitV0;
        this->_v0Weight = source._v0Weight;
        this->_estimateBaseline = source._estimateBaseline;
        this->_kernelType = source._kernelType;
        this->_maxIters = source._maxIters;
        this->_breakRatio = source._breakRatio;
        this->_useFista = source._useFista;
        this->_baselineSmoothing = source._baselineSmoothing;

        if (this->_grids != 0)
        {
            delete [] this->_grids;
        }
        this->_grids = new Grid*[this->_T];
        for (int t=0; t<this->_T; t++)
        {
            this->_grids[t] = source._grids[t];
        }
    }

    return *this;
}
