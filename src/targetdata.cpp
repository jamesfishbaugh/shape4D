//--------------------------------------------------------------------------
// EXTargetDataData:  Virtual base class representing target data for 
// regression
//--------------------------------------------------------------------------

#include <string.h>			// For memcpy

#include "targetdata.h"
#include "helper.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// DESTRUCTOR
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ~EXTargetData
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
TargetData::~TargetData()
{
    delete this->_grid;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetY
//----------------------------------------------------------------
// Inputs:
//   y - a 2D array (dim, num_pts) of target point
//
// Outputs:
//----------------------------------------------------------------
// Setter for target points
//----------------------------------------------------------------
void TargetData::SetY(const Array2D<double>& y)
{
    this->_ny = y.GetWidth();
    this->_y = y;
}

//----------------------------------------------------------------
// SetVY
//----------------------------------------------------------------
// Inputs:
//   vy - a 2D array (dim, num_pts) of target triangles
//
// Outputs:
//----------------------------------------------------------------
// Setter for target triangles
//----------------------------------------------------------------
void TargetData::SetVy(const Array2D<int>& vy)
{
    this->_nvy = vy.GetWidth();
    this->_vy = vy;
}

//----------------------------------------------------------------
// SetVx
//----------------------------------------------------------------
// Inputs:
//   vx - a 2D array (dim, num_pts) of source triangles that
//        correspond to this target
//
// Outputs:
//----------------------------------------------------------------
// Setter for source triangles corresponding to this target
//----------------------------------------------------------------
void TargetData::SetVx(const Array2D<int>& vx)
{
    this->_nvx = vx.GetWidth();
    this->_vx = vx;
}

//----------------------------------------------------------------
// SetSigmaW
//----------------------------------------------------------------
// Inputs:
//   sigmaW - the scale for the norm on currents
//
// Outputs:
//----------------------------------------------------------------
// Setter for the scale for the norm on currents
//----------------------------------------------------------------
void TargetData::SetSigmaW(double sigmaW)
{
    this->_sigmaW = sigmaW;
}

//----------------------------------------------------------------
// SetTimept
//----------------------------------------------------------------
// Inputs:
//   timept - the time corresponding to this target in whatever
//            time units used in the experiment
//
// Outputs:
//----------------------------------------------------------------
// Setter for the time corresponding to this target
//----------------------------------------------------------------
void TargetData::SetTimept(double timept)
{
    this->_timept = timept;
}

//----------------------------------------------------------------
// SetTimeIndex
//----------------------------------------------------------------
// Inputs:
//   timeIndex - the index into the time discretization array 
//               
// Outputs:
//----------------------------------------------------------------
// Setter for the time index corresponding to this target
//----------------------------------------------------------------
void TargetData::SetTimeIndex(int timeIndex)
{
    this->_timeIndex = timeIndex;
}

//----------------------------------------------------------------
// SetWeight
//----------------------------------------------------------------
// Inputs:
//   weight - the importance of this target for regression               
//
// Outputs:
//----------------------------------------------------------------
// Setter for the weight corresponding to this target
//----------------------------------------------------------------
void TargetData::SetWeight(double weight)
{
    this->_weight = weight;
}

//----------------------------------------------------------------
// SetKernelType
//----------------------------------------------------------------
// Inputs:
//   kernelType - kernel type: p3m or exact
//
// Outputs:
//----------------------------------------------------------------
// Setter for the kernel type for data matching
//----------------------------------------------------------------
void TargetData::SetKernelType(char* kernelType)
{
    this->_kernelType=  kernelType;

    // If it not a supported kernel type set it to exact
    if ((strcmp(this->_kernelType,"exact") != 0) && (strcmp(this->_kernelType,"p3m") != 0))
    {
        this->_kernelType = (char*) "exact";
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetNy
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the number of target points
//----------------------------------------------------------------
// Returns the number of target points
//----------------------------------------------------------------
int TargetData::GetNy() const
{
    return this->_ny;
}

//----------------------------------------------------------------
// GetY
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - an array (dim, num_pts) of the target points
//----------------------------------------------------------------
// Returns the target points
//----------------------------------------------------------------
const Array2D<double> TargetData::GetY() const
{
    return this->_y;
}

//----------------------------------------------------------------
// GetNvy
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the number of target triangles
//----------------------------------------------------------------
// Returns the number of target triangles
//----------------------------------------------------------------
int TargetData::GetNvy() const
{
    return this->_nvy;
}

//----------------------------------------------------------------
// GetVy
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - an array (dim, num_pts) of the target triangles
//----------------------------------------------------------------
// Returns the target triangles
//----------------------------------------------------------------
const Array2D<int> TargetData::GetVy() const
{
    return this->_vy;
}

//----------------------------------------------------------------
// GetNvx
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the number of source triangles corresponding to
//            this target
//----------------------------------------------------------------
// Returns the number of source triangles corresponding to this
// target
//----------------------------------------------------------------
int TargetData::GetNvx() const
{
    return this->_nvx;
}

//----------------------------------------------------------------
// GetVx
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - an array (dim, num_pts) of the source triangles
//            corresponding to this target
//----------------------------------------------------------------
// Returns the source triangles corresponding to this target
//----------------------------------------------------------------
const Array2D<int> TargetData::GetVx() const
{
    return this->_vx;
}

//----------------------------------------------------------------
// GetSigma
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the scale of the norm on currents
//----------------------------------------------------------------
// Returns the scale of the norm on currents
//----------------------------------------------------------------
double TargetData::GetSigmaW() const
{
    return this->_sigmaW;
}

//----------------------------------------------------------------
// GetTimept
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the time corresponding to this target in whatever
//            time units used in the experiment
//----------------------------------------------------------------
// Returns the time corresponding to this target in whatever
// time units used in the experiment
//----------------------------------------------------------------
double TargetData::GetTimept() const
{
    return this->_timept;
}

//----------------------------------------------------------------
// GetTimeIndex
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the index into the time discretization array 
//----------------------------------------------------------------
// Returns the index into the time discretization array 
//----------------------------------------------------------------
double TargetData::GetTimeIndex() const
{
    return this->_timeIndex;
}

//----------------------------------------------------------------
// GetWeight
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the importance of this target in the regression
//----------------------------------------------------------------
// Returns the importance of this target in the regression
//----------------------------------------------------------------
double TargetData::GetWeight() const
{
    return this->_weight;
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
const char* TargetData::GetKernelType() const
{
    return this->_kernelType;
}


//----------------------------------------------------------------
// GetGrid
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - the grid associated with this target which is used
//            for fast kernel computations
//----------------------------------------------------------------
// Returns the grid associated with this target
//----------------------------------------------------------------
Grid* TargetData::GetGrid()
{
    return this->_grid;
}

//----------------------------------------------------------------
// GetGrid
//----------------------------------------------------------------
// Inputs:
//               
// Outputs:
//   return - const grid associated with this target which is used
//            for fast kernel computations
//----------------------------------------------------------------
// Returns the grid (const) associated with this target
//----------------------------------------------------------------
const Grid* TargetData::GetGrid() const
{
    return this->_grid;
}
