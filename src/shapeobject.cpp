//--------------------------------------------------------------------------
// ShapeObject:  Virtual base class representing a shape object
//--------------------------------------------------------------------------

#include <string.h>			// For memcpy

#include "shapeobject.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTOR/DESTRUCTOR
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ShapeObject
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
ShapeObject::ShapeObject()
{
    this->_numPoints = 0;
    this->_numEdges = 0;
    this->_sigmaW = 0.0;
    this->_timept = 0.0;
    this->_timeIndex = 0;
    this->_weight = 1.0;
}

//----------------------------------------------------------------
// tmpLandmarks
//----------------------------------------------------------------
// Inputs:
//   shape - the shape object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
ShapeObject::ShapeObject(const ShapeObject& shape)
{
    this->_points = shape._points;
    this->_numPoints = shape._numPoints;
    this->_edges = shape._edges;
    this->_numEdges = shape._numEdges;
    this->_sigmaW = shape._sigmaW;
    this->_weight = shape._weight;
    this->_timept = shape._timept;
    this->_timeIndex = shape._timeIndex;
}

//----------------------------------------------------------------
// ~ShapeObject
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
ShapeObject::~ShapeObject()
{
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetPoints
//----------------------------------------------------------------
// Inputs:
//   shapePoints - a 2D array (dim, num_pts) of shape points
//
// Outputs:
//----------------------------------------------------------------
// Setter for shape points
//----------------------------------------------------------------
void ShapeObject::SetPoints(const Array2D<double>& shapePoints)
{
    this->_numPoints = shapePoints.GetWidth();
    this->_points = shapePoints;
}

//----------------------------------------------------------------
// SetEdges
//----------------------------------------------------------------
// Inputs:
//   edges - a 2D array (dim, num_pts) of edge connectivity
//
// Outputs:
//----------------------------------------------------------------
// Setter for shape edges
//----------------------------------------------------------------
void ShapeObject::SetEdges(const Array2D<int>& edges)
{
    this->_numEdges = edges.GetWidth();
    this->_edges = edges;
}

//----------------------------------------------------------------
// SetSigmaW
//----------------------------------------------------------------
// Inputs:
//   sigmaW - the scale for the kernel for shape matching
//
// Outputs:
//----------------------------------------------------------------
// Setter for the scale for the kernel for shape matching
//----------------------------------------------------------------
void ShapeObject::SetSigmaW(double sigmaW)
{
    this->_sigmaW = sigmaW;
}

//----------------------------------------------------------------
// SetTimept
//----------------------------------------------------------------
// Inputs:
//   timept - the time corresponding to this shape object
//
// Outputs:
//----------------------------------------------------------------
// Setter for the time corresponding to this shape
//----------------------------------------------------------------
void ShapeObject::SetTimept(double timept)
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
// Setter for the time index corresponding to this shape
//----------------------------------------------------------------
void ShapeObject::SetTimeIndex(int timeIndex)
{
    this->_timeIndex = timeIndex;
}

//----------------------------------------------------------------
// SetWeight
//----------------------------------------------------------------
// Inputs:
//   weight - the importance of this shape
//
// Outputs:
//----------------------------------------------------------------
// Setter for the weight corresponding to this shape
//----------------------------------------------------------------
void ShapeObject::SetWeight(double weight)
{
    this->_weight = weight;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetNumberOfPoints
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - the number of shape points
//----------------------------------------------------------------
// Returns the number of shape points
//----------------------------------------------------------------
int ShapeObject::GetNumberOfPoints() const
{
    return this->_numPoints;
}

//----------------------------------------------------------------
// GetPoints
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - an array (dim, num_pts) of the shape points
//----------------------------------------------------------------
// Returns the shape points
//----------------------------------------------------------------
const Array2D<double> ShapeObject::GetPoints() const
{
    return this->_points;
}

//----------------------------------------------------------------
// GetNumberOfEdges
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - the number of shape edges
//----------------------------------------------------------------
// Returns the number of shape edges
//----------------------------------------------------------------
int ShapeObject::GetNumberOfEdges() const
{
    return this->_numEdges;
}

//----------------------------------------------------------------
// GetEdges
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - an array (dim, num_pts) of the shape edges
//----------------------------------------------------------------
// Returns the shape edges
//----------------------------------------------------------------
const Array2D<int> ShapeObject::GetEdges() const
{
    return this->_edges;
}

//----------------------------------------------------------------
// GetSigma
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - the scale of the kernel for shape matching
//----------------------------------------------------------------
// Returns the scale of the kernel for shape matching
//----------------------------------------------------------------
double ShapeObject::GetSigmaW() const
{
    return this->_sigmaW;
}

//----------------------------------------------------------------
// GetTimept
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - the time corresponding to this shape
//----------------------------------------------------------------
// Returns the time corresponding to this shape
//----------------------------------------------------------------
double ShapeObject::GetTimept() const
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
double ShapeObject::GetTimeIndex() const
{
    return this->_timeIndex;
}

//----------------------------------------------------------------
// GetWeight
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - the weight for this shape
//----------------------------------------------------------------
// Returns the weight for this shape
//----------------------------------------------------------------
double ShapeObject::GetWeight() const
{
    return this->_weight;
}

