//----------------------------------------------------------------
// SurfaceCurrent:  Shape data represented as a surface current
//----------------------------------------------------------------

#include "tmpsurfacecurrent.h"

#include <stdio.h>			// For printf

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS/DESTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// tmpSurfaceCurrent
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
tmpSurfaceCurrent::tmpSurfaceCurrent()
{
}

//----------------------------------------------------------------
// tmpSurfaceCurrent
//----------------------------------------------------------------
// Inputs:
//   shape - tmpSurfaceCurrent object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
tmpSurfaceCurrent::tmpSurfaceCurrent(const tmpSurfaceCurrent& shape)
{
    this->_numPoints = shape._numPoints;
    this->_points = shape._points;
    this->_numEdges = shape._numEdges;
    this->_edges = shape._edges;
    this->_sigmaW = shape._sigmaW;
    this->_timept = shape._timept;
    this->_timeIndex = shape._timeIndex;
}

//----------------------------------------------------------------
// ~tmpSurfaceCurrent
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
tmpSurfaceCurrent::~tmpSurfaceCurrent()
{
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// INTERFACE TO COMPUTE DATA MATCHING METRIC AND GRADIENT
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// Matching
//----------------------------------------------------------------
// Inputs:
//   x - a 2D array (dim, num_pts) of points to compare
//
// Outputs:
//----------------------------------------------------------------
// Computes the data matching metric between sets of surface
// currents
//----------------------------------------------------------------
double tmpSurfaceCurrent::Matching(const ShapeObject* shape)
{
    Array2D<double> x = shape->GetPoints();
    int dim = x.GetLength();
    int nx = x.GetWidth();
    double sum = 0.0f;

    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            double curpt = x(i,j);
            double diff = (curpt - this->_points(i,j));

            sum += (diff*diff);
        }
    }

    return sum;
}

//----------------------------------------------------------------
// GradMatching
//----------------------------------------------------------------
// Inputs:
//   x - a 2D array (dim, num_pts) of points to compare
//       with the target points
//
// Outputs:
//----------------------------------------------------------------
// Computes the gradient of the data matching metric between
// sets of surface currents
//----------------------------------------------------------------
Array2D<double> tmpSurfaceCurrent::GradMatching(const ShapeObject* shape)
{
    Array2D<double> x = shape->GetPoints();
    int dim = x.GetLength();
    int nx = x.GetWidth();

    Array2D<double> g(dim,nx);
    g.FillArray(0.0f);

    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<nx; j++)
        {
            g(i,j) = 2*(x(i,j)-this->_points(i,j));
        }
    }

    return g;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// OVERLOADED OPERATORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// operator =
//----------------------------------------------------------------
// Inputs:
//   shape - tmpSurfaceCurrent object to copy
//
// Outputs:
//   return - copied tmpLandmarks object
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
tmpSurfaceCurrent& tmpSurfaceCurrent::operator = (const tmpSurfaceCurrent& shape)
{
    if (this != &shape)
    {
        this->_numPoints = shape._numPoints;
        this->_points = shape._points;
        this->_numEdges = shape._numEdges;
        this->_edges = shape._edges;
        this->_sigmaW = shape._sigmaW;
        this->_timept = shape._timept;
        this->_timeIndex = shape._timeIndex;
    }

    return *this;
}
