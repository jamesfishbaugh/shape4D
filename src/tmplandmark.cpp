//----------------------------------------------------------------
// Landmarks:  Shape data represented as landmarks
//----------------------------------------------------------------

#include "tmplandmark.h"

#include <stdio.h>			// For printf

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS/DESTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// tmpLandmarks
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
tmpLandmarks::tmpLandmarks()
{
}

//----------------------------------------------------------------
// tmpLandmarks
//----------------------------------------------------------------
// Inputs:
//   shape - tmpLandmarks object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
tmpLandmarks::tmpLandmarks(const tmpLandmarks& shape)
{
    this->_numPoints = shape._numPoints;
    this->_points = shape._points;
    this->_timept = shape._timept;
    this->_timeIndex = shape._timeIndex;
}

//----------------------------------------------------------------
// ~tmpLandmarks
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
tmpLandmarks::~tmpLandmarks()
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
// Computes the data matching metric between sets of landmarks
//----------------------------------------------------------------
double tmpLandmarks::Matching(const ShapeObject* shape)
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
//   x - a 2D array (dim, num_pts, time) of points to compare
//       with the target points
//
// Outputs:
//----------------------------------------------------------------
// Computes the gradient of the data matching metric between
// sets of landmarks
//----------------------------------------------------------------
Array2D<double> tmpLandmarks::GradMatching(const ShapeObject* shape)
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
//   shape - tmpLandmarks object to copy
//
// Outputs:
//   return - copied tmpLandmarks object
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
tmpLandmarks& tmpLandmarks::operator = (const tmpLandmarks& shape)
{
    if (this != &shape)
    {
        this->_numPoints = shape._numPoints;
        this->_points = shape._points;
        this->_timept = shape._timept;
        this->_timeIndex = shape._timeIndex;
    }

    return *this;
}
