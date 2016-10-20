//----------------------------------------------------------------
// EXLandmarks:  Target data represented as landmarks
//----------------------------------------------------------------

#include "landmarks.h"

#include <stdio.h>			// For printf

#include "gridoptimize.h"		// To make fast kernel computations

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS/DESTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXLandmarks
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
Landmarks::Landmarks()
{
	this->_ny = 0;
	this->_nvy = 0;
	this->_nvx = 0;
	this->_sigmaW = 0;
	this->_timept = 0.0;
	this->_timeIndex = -1;
	this->_weight = 0.0;
	this->_grid = 0;
}

/*
//----------------------------------------------------------------
// EXLandmarks
//----------------------------------------------------------------
// Inputs: 
//   y - a 2D array (dim, num_pts) of target points
//   vx - a 2D array (dim, num_pts) of source trianges which
//        correspond to this target
//   timept - the time corresponding to this target in whatever
//            time units used in the experiment
//   timeIndex - the index into the time discretization array
//
// Outputs:
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
EXLandmarks::EXLandmarks(const EX2DArray<double>& y, const EX2DArray<int>& vx, double timept, int timeIndex)
{
	this->_ny = y.GetWidth();
	this->_y = y;
	this->_nvx = vx.GetWidth();
	this->_vx = vx;
	this->_timept = timept;
	this->_timeIndex = timeIndex;
	this->_weight = 1.0;

}
*/

//----------------------------------------------------------------
// EXLandmarks
//----------------------------------------------------------------
// Inputs: 
//   y - a 2D array (dim, num_pts) of target points
//   vx - a 2D array (dim, num_pts) of target trianges
//   timept - the time corresponding to this target in whatever
//            time units used in the experiment
//   timeIndex - the index into the time discretization array
//   weight - the importance of this target in the regression
//
// Outputs:
//----------------------------------------------------------------
// Init constructor - recommended constructor
//----------------------------------------------------------------
Landmarks::Landmarks(const Array2D<double>& y, const Array2D<int>& vy, const Array2D<int>& vx, double sigmaW, double timept, int timeIndex, double weight)
{
	this->_ny = y.GetWidth();
	this->_y = y;
	this->_nvy = vy.GetWidth();
	this->_vy = vy;
	this->_nvx = vx.GetWidth();
	this->_vx = vx;
	this->_sigmaW = sigmaW;
	this->_timept = timept;
	this->_timeIndex = timeIndex;
	this->_weight = weight;
	
    this->_grid = new Grid();
	
	//EX2DArray<double> tempCenters(3,this->_nvy);
	//EX2DArray<double> tempNormals(3,this->_nvy);
	//this->ComputeCentersAndNormals(y, vy, tempCenters, tempNormals);
	
	//this->_centers = tempCenters;
	//this->_normals = tempNormals;
}

//----------------------------------------------------------------
// EXLandmarks
//----------------------------------------------------------------
// Inputs: 
//   shape - EXLandmarks object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
Landmarks::Landmarks(const Landmarks& shape)
{
	this->_ny = shape._ny;
	this->_y = shape._y;
	this->_nvy = shape._nvy;
	this->_vy = shape._vy;
	this->_nvx = shape._nvx;
	this->_vx = shape._vx;
	this->_sigmaW = shape._sigmaW;
	this->_timept = shape._timept;
	this->_timeIndex = shape._timeIndex;
	//this->_centers = shape._centers;
	//this->_normals = shape._normals;
	this->_grid = shape._grid;
}

//----------------------------------------------------------------
// ~EXLandmarks
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
Landmarks::~Landmarks()
{
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// INTERFACE TO COMPUTE DATA MATCHING METRIC AND GRADIENT
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// Matching
//----------------------------------------------------------------
// Inputs: 
//   x - a 3D array (dim, num_pts, time) of points to compare 
//       with the target points
//   t - the time index into x to use
//
// Outputs:
//----------------------------------------------------------------
// Computes the data matching metric between sets of landmarks
//----------------------------------------------------------------
double Landmarks::Matching(const Array3D<double>& x, int t)
{
	int dim = x.GetLength();
	int nx = x.GetWidth();
	double sum = 0.0f;

	//printf("For timept %d\n",t);

	for (int i=0; i<dim; i++)
	{
		for (int j=0; j<nx; j++)
		{
			//int curind = this->_vx(i,j);
			double curpt = x(i,j,t);
			double diff = (curpt - this->_y(i,j));
			
			sum += (diff*diff);
		}
	}

	//printf("   sum =  %f\n", sum);
	
	return sum;
}

//----------------------------------------------------------------
// GradMatching
//----------------------------------------------------------------
// Inputs: 
//   x - a 3D array (dim, num_pts, time) of points to compare 
//       with the target points
//   t - the time index into x to use
//
// Outputs:
//----------------------------------------------------------------
// Computes the gradient of the data matching metric between 
// sets of landmarks
//----------------------------------------------------------------
Array2D<double> Landmarks::GradMatching(const Array3D<double>& x, int t)
{
	int dim = x.GetLength();
	int nx = x.GetWidth();
	
    Array2D<double> g(dim,nx);
	g.FillArray(0.0f);
	
	for (int i=0; i<dim; i++)
	{
		for (int j=0; j<nx; j++)
		{
			//int curind = this->_vx(i,j);
			g(i,j) = 2*(x(i,j,t)-this->_y(i,j));
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
//   shape - EXLandmarks object to copy
//
// Outputs:
//   return - copied EXLandmarks object 
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
Landmarks& Landmarks::operator = (const Landmarks& shape)
{
if (this != &shape)
	{
		this->_ny = shape._ny;
		this->_y = shape._y;
		this->_nvy = shape._nvy;
		this->_vy = shape._vy;
		this->_nvx = shape._nvx;
		this->_vx = shape._vx;
		this->_sigmaW = shape._sigmaW;
		this->_timept = shape._timept;
		this->_timeIndex = shape._timeIndex;
		//this->_centers = shape._centers;
		//this->_normals = shape._normals;
		this->_grid = shape._grid;
	}
	
	return *this;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeCentersAndNormals
//----------------------------------------------------------------
// Inputs: 
//   pts - 2D array (dim, num_pts) of data points
//   tris - 2D array (dim, num_pts) of data triangles
//   centers - 2D array (dim, num_pts) for face centers 
//             return value
//   normals - 2D array (dim, num_pts) for face normals 
//             return value
//
// Outputs:
//----------------------------------------------------------------
// Computes face centers and normals given points and face 
// connectivity.  Returns centers and normals through reference 
// parameters.
//----------------------------------------------------------------
void Landmarks::ComputeCentersAndNormals(const Array2D<double>& pts, const Array2D<int>& tris, Array2D<double> &centers, Array2D<double> &normals)
{
	int numTris = tris.GetWidth();
	
	// Loop over all the faces
	for (int f=0; f<numTris; f++)
	{
		double v[9];
		for (int k=0; k<3; k++)
		{
			for (int j=0; j<3; j++)
			{
				v[k+3*j] = pts(k,tris(j,f)); 
			}
		}
		
		// Calculate the center of each face c = (v1+v2+v3)/3;
		centers(0,f) = ((double)(v[0]+v[3]+v[6]))/3.0f;
		centers(1,f) = ((double)(v[1]+v[4]+v[7]))/3.0f;
		centers(2,f) = ((double)(v[2]+v[5]+v[8]))/3.0f;
		
		// Calcuate the normal of each face N = [(v2-v1)a(v3-v1)]/2;
		normals(0,f) = ((double)(((v[4]-v[1])*(v[8]-v[2])-(v[5]-v[2])*(v[7]-v[1]))))/2.0f;
		normals(1,f) = ((double)(((v[5]-v[2])*(v[6]-v[0])-(v[3]-v[0])*(v[8]-v[2]))))/2.0f;
		normals(2,f) = ((double)(((v[3]-v[0])*(v[7]-v[1])-(v[4]-v[1])*(v[6]-v[0]))))/2.0f;
	}
}
